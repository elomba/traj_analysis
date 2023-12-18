!
! This program analyzes a LAMMPS generated trajectory in NETCdf format (alternative .lammpstrj 
! in the format "id type id x y z vx vy vz") and computes the pair distribtion
! functions, structure factor, and if a  rcl>0 is defined to analyze
! clusters, the cluster size distribution, and for seggregated
! clusters its pdf and s(Q). When idir > 0 idir defines a non-periodic spatial dimension. Across this 
! dimension a density profile analysis is performed. 
! In first input line, dim defines the dimensionality of the problem. 
! IMPORTANT !!!: When dim=2 it is ASSUMED that the 
! the z-direction is ommited, even it is still read in the trajectory netcdf file. 
!
!
! E. Lomba, Madrid, Dec. 2022

program trj_analysis
  use myncdf
  !
  ! Rename variables from read_nc NetCDF file reading routines
  !
  !use configuration, only : org, cell_in=>cell, r_in=>r, v_in=>v,&
  !     & tstep=>scale, ity_in=>ity, nstep_in=>step, natoms
  !
  ! Module for NetCDF standard f90 library
  !
  use netcdf
  Use Comun
  Use cudafor
  Use dev_def
  Use gpcodes
  use gpucells
!  Use cpucodes
  Use dbscan
  use clusters
  use thrust
  use cells
  use linkcell
  
  implicit none

  integer, parameter :: nclmax=1000
  Integer :: Nmol, iunit=10, myid, i, j, k, l, ntm, imff, imll, &
       & Iconf, isel, io,  jsel, istart, ind, &
       &  istat, nbcuda,nthr, ll, im1, ip1, ncl5, ndist, ncid
  Real (myprec) :: ecaver, rr2, rr,&
       & temperature, Tfact, taver, cpu0, cpu1, cpu2, cpu3,cpu4, ri, deltaV, gr,&
       & grt, xfj, error, box, const, s11, s12, s22, scc,&
       & x1, x2, valor, sumx, sumy, sumz,  high, low, zpos
  Real (kind=8) :: ekin, suma
  logical :: cpot, ncdf=.false., found
  Character ::History_file*18, error_msg*64
  real :: time, t1, t2, t3, t4, tGPU, tsQ, trdf,tread,&
       &  dist, r2, norm

  call cpu_time(t1)
  call cpu_time(cpu0)
  istat = cudaSetDevice(0)

  istat = cudaGetDeviceProperties(prop, 0)
  write(*,"(/,'Device Name: ',a)") trim(prop%name)
  write(*,"('Compute Capability: ',i0,'.',i0)") &
       prop%major, prop%minor
  myid = 0
  open(8,file='Input')
  Read(8,*) Nconf, ndim
  Read(8,*) Nmol, UnitR, tuniti, use_cell
  Read(8,*) deltar, keytrj, tstep, rcl, jmin, dcl, minclsize
  read(8,*) qmin, qmax,  sigma
  Read(8,*) nsp
  if (nsp > nspmax) then
     stop(" Error : no. of species too large ")
  end if
  !
  ! Initialize some CUDA vars
  !
  nbcuda = Nmol/nthread
  if (nbcuda*nthread .ne. Nmol) nbcuda=nbcuda+1
  bl2 = dim3(Nmol/16,Nmol/16,1)
  thr2 = dim3(16, 16,1)
  call InitVars(ndim,nsp,Nmol)
  call init_dbscan(Nmol)
  call init_gdbscan(Nmol)
  Read(8,*) bsc(1:nsp)
  do i=1, nsp
     Read(8,*) atoms(i), mat(i)
  Enddo
  Call Init(myid)
  Read(8,'(A)') History_file
  Print *, ' Opening history file ',History_file
  if (index(History_file,".nc") > 0) then
     print *, "** Trajectory in NETCDF file format"
     call check(nf90_open(path = History_file, mode=NF90_WRITE, ncid = ncid),io)
     if (io /= nf90_noerr) then
        stop(' *** Error opening Netcdf file '//History_file//' ... exiting !')
     end if
     io = 0
     ncdf = .true.
  else
     Print *, ' Reading history ..'
     Open (iunit, file=History_file)
  endif
  read(8,*) idir
  if (idir > 0) then
     if (idir > ndim) stop(" ** Error : wrong non periodic dimension")
     write(*, "(' ** reading wall positions in non periodic dimension ',a1)")npdim(idir)
     pdir(idir) = .false.
     read(8,*) pwall, pwallp
  end if
     
  if ( Nmol > nmolmax) Then
     write(error_msg,'(" **** Error: redimension nmolmax=",i7," as >= ",i7)')nmolmax,Nmol
     stop(error_msg)
  endif
  Write(*,'(A,I6)') "  Number of atoms : ",Nmol
  ecaver = 0.0
  ekclaver = 0
  ekinclsav = 0
  ! Initialize timers
  ! Start GPU timing
  istat = cudaEventCreate(startEvent)
  istat = cudaEventCreate(stopEvent)
  istat = cudaEventRecord(startEvent,0)
  sizedist(:) = 0
  if (rcl>0) then
     open(188,file='centers.lammpstrj')
     open(200,file='clusevol.dat')
  endif
  tsQ = 0.
  trdf = 0.
  tread = 0.0
  Do Iconf = 1, Nconf
     call cpu_time(cpu2)
     if (ncdf) then
        call read_nc_cfg(ncid,Iconf,io,unit=55)
        ! Transform variables read from ncdf file to appropriate
        ! format and units for processing
        call trans_ncdfinput(Nmol,Iconf)
     else
        Call  ReadCfg(iunit,io,Iconf,istart)
     endif
    
     call cpu_time(cpu3)
     tread = tread + cpu3-cpu2
     if (rcl > 0) then
        if (use_cell) then
           if (Iconf==1) call Init_Cell(rcl)
           call build_cells
        end if
     endif
     Call reset_struct
     if (Iconf .eq. 1) then
        call InitStorage(Nmol,nsp)
        if (rcl>0) call initclus(Iconf)
     end if
     side2 = (side/2)**2
     fk(:) = 2*pi/sidel(:)
     sidel_d(:) = sidel(:)
     If (natms.Ne.Nmol) Then
        write(error_msg,"(' Error: input parameter Nmol=',i7,' must be equal to ',i7)")Nmol,natms
        stop(error_msg)
     End If
     done(:) = 0
     call cpu_time(cpu2)
     !
     ! Transfer data to device
     !
     r_d(1:ndim,1:Nmol) = r(1:ndim,1:Nmol)
     fk_d(:) = fk(:)
     !
     if (rcl > 0) then
        !
        call cluster_search(Iconf,Nmol,stopEvent,startEvent,nbcuda,nthread)
     endif
     !
     bav = Sum(ntype(1:nsp)*bsc(1:nsp))/Real(natms,myprec)
     bav2 = bav**2
     b2av = ((Sum(ntype(1:nsp)*bsc(1:nsp)*bsc(1:nsp))/Real(natms,myprec)))
     If (io < 0 ) Then
        Write(*,'(/" *** Last Configuration reached; exiting main loop")')
        Nconf = Iconf-1
        Exit
     Endif
     densty= natms/volumen

     If (keytrj > 0) Then
        ekin = masa(natms)*Dot_product(vel(1:ndim,natms&
             &),vel(1:ndim,natms))
     End If
     Do i=1, natms-1
        If (keytrj > 0) Then
           ekin = ekin + masa(i)*Dot_product(vel(1:ndim,i)&
                &,vel(1:ndim,i))
        Endif
     End Do
     ekin = 0.5d0*ekin
     ecaver = ecaver + ekin
     If (keytrj > 0) Then
        Tfact = (natms-nmzero-1)*ndim
        if (UnitR) Then
           taver = 2*ecaver/(Tfact*Iconf)
           temperature = 2*ekin/(Tfact)
        Else
           taver = 2*ecaver*(aunit/tunit)**2/(Tfact*Rgas*Iconf)
           temperature = 2*ekin*(aunit/tunit)**2/(Tfact*Rgas)
        Endif
     Endif
     ! Cuda kernel calls
     !
     ! Call RDF
     !
     t2= gptime(stopEvent,startEvent)
     call RDFcomp(Nmol,Iconf,nbcuda,nthread)
     
     if (idir > 0 ) then
        !
        ! If required compute density profile along idir direction
        call  profile_comp(Nmol,Iconf,nbcuda,nthread)
     endif

     t3 = gptime(stopEvent,startEvent)
     trdf = trdf + t3-t2
     !
     ! Comput S(Q)'s
     !
     call SQcalc(Nmol)
     t4 = gptime(stopEvent,startEvent) 
     tsQ = tsQ + t4-t3
     if (rcl > 0) then
        ! Create cluster objects (reallocate memory for each
        ! configuration) 
        !
        call cluster_analysis(Iconf,Nmol)
     endif
     ! Periodic output
     If (Mod(Iconf-1,1).Eq.0) Then
        call cpu_time(cpu1)
        Write(*,"(/' ** Working on MD step no. ',i8,' cpu time=',f15.2&
             &/)") nstep, cpu1-cpu0
        cpu0 = cpu1
        write(*,"(' ** Kinetic energy=',f15.4,' Kcal/mol, average=',f15.4,'Kcal/mol')")0.00198717*ekin*(aunit/tunit)**2/Rgas ,0.00198717*ecaver*(aunit/tunit)**2/Rgas/Iconf 
        if (rcl>0) then
           write(*,"(' ** Cluster kinetic energy=',f15.4,' Kcal/mol, average=',f15.4,'Kcal/mol')")0.00198717*ekincl*(aunit/tunit)**2/Rgas, 0.00198717*ekclaver*(aunit/tunit)**2/Rgas/Iconf
           write(*,"(' ** Internal cluster kinetic energy=',f15.4,'&
                & Kcal/mol, average=',f15.4,'Kcal&
                &/mol')")0.00198717*ekincls*(aunit/tunit)**2/Rgas,&
                & 0.00198717*ekinclsav*(aunit/tunit)**2/Rgas/Iconf
        endif
        If (keytrj > 0) then
           if (rcl > 0) then
              Tfact = nint(sum(sizedist(:))/real(Iconf))*ndim
               Write(*,"(' ** Average cluster temperature =',f10.4&
                   &,' K')") 2*ekclaver*(aunit/tunit)**2/(Tfact*Rgas*Iconf)
            endif
           if (UnitR) Then
              Write(*,"(' ** Temperature*=',f10.4&
                   &,' average=',f10.4,' density*=',f10.6&
                   &)")temperature,taver,densty
           Else
              Write(*,"(' ** Temperature=',f10.4&
                   &,' K average=',f10.4,'K density=',f10.6' 1&
                   &/A^3')")temperature,taver,densty
           endif
        else
           Write(*,"(' Density*=',f10.6&
                &)")natms/volumen
        endif
        if (rcl>0) then
           write(*,"(' ** Average cluster radius',f8.3,' average &
             &cluster density ',f10.7)")avradio/iconf, averdens/iconf
            write(*,"(' ** Internal cluster density  ',f10.7)")&
                 & sum(densav(:))/Nu_clus
            write(*,"(' ** No. of clusters for this configuration :',i5)") Nu_clus
            print *, " ··Time for graph construction", tgraph/iconf
            print *, " ··Thrust time ",tthrus/iconf
            print *, " ··Time adjacency list construction =", tadj/iconf
            print *, " ··Time for BFS cluster search =", tbfs/iconf
        endif
        print *, " ··Time for rdf ", trdf/iconf
        print *, " ··Time por S(Q) ", tsQ/iconf
        print *, " ··Time config in/out  ", tread/iconf
     Endif
  end Do
  ! Normalize density profiles computed along the non-periodic dimension
   if (idir>0) then
     found = .false.
     i=0
     do while (.not. found)
        if (maxval(densprof(i,1:nsp)) < 1.0e-6) then 
           low = (2*i+1)*deltar+pwall
           i=i+1
        else
           found=.true.
        endif
     end do
     found = .false.
     i=width
     do while (.not. found) 
        if (maxval(densprof(i,1:nsp)) < 1.0e-6) then 
           high = (2*i+1)*deltar+pwall
           i=i-1
        else
           found=.true.
        endif
     end do
     high = (pwallp+high)/2
     low = (pwall+low)/2
     write(*,"('   Highest accesible constrained coordinate for fluid',f15.7,/'   Lowest accessible constrained coordinate for fluid',f15.7)")high,low
     volumen_r = (high-low)*product(sidel,pdir)
     rdenst = natms/volumen_r
     write(*,'(/" ** Total renormalized density =",f10.6)') rdenst
     do i = 1, nsp
        write(*,'(/"     - Renormalized density per species",i1," =",f10.6)') i, ntype(i)/volumen_r
     end do
     open(222,file='densprof.dat')
     !
     ! One dimensional rho(x) normalized such that rho_i(L/2) approximates the bulk average rho_i 
     !
     write(222,'("#       r ",7x,8(8x,"rho(",i1,")":))')((i),i=1,nsp)
     do i=0, width
        zpos = (2*i+1)*deltar+pwall
        write(222,'(10f15.7)')zpos,(high-low)*densprof(i,1:nsp)/(Nconf*2*deltar*volumen_r)
     enddo
     close(222)
  endif

  !
  ! Recover device values for S(Q) and particle histograms
  !
  sqf(:) = sqf_d(:)
  if (rcl>0) sqfcl(:) = sqfcl_d(:)
  sqfp(1:nqmax,1:nsp) = sqfp_d(1:nqmax,1:nsp)
  write(*,"(/60('-')/' ** End: Total GPU time ',f15.7,' secs **')") gptime(stopEvent,startEvent)
  !
  ! destroy timing events
  istat = cudaEventDestroy(startEvent)
  istat = cudaEventDestroy(stopEvent)
  Open(99,file=fname99)
  !
  ! Printout S(Q)'s 
  ! 

  open(100,file='sq.dat')
  open(110,file='sqmix.dat')
  x1 = (real(ntype(1))/real(Nmol))
  x2 = (real(ntype(2))/real(Nmol))
  do i=1, nqmax
     if (nsp == 2) then
        s11 = x1*sqfp(i,1)/(ntype(1)*Nconf*real(nq(i)))
        s22 = x2*sqfp(i,2)/(ntype(2)*Nconf*real(nq(i)))
        s12 = 0.5*(sqf(i)/(Nmol*Nconf*real(nq(i))) - s11 - s22)
        scc = x2**2*s11+x1**2*s22-2*x1*x2*s12
        write(100,'(6f15.7,i12)')i*dq,sqf(i)/(Nmol*Nconf*real(nq(i)))&
             &,scc,s11,s22,s12,nq(i)
     else
        write(100,'(2f15.7,i12)')i*dq,sqf(i)/(Nmol*Nconf*real(nq(i))),nq(i)
     endif
     write(110,'(7f15.7)')i*dq,(sqfp(i,j)/(ntype(j)*Nconf&
          &*real(nq(i))),j=1,nsp)
  end do
  do i = 1, nsp
     write(*,'(" ** ",i6," atoms of type ",i2)')ntype(i), i
  end do
  !
  ! Printout partial pair distribution functions
  !
  call printrdf(rcl,lsmax)
  if (rcl > 0) then
     call print_clusinfo(nqmin,Nmol)
  end if
  call cpu_time(t2)
  write(*,"(' ** CPU time ',f15.7,' secs **')")t2-t1

end program trj_analysis

