subroutine trans_ncdfinput(Nmol, Iconf)
  !
  ! Rename variables from read_nc NetCDF file reading routines
  !
  use configuration, only : org, cell_in=>cell, r_in=>r, v_in=>v,&
       & tstep=>scale, ity_in=>ity, nstep_in=>step, natoms
  use comun, only : Vel, r, cell, sidel, side, volumen, itype, bscat, tunit&
       &, tuniti,  ntype, masa, bsc, mat, nstep, natms, ndim,&
       & vector_product, UnitR, keytrj, nsp
  implicit none
  integer, intent(in) :: Iconf, Nmol
  Integer :: i, j, it
  Integer, dimension(:), allocatable :: counter, nct
  if (natoms .ne. Nmol) then
     stop(' *** Error : no. of atoms in input file and trajectory file differ !!!!')
  end if
  allocate(nct(nsp),counter(nsp))
  nstep = nstep_in(1)
  if (Iconf == 1) then
     !Define inital vars and checks
     natms = natoms
     cell(:,:) = 0
  endif
  sidel(:) = cell_in(:,1)
  side = Minval(sidel(1:ndim))
  volumen = 1.0
  do i=1,ndim
     cell(i,i) = cell_in(i,1)
     volumen = volumen*cell(i,i)
  end do
  if (allocated(v_in)) then
   keytrj = 1
  else
   keytrj = 0
  endif
   !
  ! Remap coordinates so as to have first ntype(1) particles in contiguous positions, followed  
  ! by ntype(2) particles and so on .. and transform from netcdf format to local vars
  ! 
  ntype(:) = 0
  do i=1,Nmol
   ntype(ity_in(i,1)) = ntype(ity_in(i,1))+1
  end do 
  nct(1)=0
  do i=2,nsp
   nct(i) = nct(i-1)+ntype(i-1)
  enddo
  counter(:) = nct(:)
  do i=1, Nmol
   it = ity_in(i,1)
   counter(it) = counter(it)+1
   j= counter(it)
   !
   ! NOTE: when ndim=2, z component of r_in,v_in is discarded 
   !
   if (keytrj>0) vel(1:ndim,j) = v_in(1:ndim,i,1)*tunit/tuniti
     r(1:ndim,j) = r_in(1:ndim,i,1)
     r(1:ndim,j) = r(1:ndim,j) -sidel(1:ndim)*int(r(1:ndim,j)&
          &/sidel(1:ndim))
     itype(j) = it
  enddo 
  do i=1,nsp
     masa(nct(i)+1:nct(i)+ntype(i)) = mat(i)
     bscat(nct(i)+1:nct(i)+ntype(i)) = bscat(i)
  end do
  
end subroutine trans_ncdfinput

subroutine reset_struct
  use dbscan, only : visited, core, neighbors, border
  use gpcodes, only :  gvisited, gcore, gneighbors, gborder
  implicit none
  gvisited(:) = .false.
  gcore(:) = .false.
  gborder(:) = .false.
  gneighbors(:) = 0
  visited(:) = .false.
  core(:) = .false.
  border(:) = .false.
  neighbors(:) = 0
end subroutine reset_struct

Subroutine InitVars(ndim,nsp,Nmol)
  Use dev_def
  Use clusters, only : sizedist, rcl
  Use Comun, only : bsc, atoms, mat, ntype, lty, itype, ipos, bscat,&
       & masa, vel, force, r, list, done
  Implicit none
  integer, intent(In) :: ndim, nsp, Nmol
  Allocate(bsc(nsp),atoms(nsp),mat(nsp),ntype(nsp),lty(nsp),itype(Nmol),ipos(nsp,Nmol)&
       &,bscat(Nmol),masa(Nmol))
  if (rcl>0) Allocate(sizedist(Nmol))
  Allocate(vel(ndim,Nmol))
  Allocate(force(ndim,Nmol))
  Allocate(r(ndim,Nmol))
  !
  ! Device code: allocation and data transfer
  !
  Allocate(r_d(ndim,Nmol),itype_d(Nmol),ntype_d(nsp),ipos_d(nsp,Nmol))
  Allocate( list(Nmol), done(Nmol))
  allocate(gdone(Nmol))

end Subroutine InitVars

Subroutine InitStorage(Nmol,nsp)
  Use dev_def
  Use cudafor
  Use Comun
  Use gpcodes, only : tBlock, grid
  Implicit none
  integer, intent(IN) :: Nmol, nsp
  integer :: i,j,k, ind, indm
  real :: rr
  if (idir > 0) then
     width =  nint((pwallp-pwall)/(2*deltar))
     allocate(densprof(0:width,nsp),densprofi(0:width,nsp))
     allocate(densprof_d(0:width,nsp))
     densprof(:,:) = 0.0d0
     densprof_d(:,:) = 0.0d0
  endif
  dq = 2*pi/sidel(1)
  nqmax = nint(qmax/dq)
  lmaxx = nint(qmin/(2*pi/sidel(1)))
  lmaxy = nint(qmin/(2*pi/sidel(2)))
  if (ndim ==3) lmaxz = nint(qmin/(2*pi/sidel(3)))
  lsmax = nint(minval(sidel(1:ndim))/2.0/deltar)+1
  ! Number of threads for S(Q) calculation (64)
  if (ndim == 3) then
     grid = dim3(4,4,4)
     if (mod(lmaxx,4).ne.0) lmaxx = lmaxx -4 + mod(lmaxx,4)
     if (mod(lmaxy,4).ne.0) lmaxy = lmaxy -4 + mod(lmaxy,4)
     if (mod(lmaxz,4).ne.0) lmaxz = lmaxz -4 + mod(lmaxz,4)
     tBlock = dim3(lmaxx/grid%x+1,lmaxy/grid%y+1,lmaxz/grid%z+1)
     nqmin = min(lmaxx,lmaxy,lmaxz)
  else
     grid = dim3(16,16,1)
     if (mod(lmaxx,16).ne.0) lmaxx = lmaxx -16 + mod(lmaxx,16)
     if (mod(lmaxy,16).ne.0) lmaxy = lmaxy -16 + mod(lmaxy,16)
     tBlock = dim3(lmaxx/grid%x+1,lmaxy/grid%y+1,1)
     nqmin = min(lmaxx,lmaxy)
  endif
  fk(:) = 2*pi/sidel(:)
  if (nqmin == 0) then
     nqmin = 1
     lmaxx = 1
     lmaxy = 1
     lmaxz = 1
  endif
  qmin = lmaxx*fk(1)
  Allocate(sqf(nqmax),sqfcl(nqmax),sqfp(nqmax,nsp),nq(nqmax))
  Allocate(histomix(lsmax,nsp,nsp),histomixi(lsmax,nsp,nsp)&
       &,gcluster(lsmax),gclustav(lsmax),gclcl(lsmax),rhoclus(0:lsmax),&
       & rhoclusav(0:lsmax),sqcl(nqmax)) 
  Allocate (histomix_d(lsmax,nsp,nsp),sqf_d(nqmax),sqfcl_d(nqmax),sqfp_d(nqmax&
       &,nsp))
!  print *, ' bien2'
  nq_d(:) = nq(:)
  gclustav(:) = 0.0
  rhoclus(:) = 0.0
  rhoclusav(:) = 0.0
  gclcl(:) = 0
  histomix(:,:,:) = 0
  histomix_d(:,:,:) = 0
  sqf_d(:) = 0.0
  sqfcl_d(:) = 0.0
  sqfp_d(:,:) = 0
  sqf(:) = 0.0d0
  sqfcl(:) = 0.0d0
  sqfp(:,:) = 0
  nq(:) = 0
  if (nqmin > 10) then
     indm = 0
     if (ndim==3) then
        do i=0,lmaxx-1
           do j=0,lmaxy-1
              do k=0, lmaxz-1
                 if (i+j+k.ne.0) then
                    rr = sqrt((i*fk(1))**2+(j*fk(2))**2+(k*fk(3))**2)
                    if (rr .lt. qmin) then
                       ind = nint(rr/dq)
                       if (ind > indm) indm=ind
                       if (ind <=nqmax) nq(ind)=nq(ind)+1
                    endif
                 endif
              end do
           end do
        end do
        write(*,"(' *** Using ',i4,'x',i4,'x',i4,' vectors **'/)")lmaxx,lmaxy,lmaxz
        nq(indm+1:nqmax) = 3
     else
        do i=0,lmaxx-1
           do j=0,lmaxy-1
              if (i+j.ne.0) then
                 rr = sqrt((i*fk(1))**2+(j*fk(2))**2)
                 if (rr .lt. qmin) then
                    ind = nint(rr/dq)
                    if (ind > indm) indm=ind
                    if (ind <=nqmax) nq(ind)=nq(ind)+1
                 endif
              endif
           end do
        end do
        nq(indm+1:nqmax) = 2
        write(*,"(' *** Using ',i4,'x',i4,' vectors **'/)")lmaxx,lmaxy
     end if
     qmin2 = qmin**2
     write(*,"(/' *** Averaging over all q-vectors for q <',f10.6,' A^-1 **')") qmin
  else
     if (ndim==3) then
         nq(1:nqmax) = 3
     else
      nq(1:nqmax) = 2
     endif
     nqmin = 1
     write(*,"(/' *** Averaging over (100), (010) and (001)  q-vectors  **'/)") 
  endif

  ndr = nint(sidel(1)/2.0/deltar)
  lty(:) = 0
  do i=1, Nmol
     lty(itype(i))=lty(itype(i))+1
     ipos(itype(i),lty(itype(i))) = i
  enddo
  !
  ! Device data transfer
  !
  ntype_d(:) = ntype(:)
  ipos_d(:,:) = ipos(:,:)
  itype_d(1:Nmol) = itype(1:Nmol)
end Subroutine InitStorage

Subroutine ReadCfg(iunit,io,iconf,istart)
  Use Comun
  Implicit None
  Integer, Intent(IN) :: iunit,iconf,istart
  Integer, Intent(OUT) :: io
  Integer :: i, l, j, iatm, idum1, id, ityp, imt
  Real(kind=4) :: dumy,  rtemp(3), veltemp(3), forcetemp(3)
  logical, save :: first=.true.
  !  real(kind=4), external :: surface
  io = 0
  !
  ! LAMMPS trajectory files are disordered regarding atom id's. This must be taken care of.
  ntype(:) = 0
  Read(iunit,'(1x)',iostat=io)
  Read(iunit,*)nstep
  Read(iunit,'(1x)')
  Read(iunit,*)natms
  Read(iunit,'(1x)')
  cell(:,:) = 0.0
  Do i=1, ndim
     Read(iunit,*)rlow(i),rup(i)
     cell(i,i) = rup(i)-rlow(i)
  End Do
  if (ndim==3) then
     Read(iunit,'(1x)')
  elseif (ndim==2) then
     Read(iunit,'(/1x)')
  endif

  Do i=1, natms
     If (keytrj.Eq.0) Then
        Read(iunit,*) id, ityp,  imt, Rtemp(1:3)
     Else If (keytrj.Eq.1) Then
        Read(iunit,*) id, ityp,  imt, Rtemp(1:3), Veltemp(1:3)
     Else If (keytrj.Eq.2) Then
        Read(iunit,*) id, ityp,  imt, Rtemp(1:3), Veltemp(1:3), forcetemp(1:3)
     Endif
     if (ndim == 3) then
        R(1:3,id) = Rtemp(1:3)
        Vel(1:3,id) = Veltemp(1:3)
        Force(1:3,id) = forcetemp(1:3)
     else
        R(1:2,id) = Rtemp(1:2)
        Vel(1:2,id) = Veltemp(1:2)
        Force(1:2,id) = forcetemp(1:2)
     endif
     !
     ! Rescale cell origin to (0,0,0)
     R(:,id) = R(:,id)-rlow(:)
     !
     ! Wrap coordinates inside box
     !
     where (R(:,id)>= rup(:)) R(:,id) = R(:,id) - (rup(:)-rlow(:))
     !
     ! defaul lmp trajectory file uses box units (better use custom dump)
     !
     ! forall (j=1:ndim) R(i,j) = cell(j,j)*R(i,j)
     itype(id) = ityp
     masa(id) = mat(itype(id))
     bscat(id) = bsc(itype(id))

     ntype(itype(id)) = ntype(itype(id))+1
  End Do

  !
  ! Orthogonal cells
  !
  Forall (i=1:ndim) sidel(i) = cell(i,i)
  side = Minval(sidel(1:ndim))
  !
  !
  !
  If (ndim == 3) Then
     volumen = Abs(Dot_product(cell(1:ndim,1),vector_product(cell(1:ndim,2),cell(1:ndim,3))))
  Else
     volumen = surface(cell(1:ndim,1),cell(1:ndim,2))
  End If
  if (.not. UnitR) then
     Vel(:,:) =  Vel(:,:)*tunit/tuniti
     Force(:,:) = Force(:,:)*(tunit/tuniti)**2
  Endif
  !
  !
  ! count atoms of zero mass (do not contribute to temperature calculations)
  !
  nmzero = 0
  do i=1,natms
     if (masa(i)< 1.d-6) nmzero=nmzero+1
  Enddo
  if (first .and. nmzero > 0) Then
     Print *, ' ** a total of ', nmzero,' atoms of zero mass detected'
     first = .false.
  Endif
End Subroutine ReadCfg

Subroutine Init(myid)
  Use comun
  Use Fourier
  Implicit None
  Integer, intent(IN) :: myid
  Character :: form*2
  form = "i1"
  If (myid > 9) form="i2"
  Write(fname72,'("sqsamp",'//form//',".dat")')myid
  Write(fname99,'("gmixsim",'//form//',".dat")')myid
End Subroutine Init
