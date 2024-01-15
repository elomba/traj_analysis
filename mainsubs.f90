!
! Basic subroutines to control calculations
!

subroutine initclus(Iconf)
   !
   ! Init cluster arrays
   !
  use clusters
  Use comun, only : ndr
  implicit none
  integer, intent(IN) :: Iconf
  drho =0.05
  ndrho = nint(1.4/drho)
  allocate(densclus(0:ndrho),radii(0:ndr))
  densclus(:) = 0
  radii(:) = 0
end subroutine initclus
subroutine cluster_search(iconf,nmol,stopEvent,startEvent,nbcuda,nthread)
   !
   ! Main cluster search call
   !
  use cudafor
  use comun
  use clusters
  use gpucells
  use gpcodes
  use dev_def
  use dbscan
  use thrust
  use linkcell, only : use_cell
  implicit none
  type(cudaEvent), intent(INOUT) :: startEvent, stopEvent
  integer, intent(IN) :: iconf, nmol, nbcuda, nthread
  integer :: cluster_id, i, j
  real :: t2, t3, t4
  Character :: error_msg*64
  t2 = gptime(stopEvent,startEvent)
  !
  ! Build neighbor tables for graph constructtion
  !
  if (use_cell) then
     call gpu_graph_cell<<<nbcuda,nthread>>>(r_d,Nmol,ndim,rcl**2,sidel_d,jmin)
  else
     call gpu_graf<<<nbcuda,nthread>>>(r_d,Nmol,ndim,rcl**2,sidel_d,jmin)
  endif
  goffset(:) = gneighbors(:)
  neighbors(:) = gneighbors(:)
  core(:) = gcore(:)
  t3 = gptime(stopEvent,startEvent) 
  tgraph = tgraph + t3-t2
  call thrustscan(goffset,Nmol,1)
  offset(:) = goffset(:)
  t2 =  gptime(stopEvent,startEvent) 
  tthrus = tthrus + t2-t3
  !
  ! Determine adjacency between nodes
  !
  if (use_cell) then
     call gpu_adj_cell<<<nbcuda,nthread>>>(r_d,Nmol,ndim,rcl**2,sidel_d)
  else
     call gpu_adj<<<nbcuda,nthread>>>(r_d,Nmol,ndim,rcl**2,sidel_d)
  endif
  adjacency(:) = gadjacency(:)
  gdone(:) = 0
  done(:) = 0
  cluster_id = 0
  visited(:) = .false.
  t4 = gptime(stopEvent,startEvent) 
  tadj = tadj + t4-t2
  !
  !
  ! Recursive BFS search of connected clusters
  do i = 1, nmol
     if (core(i) .and. done(i)==0) then
        gvisited(:) = .false.
        cluster_id = cluster_id+1
        call  breadfirst(i,nbcuda, nthread, Nmol,jmin,cluster_id, border, visited, done)
     endif
  end do
  t3 = gptime(stopEvent,startEvent) 
  !
  ! Timinh resuls
  tbfs = tbfs + t3-t4
!  print *, " ** Total number of clusters found : ", cluster_id, " Conf no. ", iconf
  if (maxval(neighbors(:))>50) then
     write(error_msg,'(" Warning:redimension next to at least ", i5)') maxval(neighbors(:))
     stop(error_msg)
  end if
  Nu_clus = cluster_id
end subroutine cluster_search

subroutine cluster_analysis(Iconf, Nmol)
   !
   ! Compute cluster properties
   !
  use comun
  use clusters
  use gpucells
  use gpcodes
  use dev_def
  use cpucodes, only : minj, rdfcl, cldens, rdfclcl, denspcl
  implicit none
  integer, intent(IN) :: Iconf, Nmol
  integer :: i, j, k, cluster_id, icl, ip, ind, indr, maxcln, maxcl,&
       & imaxcl, clsz, nbigcl, clussize, cmin 
  Real (myprec), Dimension(3) ::  c0, rv, rv2, rvi, rvj, desp, xi,si, xss,sss, theta, clcent, vcl
  Real (myprec) :: denscl
  Real (kind=8) ::  mcl, ekcls, suma, rclus, rcluster, volcl
  if (Nu_clus < 1) then
     maxcln =1
  else
     maxcln = Nu_clus
  endif
  if ( Iconf > 1) then
     deallocate(cluster,contador)
  endif
  allocate(cluster(maxcln),contador(maxcln))
  Cluster_id = maxcln
  contador(:) = 0
  maxcl = 0
  do i = 1, cluster_id
     clsz = count(done==i)
     if (clsz > maxcl) then
        maxcl = clsz
        imaxcl = i
     endif
     cluster(i)%clsize = clsz
     allocate(cluster(i)%members(clsz))
  end do
  contador(:) = 0
  do i = 1, Nmol
     j = done(i)
     if (j> 0) then
        contador(j) = contador(j)+1
        cluster(j)%members(contador(j)) = i
     endif
  enddo
  nbigcl = 0
  do i=1, cluster_id
     if (cluster(i)%clsize >= minclsize) nbigcl=nbigcl+1
  end do
  !
  ! Write trajectory file with centers of mass of clusters with size > minclsize
  !
  write(188,"('ITEM: TIMESTEP'/I12/'ITEM: NUMBER OF ATOMS'/I12/'ITE&
       &M: BOX BOUNDS pp pp pp')")nstep, nbigcl
  write(188,"(2f15.7)")(0.0,sidel(i),i=1,ndim)
  if (ndim == 2) write(188,"('-0.5 0.5')")
  write(188,"('ITEM: ATOMS id type id x y z')")
  !        write(189,*) " iconf=", iconf
  icl = 0
  do I = 1, cluster_id
     j = cluster(i)%clsize
     ! Calculate cluster size distro
     ind = nint(real(j)/real(dcl))
     sizedist(ind) =  sizedist(ind)+1
     if (j>=jmin) then
        ! Determine center of mass with PBC
        ! https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
        xss(:) = 0
        sss(:) = 0
        vcl(:) = 0
        mcl = 0
        ekcls = 0
        do k = 1, j
           ip = cluster(i)%members(k)
           rv(1:ndim) = r(1:ndim,ip)
           theta(1:ndim) = 2*pi*rv(1:ndim)/sidel(1:ndim)
           xss(1:ndim) = xss(1:ndim)+cos(theta(1:ndim))
           sss(1:ndim) = sss(1:ndim)+sin(theta(1:ndim))
           vcl(1:ndim) = vcl(1:ndim)+masa(ip)*vel(1:ndim,ip)
           mcl = mcl+masa(ip)
        end do
        xss(1:ndim)=xss(1:ndim)/j
        sss(1:ndim)=sss(1:ndim)/j
        theta(1:ndim) = atan2(-sss(1:ndim),-xss(1:ndim))+pi
        !
        ! Restore origin for centers of mass position to the original coordinate
        !
        cluster(i)%center(1:ndim) =  sidel(1:ndim)*theta(1:ndim)/(2*pi)+rlow(1:ndim)
        !
        ! Store cluster's center of mass velocity and net mass
        !
        cluster(i)%vl(1:ndim) = vcl(1:ndim)/mcl
        cluster(i)%mass = mcl
        ekcls = 0
        do k = 1, j
           ip = cluster(i)%members(k)
           ekcls = ekcls + 0.5*masa(ip)*dot_product(vel(:,ip)&
                &-cluster(i)%vl(:),vel(:,ip)-cluster(i)%vl(:))
        enddo
        cluster(i)%ekin = ekcls
     endif
     if (j>=minclsize)then
        icl = icl+1
        if (ndim ==3) then
         write(188,'(3i10,3f15.7)')icl,nsp+2,icl,cluster(i)%center(1:ndim)
        else
         write(188,'(3i10,3f15.7)')icl,nsp+2,icl,cluster(i)%center(1:ndim),0.0
        endif 
     end if
  end do
  ekincl = 0
  ekincls = 0
  suma = 0
  !
  ! Determine single cluster kinetic energies and inter-cluster kinetic contribution
  !
  do i = 1, cluster_id
     suma = suma+cluster(i)%clsize
     if (cluster(i)%clsize >= minclsize) then
        ekincl = ekincl + 0.5*cluster(i)%mass&
             &*dot_product(cluster(i)%vl(:),cluster(i)%vl(:))
        ekincls = ekincls + cluster(i)%ekin
     endif
  enddo
  ekclaver = ekclaver + ekincl
  ekinclsav = ekinclsav + ekincls
  write(200,'(i7,i4,f15.7)')iconf, cluster_id, 1.0-(suma/natms)
  if (iconf > 1) then
     deallocate(rclxyz)
     deallocate(rclxyzd)
     deallocate(densav)
  endif
  allocate(rclxyz(ndim,cluster_id),densav(cluster_id))
  allocate(rclxyzd(ndim,cluster_id))
  gcluster(:) = 0
  rhoclus(:) = 0
  denscl = 0.0
  rclus =  0.0
  !
  ! Determine cluster radii, density ..
  !
  do i=1, cluster_id
     clussize = cluster(i)%clsize
     clcent(:) =  cluster(i)%center(1:ndim)
     rclxyz(1:ndim,i) =  clcent(1:ndim)
     rcluster = radius(cluster(i)%members(:),clussize,clcent,sidel)
     cluster(i)%radio = rcluster
     if (ndim==3) then
        volcl = (4*pi*rcluster**3/3.0)
     else
        volcl = pi*rcluster**2
     endif
     cluster(i)%cldens = real(clussize)/volcl
     if (ndim==3) then
      ind = nint(cluster(i)%cldens*sigma**3/drho)
     else
      ind = nint(cluster(i)%cldens*sigma**2/drho)
     end if
     indr = nint(cluster(i)%radio/deltar)
     if (ind<=ndrho) densclus(ind) = densclus(ind)+1
     if (indr<=ndr)radii(indr) = radii(indr)+1
     denscl = denscl + cluster(i)%cldens
     rclus = rclus +  rcluster
  end do
  rclxyzd(:,:) = rclxyz(:,:)
  !
  ! Compute cluster-cluster structure factor
  ! 
  if (ndim == 3) then
     call sqf3Dcl<<<tBlock,grid>>>(rclxyzd,cluster_id,fk_d,sqfcl_d,nqmax,qmin2,dq)
  else
     call sqf2Dcl<<<tBlock,grid>>>(rclxyzd,cluster_id,fk_d,sqfcl_d,nqmax,qmin2,dq)
  end if
  averdens = averdens + denscl/cluster_id
  avradio = avradio + rclus/cluster_id
  ! Determine average cluster density, cluster g(r) and density profiles
  do i = 1, cluster_id
     densav(i) = cldens(ndim,cluster(i)%members(:),cluster(i)%clsize,sidel&
          &,side2, cluster(i)%center(1:ndim),0.8*cluster(i)%radio)
     if (cluster(i)%clsize > minclsize) then
        cmin  = minj(ndim,cluster(i)%members,cluster(i)%clsize&
             &,cluster(i)%center,sidel,side2) 
        call rdfcl(ndim,cluster(i)%members(:),cluster(i)%clsize,sidel&
             &,side2,gcluster,lsmax,deltar,densav(i),cmin)
        call denspcl(ndim,cluster(i)%members(:),cluster(i)%clsize,sidel&
             &,side2,rhoclus(0:lsmax),lsmax,deltar,cluster(i)%center(1:ndim))
     endif
  end do
  gclustav(:) = gclustav(:)+gcluster(:)/nbigcl
  rhoclusav(:) = rhoclusav(:)+rhoclus(:)/nbigcl
  !
  ! Cluster-cluster pdf
  !
  call rdfclcl(ndim,rclxyz(1:ndim,1:cluster_id),cluster_id,sidel&
       &,side2,gclcl,lsmax,deltar,volumen)
  NTclus = suma
end subroutine cluster_analysis

subroutine SQcalc(Nmol)
  use comun
  use gpucells
  use gpcodes
  use dev_def
  implicit none
  integer, intent(IN) :: Nmol

  !
  ! Sample all 3D vectors for q<qmin and
  ! (1,0,0),(0,1,0),(0,0,1) beyond
  !
  if (nqmin>1) then
     if (ndim == 3) then
        call sqf3D<<<tBlock,grid>>>(r_d,Nmol,ndim,itype_d,ntype_d,ipos_d,fk_d,sqf_d,sqfp_d,&
             & nsp,nqmax,qmin2,dq)
        call sqfact<<<nqmax/nthread+1,nthread>>>(r_d,Nmol,ndim,itype_d,ntype_d,ipos_d,fk_d,&
             & sqf_d,sqfp_d,nsp,nqmax,nqmin+1)
     else
        call sqf2D<<<tBlock,grid>>>(r_d,Nmol,ndim,itype_d,ntype_d,ipos_d,fk_d,sqf_d,sqfp_d,&
             & nsp,nqmax,qmin2,dq)
        call sqfact2<<<nqmax/nthread+1,nthread>>>(r_d,Nmol,ndim,itype_d,ntype_d,ipos_d,fk_d,&
             & sqf_d,sqfp_d,nsp,nqmax,nqmin+1)
     endif
  else
     !
     ! for very Uniform fluids and large no. of configurations
     ! just sample   (1,0,0),(0,1,0),(0,0,1) directions
     !
     if (ndim == 3) then
        call sqfact<<<nqmax/nthread+1,nthread>>>(r_d,Nmol,ndim,itype_d,ntype_d,ipos_d,fk_d,sqf_d,sqfp_d,nsp,nqmax,nqmin)
     else
        call sqfact2<<<nqmax/nthread+1,nthread>>>(r_d,Nmol,ndim,itype_d,ntype_d,ipos_d,fk_d,sqf_d,sqfp_d,nsp,nqmax,nqmin)
     end if
  end if
   end subroutine SQcalc

subroutine printrdf(rcl, lsmax)
   !
   !   Print out pdf's
   !
  use comun
  implicit none
  real(myprec), intent(IN) :: rcl
  Real (myprec) :: gmix(nspmax,nspmax), deltav, ri, xfj
  integer, intent(in) :: lsmax
  integer :: i, j, l, k
  Open(99,file=fname99)
  if (rcl>0) then
   write(99,"('#       r',16x,'g_cl(r)        g_cl-cl(r)',5x,6('g_',2i1,'(r)',8x:))")(((j,k),k=j,nsp),j=1,nsp)
  else   
   write(99,"('#       r',16x,6('g_',2i1,'(r)',9x:))")(((j,k),k=j,nsp),j=1,nsp)
  end if
  Do i = 1, lsmax-2
     ri = i*deltar
     !
     ! Compute 3d ad 2d normalizations factors
     !
     if (ndim == 3) then
        deltaV = 4*pi*((ri+deltar/2)**3-(ri-deltar/2)**3)/3.0
     else
        deltaV = pi*((ri+deltar/2)**2-(ri-deltar/2)**2)
     endif
     !
     !
     Do j=1,nsp

        Do l=j,nsp
           xfj = real(ntype(j),kind=8)/Real(natms,kind=8)
           gmix(j,l) = (j/l+1)*volumen*histomix(i,j,l)/(deltaV*ntype(l)*ntype(j)*Nconf)
           if (idir >0) gmix(j,l) = densty*gmix(j,l)/rdenst
           
        End Do
     End Do
     if (rcl > 0) then
        Write(99,'(18f16.7)')i*deltar,&
             & gclustav(i)/(deltaV*Nconf),2*gclcl(i)/(deltaV*Nconf), (gmix(j,j:nsp),j=1,nsp)
     else
        Write(99,'(18f16.7)')i*deltar,(gmix(j,j:nsp),j=1,nsp)
     endif
  End Do
  close(99)
end subroutine printrdf

subroutine print_clusinfo(nqmin, Nmol)
  use comun 
  use clusters
  implicit none
  integer, intent(IN) :: nqmin, Nmol
  integer :: i, ndist
  real(myprec) :: avcldens, deltaV, ri, suma,  norm
  Write(*,"(' ** Average total number of particles in clusters ', f10.2)")NTclus/nconf
  Write(*,"(' ** Average total number of clusters ', I5)") nint(sum(sizedist(:))/real(nconf))
  avcldens =  sum(sizedist(:)/real(nconf))/volumen
  Write(*,"(' ** Average cluster density ', f15.9)") avcldens
  open(125,file='dens.dat')
  open(126,file='radii.dat')
  open(999,file='rhoprof.dat')
  open(1001,file='clustdistr.dat')
  write(125,"('#      rho_cl        rho_cl*        N(rho_cl)')")
  do i = 1, ndrho
     write(125,"(5f15.7)")i*drho,i*drho/sigma**3,(real(densclus(i))/drho/Nconf)
  enddo
  write(126,"('#      r                  R_cl(r)')")
  write(999,"('#      r                  rho_cl(r)')")
  write(999,'(2f15.6)') 0.0, rhoclusav(0)/(4*pi*((deltar&
  &/2)**3)/3.0*Nconf)
  do i = 1, ndr
     ri = i*deltar
     if (ndim == 3) then
        deltaV = 4*pi*((ri+deltar/2)**3-(ri-deltar/2)**3)/3.0
     else
        deltaV = pi*((ri+deltar/2)**2-(ri-deltar/2)**2)
     endif
     if (radii(i).ne.0) write(126,"(2f15.7)")ri,(real(radii(i))/sum(radii(:))/deltar)
     write(999,'(2f15.6)') ri, rhoclusav(i)/(deltaV*Nconf)
  end do

  open(102,file='sqcl.dat')
  write(102,"('#      Q            S_cl-cl(Q)')")
  do i = 1, nqmin
     write(102,'(2f15.7)')i*dq,sqfcl(i)/(Nconf*real(nq(i)))
  enddo

  !
  ! Normalize cluster size distribution a s*n(s)/Ntotal
  !
  write(1001,"('#   N                %clus.               rho_cl ')")
  suma = 0
  ndist = nint(real(Nmol)/real(dcl))
  norm = (0.5*(sizedist(1)+sizedist(ndist))+sum(sizedist(2:ndist-1)))*dcl
  do i = 1, ndist
     suma = i*sizedist(i)*dcl+suma
     if (sizedist(i) > 0) write(1001,'(4f15.7)')(i-0.5)*dcl, real(i*sizedist(i))/(real(Nconf*dcl*Nmol)),real(sizedist(i))/norm
  enddo
  close(1001)
  close(999)
  close(102)
  close(125)
  close(126)
end subroutine print_clusinfo

subroutine RDFcomp(Nmol,Iconf,nbcuda,nthread)
   !
   ! Main call to compute pdfs and move data back and forth from device
   !
  use comun
  use dev_def
  use gpcodes, only : rdf_sh, rdf2_sh
  implicit none
  integer, intent(IN) :: Nmol, Iconf, nbcuda, nthread
  integer :: istart, i, ind, k, l
  histomix_d(:,:) = 0
  if (ndim == 3) then
     call rdf_sh<<<nbcuda/8+1,8*nthread>>>(r_d,Nmol,ndim,histomix_d,nit&      
          &,nsp,nbcuda*lsmax,lsmax,itype_d,side2,sidel_d,deltar)
  else
     call rdf2_sh<<<nbcuda/16+1,16*nthread>>>(r_d,Nmol,ndim,histomix_d,nit&
          &,nsp,nbcuda*lsmax,lsmax,itype_d,side2,sidel_d,deltar)
  endif
  ! Block histograms back from device to host
  histomixi(:,:) = histomix_d(:,:)
  !
  ! Gather histograms from each block
  !
  do i = 0, nbcuda-1
     istart = i*lsmax
     ind = 1
     do k = 1, nsp
      do l = k, nsp
       histomix(1:lsmax,k,l) = histomix(1:lsmax,k,l)+real(histomixi(istart+1:istart+lsmax,ind))
       ind = ind+1
      end do
   end do 
  enddo
end subroutine RDFcomp


subroutine profile_comp(Nmol,Iconf,nbcuda,nthread)
   !
   ! Compute density profile along non periodic dimension
   !
  use comun
  use dev_def
  use gpcodes, only : dprof
  implicit none
  integer, intent(IN) :: Nmol, Iconf, nbcuda, nthread
  densprof_d(:,:)= 0
  call dprof<<<nbcuda,nthread>>>(r_d,Nmol,ndim,densprof_d,width,nsp,idir,itype_d&
       &,pwall,deltar)
  densprofi(:,:) = densprof_d(:,:)
  densprof(:,:) = densprof(:,:) + real(densprofi(:,:))
end subroutine profile_comp

subroutine normdenspr
   use comun
   implicit none
   integer :: I
   real(myprec) :: high, low, zpos
   logical :: found
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
   end subroutine normdenspr

   subroutine printSQ(Nmol)
      use comun
      implicit none
      !
      ! Printout S(Q)'s  
      !   
      real(myprec) :: x1, x2, s11, s22, s12, scc
      integer, intent(IN) :: Nmol
      integer :: i, j
  open(100,file='sq.dat')
  open(110,file='sqmix.dat')
  x1 = (real(ntype(1))/real(Nmol))
  x2 = (real(ntype(2))/real(Nmol))
  if (nsp==2) then
   write(100,"('#           Q        S_NN(Q)          S_cc(Q)       S_11(Q)        S_22(Q)           S_12(Q)         n(Q)')")
  else
   write(100,"('#           Q       S_NN(Q)          n(Q)')")
  end if
  write(110,"('#       Q',14x,6('S_',2i1,'(Q)',9x:))")((j,j),j=1,nsp)
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
  close(100)
  close(110)
  end subroutine printSQ
