module thrust

interface thrustscan
subroutine scan_int( input,N, start) bind(C,name="scan_int_wrapper")
use iso_c_binding
integer(c_int),device:: input(*)
integer(c_int),value:: N
integer(c_int),value:: start
end subroutine

subroutine scan_float( input,N, start) bind(C,name="scan_float_wrapper")
use iso_c_binding
real(c_float),device:: input(*)
integer(c_int),value:: N
real(c_float),value:: start
end subroutine

subroutine scan_double( input,N, start) bind(C,name="scan_double_wrapper")
use iso_c_binding
real(c_double),device:: input(*)
integer(c_int),value:: N
real(c_double), value:: start
end subroutine

end interface

end module thrust

module gpcodes
  use cudafor
  use comun, only: nmgrid, nspmax
  implicit none
  integer, parameter :: nthread=64
  integer, device, allocatable :: goffset(:), gneighbors(:), gadjacency(:)
  logical (kind=1), device, allocatable :: gvisited(:), gcore(:), gborder(:)
  type (dim3) :: grid, tBlock, thr2, bl2
contains
 
   attributes (device) integer function fij(i,j,nsp)
      implicit none
      integer, value :: i, j, nsp
      if (i<=j) then
        fij = (i-1)*(nsp-1)-(i-1)*(i-2)/2+j
       else
         fij = (j-1)*(nsp-1)-(j-1)*(j-2)/2+i
       endif
end function fij

   subroutine init_gdbscan(Nmol)
    integer, intent(IN) :: Nmol
    allocate(goffset(Nmol), gneighbors(Nmol), gvisited(Nmol),gcore(Nmol), gborder(Nmol), gadjacency(20*Nmol))
    gvisited(:) = .false.
    gcore(:) = .false.
    gborder(:) = .false.
    gneighbors(:) = 0
  end subroutine init_gdbscan

  attributes (global) subroutine gpu_graf(r,Nmol,dim, rcl2,sidel,jmin)
     integer, value, intent(IN) :: Nmol, dim
     integer, value :: jmin
     real, value :: rcl2
    real(kind=4), intent(IN) :: r(dim,Nmol), sidel(dim)
    real (kind=4) :: rv(3), rr2
    integer :: i, j, counter
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol) then
    gneighbors(i) = 0
       counter = 0
       do j=1, Nmol
          if (i.ne.j) then
             rv(1:dim) = r(1:dim,i)-r(1:dim,j)
             ! PBC 
             rv = rv-sidel*nint(rv/sidel)
             rr2= Dot_product(rv,rv)
             if (rr2 < rcl2) then
                counter = counter+1
             endif
          end if
       end do
       gneighbors(i) = counter
       if (counter >= jmin-1)  gcore(i)= .true.
    endif
  end subroutine gpu_graf


  attributes (global) subroutine gpu_adj(r,Nmol,dim, rcl2,sidel)
    integer, value, intent(IN) :: Nmol, dim
    real, value :: rcl2
    real(kind=4), intent(IN) :: r(dim,Nmol), sidel(dim)
    real (kind=4) :: rv(3), rr2
    integer :: i, j, counter
    integer, device :: count
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol) then
       counter = 0
       do j=1, Nmol
          if (i.ne.j) then
             rv(1:dim) = r(1:dim,i)-r(1:dim,j)
             ! PBC 
             rv = rv-sidel*nint(rv/sidel)
             rr2= Dot_product(rv,rv)
             if (rr2 < rcl2) then
                gadjacency(goffset(i)+counter) = j
                counter = counter+1
             endif
          end if
       end do
    endif
  end subroutine gpu_adj


  
  subroutine breadfirst(i,nbcuda, nthread, Nmol,jmin,cluster_id, border, visited, done)
    integer, intent(IN) :: Nmol, cluster_id, jmin, i, nbcuda, nthread
    integer, intent(INOUT) :: done(Nmol)
    logical (kind=1), intent(INOUT) :: visited(Nmol), border(Nmol)
    logical (kind=1) :: gb 
    integer :: j
    border(i) = .true.
    gborder(i) = .true.
    do while (count(border==.true.) .ne. 0)
       call gbfs<<<nbcuda,nthread>>>(Nmol, gvisited, gborder, cluster_id, jmin)
       border(:) = gborder(:)
    end do
    visited(:) = gvisited(:)
    do j = 1, nmol
       if (visited(j)) then
          done(j) = cluster_id
       endif
    enddo
  end subroutine breadfirst
  
  attributes (global) subroutine gbfs(Nmol, visited, border, cluster_id, jmin)
    integer, value, intent(IN) :: Nmol, cluster_id, jmin
    logical (kind=1), intent(INOUT) :: visited(Nmol), border(Nmol)
    integer :: i, j, k
    j = (blockidx%x-1) * blockdim%x + threadidx%x
    if (j<=Nmol) then
       if (border(j)) then
          border(j) = .false.
          visited(j) = .true.
          do i = goffset(j),goffset(j)+gneighbors(j)-1
             k =  gadjacency(i)
             if (.not.visited(k)) then
                border(k) = .true.
             endif
          enddo
       end if
    end if
  end subroutine gbfs




  attributes (global) subroutine rdf_sh(r,Nmol,dim,histomix,nit,nsp,hdim,lsmax,itype&
       &,side2,sidel,deltar)
       !
       ! Mossively parallel rdf calculation using shared memory
       !
       use comun, only : dimsh, nitmax
    integer, value, intent(IN) :: Nmol, dim, nit, lsmax, hdim, nsp
    integer, intent(IN) :: itype(Nmol)
    integer i, j, ind, ia, fact, iti, itj, istart, ij
    real ::  rr2, rr, xi, yi, zi, xd, yd, zd
    real, dimension(3) :: rv
    real, value, intent(IN) :: side2, deltar
    real, intent(IN) :: sidel(3)
    real, dimension(dim,Nmol), intent(IN) :: r
    integer, intent(INOUT) :: histomix(hdim,nit)
    integer, shared :: histomix_s(dimsh,nitmax)
    histomix_s(:,:) = 0
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol-1) then
      iti = itype(i)
       xi = r(1,i)
       yi = r(2,i)
       zi = r(3,i)
       Do j=i+1,Nmol
          xd = r(1,j)-xi
          yd = r(2,j)-yi
          zd = r(3,j)-zi
          xd = xd -sidel(1)*nint(xd/sidel(1))
          yd = yd -sidel(2)*nint(yd/sidel(2))
          zd = zd -sidel(3)*nint(zd/sidel(3))
          rr2= xd*xd+yd*yd+zd*zd
          If (rr2.Lt.side2) Then
             itj = itype(j)
             rr = __fsqrt_rn(rr2)
             ind = Nint(rr/deltar)
             !
             ! Use shared memory histogram to speed calculations
             !
             ij = fij(iti,itj,nsp)
             ia = atomicadd(histomix_s(ind,ij),1)
          !   if (iti /= itj) ia = atomicadd(histomix_s(ind,ij),1)
          endif
       Enddo
       call syncthreads()
       if (threadidx%x == 1) then
         istart = (blockidx%x-1)*lsmax
         ! 
         ! Store each block's shared memory histogram on different positions
         ! in global memory
         !
         histomix(istart+1:istart+lsmax,:) = histomix_s(1:lsmax,:)+histomix(istart+1:istart+lsmax,:)
       endif 
    end if

  end subroutine rdf_sh

  attributes (global) subroutine rdf2_sh(r,Nmol,dim,histomix,nit,nsp,hdim,lsmax,itype&
       &,side2,sidel,deltar)
       !
       ! Masively parallel rdf calculation using shared memory (2D)
       !
       use comun, only : dimsh, nitmax
    integer, value, intent(IN) :: Nmol, dim, nit, lsmax, hdim, nsp
    integer, intent(IN) :: itype(Nmol)
    integer i, j, ind, ia, fact, iti, itj, istart, ij
    real ::  rr2, rr, xi, yi, zi, xd, yd, zd
    real, dimension(3) :: rv
    real, value, intent(IN) :: side2, deltar
    real, intent(IN) :: sidel(3)
    real, dimension(dim,Nmol), intent(IN) :: r
    integer, intent(INOUT) :: histomix(hdim,nit)
    integer, shared :: histomix_s(dimsh,nitmax)
    histomix_s(:,:) = 0
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol-1) then
      iti = itype(i)
       xi = r(1,i)
       yi = r(2,i)
       Do j=i+1,Nmol
          xd = r(1,j)-xi
          yd = r(2,j)-yi
          xd = xd -sidel(1)*nint(xd/sidel(1))
          yd = yd -sidel(2)*nint(yd/sidel(2))
          rr2= xd*xd+yd*yd
          If (rr2.Lt.side2) Then
             itj = itype(j)
             rr = __fsqrt_rn(rr2)
             ind = Nint(rr/deltar)
             ij = fij(iti,itj,nsp)
             !
             ! Use atomics over shared memory to minimize collisions
             !
             ia = atomicadd(histomix_s(ind,ij),1)
          !   if (iti /= itj) ia = atomicadd(histomix_s(ind,ij),1)
          endif
       Enddo
       call syncthreads()
       if (threadidx%x == 1) then
         istart = (blockidx%x-1)*lsmax
         !
         ! Store each block's shared memory histogram on different positions
         ! in global memory
         !
         histomix(istart+1:istart+lsmax,:) = histomix_s(1:lsmax,:)+histomix(istart+1:istart+lsmax,:)
       endif 
    end if
  end subroutine rdf2_sh



  attributes (global) subroutine dprof(r,Nmol,dim,densprof,width,nsp,idir,itype&
       &,pwall,deltar)
    integer, value, intent(IN) :: Nmol, dim, nsp, width, idir
    integer, intent(IN) :: itype(Nmol)
    integer i, ind, ia, ity
    real, dimension(3) :: rv
    real, value, intent(IN) :: deltar, pwall
    real, dimension(dim,Nmol), intent(IN) :: r
    integer, intent(INOUT) :: densprof(width,nsp)
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol) then
       !
       ! Warning: density profile must be computed along a confined direction or using wrapped coordinates
       !
       ind = int((r(idir,i)-pwall)/(2*deltar))
       if (ind < 0 .or. ind > width) then
          print *, ' Error: density profile out of boundaries ', ind, width
          stop
       endif
       ity = itype(i)
       ia = atomicadd(densprof(ind,ity),1)
    endif
  end subroutine dprof

  
  
  attributes (global) subroutine sqfact(r_d,Nmol,dim,itype,ntype,ipos,fk,sqf,sqfp,nsp,kmax,kmin)
    integer, value, intent(IN) :: Nmol, dim, nsp,  kmax, kmin
    integer, intent(IN) :: itype(Nmol), ntype(nsp),ipos(nsp,Nmol)
    real, dimension(3,Nmol), intent(IN) :: r_d
    real (kind=4), intent(INOUT) :: sqf(kmax),sqfp(kmax,nsp)
    real (kind=4) :: rkx, rky, rkz, fk1, fk2, fk3, kf1, kf2, kf3, suma
    real (kind=8) :: sumcx, sumsx, sumcy, sumsy, sumcz, sumsz,  &
         & sumcpx, sumspx, sumcpy, sumspy, sumcpz, sumspz
    real (kind=4), intent(IN) :: fk(dim)
    integer :: i,j,k, ind, it, ia
    k = (blockidx%x-1) * blockdim%x + threadidx%x
    if (k .le. kmax .and. k .ge.kmin) then
       fk1 = fk(1)
       fk2 = fk(2)
       fk3 = fk(3)
       sumcx = 0
       sumsx = 0
       sumcy = 0
       sumsy = 0
       sumcz = 0
       sumsz = 0
       kf1 = k*fk1
       kf2 = k*fk2
       kf3 = k*fk3
       i = 0
       do j=1, nsp
          sumcpx = 0
          sumspx = 0
          sumcpy = 0
          sumspy = 0
          sumcpz = 0
          sumspz = 0
          do it=1, ntype(j)
             i = i + 1
             rkx = r_d(1,i)*kf1
             rky = r_d(2,i)*kf2
             rkz = r_d(3,i)*kf3
             sumcpx = sumcpx + __cosf(rkx)
             sumspx = sumspx + __sinf(rkx)
             sumcpy = sumcpy + __cosf(rky)
             sumspy = sumspy + __sinf(rky)
             sumcpz = sumcpz + __cosf(rkz)
             sumspz = sumspz + __sinf(rkz)
          end do
          sumcx = sumcx+sumcpx
          sumsx = sumsx+sumspx
          sumcy = sumcy+sumcpy
          sumsy = sumsy+sumspy
          sumcz = sumcz+sumcpz
          sumsz = sumsz+sumspz
          suma = (sumcpx**2+sumspx**2+sumcpy**2+sumspy**2+sumcpz**2+sumspz**2)
          ia = atomicadd(sqfp(k,j),suma)
         end do
       sqf(k) = sqf(k)+(sumcx**2+sumsx**2+sumcy**2+sumsy**2+sumcz**2+sumsz**2)
    endif
  end subroutine sqfact

  attributes (global) subroutine sqfact2(r_d,Nmol,dim,itype,ntype,ipos,fk,sqf,sqfp,nsp,kmax,kmin)
    integer, value, intent(IN) :: Nmol, dim, nsp,  kmax, kmin
    integer, intent(IN) :: itype(Nmol), ntype(nsp),ipos(nsp,Nmol)
    real, dimension(2,Nmol), intent(IN) :: r_d
    real (kind=4), intent(INOUT) :: sqf(kmax),sqfp(kmax,nsp)
    real (kind=4) :: rkx, rky, rkz, fk1, fk2, fk3, kf1, kf2, kf3, suma
    real (kind=8) :: sumcx, sumsx, sumcy, sumsy, sumcpx, sumspx, sumcpy, sumspy
    real (kind=4), intent(IN) :: fk(dim)
    integer :: i,j,k, ind, it, ia
    k = (blockidx%x-1) * blockdim%x + threadidx%x
    if (k .le. kmax .and. k .ge.kmin) then
       fk1 = fk(1)
       fk2 = fk(2)
       sumcx = 0
       sumsx = 0
       sumcy = 0
       sumsy = 0
       kf1 = k*fk1
       kf2 = k*fk2
       i = 0
       do j=1, nsp
          sumcpx = 0
          sumspx = 0
          sumcpy = 0
          sumspy = 0
          do it=1, ntype(j)
             i = i + 1
             rkx = r_d(1,i)*kf1
             rky = r_d(2,i)*kf2
             sumcpx = sumcpx + cos(rkx)
             sumspx = sumspx + sin(rkx)
             sumcpy = sumcpy + cos(rky)
             sumspy = sumspy + sin(rky)
          end do
          sumcx = sumcx+sumcpx
          sumsx = sumsx+sumspx
          sumcy = sumcy+sumcpy
          sumsy = sumsy+sumspy
          suma = (sumcpx**2+sumspx**2+sumcpy**2+sumspy**2)
          ia = atomicadd(sqfp(k,j),suma)
       end do
       sqf(k) = sqf(k)+(sumcx**2+sumsx**2+sumcy**2+sumsy**2)
    endif
  end subroutine sqfact2

  attributes (global) subroutine sqf3D(r_d,Nmol,dim,itype,ntype,ipos,fk,sqf,sqfp,nsp,kmax,q2max,dq)
    integer, intent(IN) :: itype(Nmol),  ntype(nsp), ipos(nsp,Nmol)
    integer, value, intent(IN) :: Nmol, dim, nsp, kmax
    real, dimension(3,Nmol), intent(IN) :: r_d
    real (kind=4), intent(INOUT) :: sqf(kmax),sqfp(kmax,nsp)
    real (kind=4) :: rkx, rky, rkz, fk1, fk2, fk3, kf1, kf2, kf3, tsum, q, suma
    real (kind=8) :: sumcx, sumsx, sumcpx, sumspx, q2
    real (kind=4), intent(IN) :: fk(dim)
    real (kind=4), value, intent(IN) ::  dq, q2max
    
    integer :: i ,j, it, kx, ky, kz, ind, ia
    kx = (blockidx%x-1)*blockdim%x + threadidx%x-1
    ky = (blockidx%y-1)*blockdim%y + threadidx%y-1
    kz = (blockidx%z-1)*blockdim%z + threadidx%z-1
    fk1 = fk(1)
    fk2 = fk(2)
    fk3 = fk(3)
    kf1 = kx*fk1
    kf2 = ky*fk2
    kf3 = kz*fk3
    q2 = kf1*kf1+kf2*kf2+kf3*kf3
    if (q2>0.0 .and. q2 <= q2max) then
       sumcx = 0
       sumsx = 0
       q = sqrt(q2)
       ind = nint(q/dq)
       i = 0
       do j=1, nsp
          sumcpx = 0
          sumspx = 0
          do it=1, ntype(j)
             i = i + 1
             rkx = r_d(1,i)*kf1
             rky = r_d(2,i)*kf2
             rkz = r_d(3,i)*kf3
             sumcpx = sumcpx + __cosf(rkx+rky+rkz)
             sumspx = sumspx + __sinf(rkx+rky+rkz)
          end do
          sumcx = sumcx+sumcpx
          sumsx = sumsx+sumspx
          suma = (sumcpx**2+sumspx**2)
          ia = atomicadd(sqfp(ind,j),suma)
        end do
       tsum = (sumcx**2+sumsx**2)
       ia = atomicadd(sqf(ind),tsum)
    endif
  end subroutine sqf3D

    attributes (global) subroutine sqf2D(r_d,Nmol,dim,itype,ntype,ipos,fk,sqf,sqfp,nsp,kmax,q2max,dq)
    integer, intent(IN) :: itype(Nmol),  ntype(nsp), ipos(nsp,Nmol)
    integer, value, intent(IN) :: Nmol, dim, nsp, kmax
    real, dimension(2,Nmol), intent(IN) :: r_d
    real (kind=4), intent(INOUT) :: sqf(kmax),sqfp(kmax,nsp)
    real (kind=4) :: rkx, rky, rkz, fk1, fk2, fk3, kf1, kf2, kf3, tsum, q, suma
    real (kind=8) :: sumcx, sumsx, sumcpx, sumspx, q2
    real (kind=4), intent(IN) :: fk(dim)
    real (kind=4), value, intent(IN) ::  dq, q2max
    
    integer :: i ,j, it, kx, ky, kz, ind, ia
    kx = (blockidx%x-1)*blockdim%x + threadidx%x-1
    ky = (blockidx%y-1)*blockdim%y + threadidx%y-1
    fk1 = fk(1)
    fk2 = fk(2)
    kf1 = kx*fk1
    kf2 = ky*fk2
    q2 = kf1*kf1+kf2*kf2
    if (q2>0.0 .and. q2 <= q2max) then
       sumcx = 0
       sumsx = 0
       q = sqrt(q2)
       ind = nint(q/dq)
       i = 0 
       do j=1, nsp
          sumcpx = 0
          sumspx = 0
          do it=1, ntype(j)
             i = i +1
             rkx = r_d(1,i)*kf1
             rky = r_d(2,i)*kf2
             sumcpx = sumcpx + sin(rkx+rky+1.5707963267948966)
             sumspx = sumspx + sin(rkx+rky)
          end do
          sumcx = sumcx+sumcpx
          sumsx = sumsx+sumspx
          suma = (sumcpx**2+sumspx**2)
          ia = atomicadd(sqfp(ind,j),suma)
        end do
       tsum = (sumcx**2+sumsx**2)
       ia = atomicadd(sqf(ind),tsum)
    endif
  end subroutine sqf2D

    attributes (global) subroutine sqf3Dcl(r_d,Nmol,fk,sqf,kmax,q2max,dq)
    integer, value, intent(IN) :: Nmol, kmax
    real, dimension(3,Nmol), intent(IN) :: r_d
    real (kind=4), intent(INOUT) :: sqf(kmax)
    real (kind=4) :: rkx, rky, rkz, fk1, fk2, fk3, kf1, kf2, kf3, tsum, q, suma
    real (kind=8) :: sumcx, sumsx, sumcpx, sumspx, q2
    real (kind=4), intent(IN) :: fk(3)
    real (kind=4), value, intent(IN) ::  dq, q2max
    
    integer :: i ,j, it, kx, ky, kz, ind, ia
    kx = (blockidx%x-1)*blockdim%x + threadidx%x-1
    ky = (blockidx%y-1)*blockdim%y + threadidx%y-1
    kz = (blockidx%z-1)*blockdim%z + threadidx%z-1
    fk1 = fk(1)
    fk2 = fk(2)
    fk3 = fk(3)
    kf1 = kx*fk1
    kf2 = ky*fk2
    kf3 = kz*fk3
    q2 = kf1*kf1+kf2*kf2+kf3*kf3
    if (q2>0.0 .and. q2 <= q2max) then
       sumcx = 0
       sumsx = 0
       q = sqrt(q2)
       ind = nint(q/dq)
       sumcpx = 0
       sumspx = 0
       do i=1, Nmol
          rkx = r_d(1,i)*kf1
          rky = r_d(2,i)*kf2
          rkz = r_d(3,i)*kf3
          sumcpx = sumcpx + cos(rkx+rky+rkz)
          sumspx = sumspx + sin(rkx+rky+rkz)
       end do
       suma = (sumcpx**2+sumspx**2)/Nmol
       ia = atomicadd(sqf(ind),suma)
    endif
  end subroutine sqf3Dcl

  attributes (global) subroutine sqf2Dcl(r_d,Nmol,fk,sqf,kmax,q2max,dq)
    integer, value, intent(IN) :: Nmol, kmax
    real, dimension(2,Nmol), intent(IN) :: r_d
    real (kind=4), intent(INOUT) :: sqf(kmax)
    real (kind=4) :: rkx, rky, rkz, fk1, fk2, fk3, kf1, kf2, kf3, tsum, q, suma
    real (kind=8) :: sumcx, sumsx, sumcpx, sumspx, q2
    real (kind=4), intent(IN) :: fk(3)
    real (kind=4), value, intent(IN) ::  dq, q2max
    
    integer :: i ,j, it, kx, ky, kz, ind, ia
    kx = (blockidx%x-1)*blockdim%x + threadidx%x-1
    ky = (blockidx%y-1)*blockdim%y + threadidx%y-1
    fk1 = fk(1)
    fk2 = fk(2)
    kf1 = kx*fk1
    kf2 = ky*fk2
    q2 = kf1*kf1+kf2*kf2
    if (q2>0.0 .and. q2 <= q2max) then
       sumcx = 0
       sumsx = 0
       q = sqrt(q2)
       ind = nint(q/dq)
       sumcpx = 0
       sumspx = 0
       do i=1, Nmol
          rkx = r_d(1,i)*kf1
          rky = r_d(2,i)*kf2
          sumcpx = sumcpx + cos(rkx+rky)
          sumspx = sumspx + sin(rkx+rky)
       end do
       suma = (sumcpx**2+sumspx**2)/Nmol
       ia = atomicadd(sqf(ind),suma)
    endif
  end subroutine sqf2Dcl

  function gptime(stopEvent,startEvent)
    type(cudaEvent), intent(IN) :: startEvent, stopEvent
    integer :: istat
    real :: time, gptime
    istat = cudaEventRecord(stopEvent,0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time, startEvent, stopEvent)
    gptime = time/1000.0
  end function gptime
!   attributes (global) subroutine rdf(r,Nmol,dim,histomix,nsp,lsmax,itype&
!        &,side2,sidel,deltar)
!     integer, value, intent(IN) :: Nmol, dim, nsp, lsmax
!     integer, intent(IN) :: itype(Nmol)
!     integer i, j, ind, ia, fact, iti, itj
   !  real ::  rr2, rr, xi, yi, zi, xd, yd, zd
   !  real, dimension(3) :: rv
   !  real, value, intent(IN) :: side2, deltar
   !  real, intent(IN) :: sidel(3)
   !  real, dimension(dim,Nmol), intent(IN) :: r
   !  integer, intent(INOUT) :: histomix(lsmax,nsp,nsp)
   !  i = (blockidx%x-1) * blockdim%x + threadidx%x
!     if (i<=Nmol-1) then
!       iti = itype(i)
!        xi = r(1,i)
!        yi = r(2,i)
!        zi = r(3,i)
!        Do j=i+1,Nmol
!           xd = r(1,j)-xi
!           yd = r(2,j)-yi
!           zd = r(3,j)-zi
!           xd = xd -sidel(1)*nint(xd/sidel(1))
!           yd = yd -sidel(2)*nint(yd/sidel(2))
!           zd = zd -sidel(3)*nint(zd/sidel(3))
!           rr2= xd*xd+yd*yd+zd*zd
!           If (rr2.Lt.side2) Then
!              itj = itype(j)
!              rr = __fsqrt_rn(rr2)
!              ind = Nint(rr/deltar)
!              ia = atomicadd(histomix(ind,itj,iti),1)
!              if (iti /= itj) ia = atomicadd(histomix(ind,iti,itj),1)
!           endif
!        Enddo
!     end if
!   end subroutine rdf

!   attributes (global) subroutine rdf2(r,Nmol,dim,histomix,nsp,lsmax,itype&
!        &,side2,sidel,deltar)
!     integer, value, intent(IN) :: Nmol, dim, nsp, lsmax
!     integer, intent(IN) :: itype(Nmol)
!     integer i, j, ind, ia, fact, iti, itj
!     real ::  rr2, rr, xi, yi,  xd, yd
!     real, dimension(3) :: rv
!     real, value, intent(IN) :: side2, deltar
!     real, intent(IN) :: sidel(3)
!     real, dimension(dim,Nmol), intent(IN) :: r
!     integer, intent(INOUT) :: histomix(lsmax,nsp,nsp)
!     i = (blockidx%x-1) * blockdim%x + threadidx%x
!     if (i<=Nmol-1) then
!     iti = itype(i)
!        xi = r(1,i)
!        yi = r(2,i)
!        Do j=i+1,Nmol
!           xd = r(1,j)-xi
!           yd = r(2,j)-yi
!           xd = xd -sidel(1)*nint(xd/sidel(1))
!           yd = yd -sidel(2)*nint(yd/sidel(2))
!           rr2= xd*xd+yd*yd
!           If (rr2.Lt.side2) Then
!              itj = itype(j)
!              rr = __fsqrt_rn(rr2)
!              ind = Nint(rr/deltar)
!              ia = atomicadd(histomix(ind,itj,iti),1)
!              if (iti /= itj) ia = atomicadd(histomix(ind,iti,itj),1)
!           endif
!        Enddo
!     end if
!   end subroutine rdf2
end module gpcodes
