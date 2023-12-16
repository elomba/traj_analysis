Module linkcell
  use g_types
  !
  ! Shared components in the link cell method
  !
  Real (myprec) :: cellx, celly, cellz
  Integer, Dimension(:,:), Allocatable:: neigh
  Integer, Dimension(:), Allocatable :: head, lista
  Integer :: ncell, nn, maxi,maxj,maxk
  logical :: use_cell=.true.
End Module linkcell

Module linkcell_d
  use g_types
  !
  ! Shared components in the link cell method
  !
  Real(myprec), device :: cellxd, cellyd, cellzd
  Integer, device, Dimension(:,:), Allocatable:: neighd
  Integer, device,  Dimension(:), Allocatable :: headd, listad
  Integer, device :: ncelld, nnd, maxid,maxjd,maxkd
End Module linkcell_d


module cells
  use comun, only : r, sidel, natms, ndim
  implicit none
contains
  Subroutine Init_cell(rcl)
    use g_types
    Use linkcell
    Use linkcell_d
    Implicit None
    !
    !   Only for orthogonal cells
    !
    Integer :: i,j,k,ii,jj,kk,l,in
    real(myprec), intent(IN) :: rcl
    !
    ! Cell size extends up to rcut
    !
    maxi= Int(sidel(1)/rcl)
    maxj= int(sidel(2)/rcl)
    maxk = 1
    if (ndim==3) maxk= Int(sidel(3)/rcl)
    ncell = maxi*maxj*maxk
    nn = 3**ndim
    maxid = maxi
    maxjd = maxj
    maxkd = maxk
    ncelld = ncell
    nnd = nn
    Allocate(neigh(0:ncell-1,nn),head(0:ncell-1),lista(natms))
    Allocate(neighd(0:ncell-1,nn),headd(0:ncell-1),listad(natms))
    cellx = sidel(1)/maxi
    celly = sidel(2)/maxj
    cellz = sidel(3)/maxk
    cellxd = cellx
    cellyd = celly
    cellzd = cellz
    if (min(maxi,maxj,maxk) < 3) then
       print *, " *** Error: box too small for link cell method "
       use_cell = .false.
       return
    end if
    l=0
    if (ndim ==3) then
       Do i=0,maxi-1
          Do j=0,maxj-1
             Do k=0,maxk-1
                in = 1
                Do ii=-1,1
                   Do jj=-1,1
                      Do kk=-1,1
                         neigh(l,in) = fijk(i+ii,j+jj,k+kk)
                         in = in+1
                      End Do
                   End Do
                End Do
                l=l+1
             End Do
          End Do
       End Do
    else
       Do i=0,maxi-1
          Do j=0,maxj-1
             in = 1
             Do ii=-1,1
                Do jj=-1,1
                   neigh(l,in) = fij(i+ii,j+jj)
                   in = in+1
                End Do
             End Do
             l=l+1
          End Do
       End Do
    end if
    neighd(:,:) = neigh(:,:)
  end Subroutine Init_cell

  Integer Function fijk(ix,jx,kx)
    use linkcell, only : maxi, maxj, maxk
    Implicit None
    Integer :: ix, jx, kx, i, j, k
    i = ix
    j = jx
    k = kx
    !
    !  Use periodic boundary conditions
    !
    if (i < 0) i=i+maxi
    if (j < 0) j=j+maxj
    if (k < 0) k=k+maxk
    if (i >= maxi) i=i-maxi
    if (j >= maxj) j=j-maxj
    if (k >= maxk) k=k-maxk
    fijk=(i*maxj+j)*maxk+k
  end function fijk

  Integer Function fij(ix,jx)
    Use linkcell, Only : maxi, maxj
    Implicit None
    Integer :: ix, jx, i, j
    i = ix
    j = jx
    !
    !  Use periodic boundary conditions
    !
    if (i < 0) i=i+maxi
    if (j < 0) j=j+maxj
    if (i >= maxi) i=i-maxi
    if (j >= maxj) j=j-maxj
    fij=(i*maxj+j)
  end function fij


  attributes(device) Integer Function fijkd(ix,jx,kx)
    use linkcell_d, only : maxid, maxjd, maxkd
    Implicit None
    Integer :: ix, jx, kx, i, j, k
    i = ix
    j = jx
    k = kx
    !
    !  Use periodic boundary conditions
    !
    if (i < 0) i=i+maxid
    if (j < 0) j=j+maxjd
    if (k < 0) k=k+maxkd
    if (i >= maxid) i=i-maxid
    if (j >= maxjd) j=j-maxjd
    if (k >= maxkd) k=k-maxkd
    fijkd=(i*maxjd+j)*maxkd+k
  end function fijkd

  attributes(device) Integer Function fijd(ix,jx)
    use linkcell_d, only : maxid, maxjd
    Integer :: ix, jx, i, j
    i = ix
    j = jx
    !
    !  Use periodic boundary conditions
    !
    if (i < 0) i=i+maxid
    if (j < 0) j=j+maxjd
    if (i >= maxid) i=i-maxid
    if (j >= maxjd) j=j-maxjd
    fijd=(i*maxjd+j)
  end function fijd

  subroutine build_cells
    use linkcell
    use linkcell_d
    implicit none
    integer :: n, i, j, k, icell
    head(:) = 0
    lista(:) = 0

    do n = 1, natms
       i=int((R(1,n))/cellx)
       if (i<0 .or. i>= maxi) then
          print *, ' error in cells in x direction for ', n, R(1,n),cellx
          stop
       endif
       j=int((R(2,n))/celly)
       if (j<0 .or. j>= maxj) then
          print *, ' error in cells in y direction for ', n, R(2,n),celly
          stop
       endif
       if (ndim == 3) then
          k=int((R(3,n))/cellz)
          if (k<0 .or. k>= maxk) then
             print *, ' error in cells in z direction for ', n, R(3,n),cellz
             stop
          endif
          icell=(i*maxj+j)*maxk+k
       else
          icell=i*maxj+j
       endif
       lista(n) = head(icell)
       head(icell) = n
    end do
    listad(1:natms) = lista(1:natms)
    headd(0:ncell-1) =  head(0:ncell-1)
  end subroutine build_cells


end module cells


module gpucells

  use cudafor
  use linkcell_d
  use gpcodes
  implicit none
contains

  attributes (global) subroutine gpu_graph_cell(r,Nmol,dim, rcl2,sidel,jmin)
    use g_types
    integer, value, intent(IN) :: Nmol, dim, jmin
    real(myprec), value :: rcl2
    real(myprec), intent(IN) :: r(dim,Nmol), sidel(dim)
    real(myprec)  :: rv(3), rr2
    integer :: i, j, counter, ia, ix, jy, kz, cell, icell
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol) then
       gneighbors(i) = 0
       counter = 0
       !
       ! Locate cell
       !
       ix=int((r(1,i))/cellxd)
       jy=int((r(2,i))/cellyd)
       if (dim==3) then
          kz=int((r(3,i))/cellzd)
          cell=(ix*maxjd+jy)*maxkd+kz
       else
          cell = ix*maxjd+jy
       endif
       !
       ! Loop over neighbouring cells
       !
       do icell = 1, nnd
          j = headd(neighd(cell,icell))
          do while (j .ne. 0)
             if (j.ne.i) then
                rv(1:dim) = r(1:dim,i)-r(1:dim,j)
                rv = rv-sidel*nint(rv/sidel)
                rr2= Dot_product(rv,rv)

                if (rr2 < rcl2) then
                   counter = counter+1
                endif
             end if
             j = listad(j)
          end do
       end do
       gneighbors(i) = counter
       if (counter >= jmin-1)  gcore(i)= .true.
    endif
  end subroutine gpu_graph_cell

  attributes (global) subroutine gpu_adj_cell(r,Nmol,dim, rcl2,sidel)
    use g_types
    integer, value, intent(IN) :: Nmol, dim
    real(myprec), value :: rcl2
    real(myprec), intent(IN) :: r(dim,Nmol), sidel(dim)
    real(myprec)  :: rv(3), rr2
    integer :: i, j, counter, ia, ix, jy, kz, cell, icell
    i = (blockidx%x-1) * blockdim%x + threadidx%x
    if (i<=Nmol) then
       counter = 0
       !
       ! Locate cell
       !
       ix=int((r(1,i))/cellxd)
       jy=int((r(2,i))/cellyd)
       if (dim==3) then
          kz=int((r(3,i))/cellzd)
          cell=(ix*maxjd+jy)*maxkd+kz
       else
          cell=(ix*maxjd+jy)
       end if
       !
       ! Loop over neighbouring cells
       !
       do icell = 1, nnd
          j = headd(neighd(cell,icell))
          do while (j .ne. 0)
             if (j.ne.i) then
                rv(1:dim) = r(1:dim,i)-r(1:dim,j)
                rv = rv-sidel*nint(rv/sidel)
                rr2= Dot_product(rv,rv)
                if (rr2 < rcl2) then
                   gadjacency(goffset(i)+counter) = j
                   counter = counter+1
                endif
             end if
             j = listad(j)
          end do
       end do
    endif
  end subroutine gpu_adj_cell

end module gpucells
