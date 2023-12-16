module dbscan
  implicit none
  integer, allocatable :: offset(:), neighbors(:), adjacency(:)
  logical (kind=1), allocatable :: visited(:), core(:), border(:)
  
contains
  subroutine init_dbscan(Nmol)
    integer, intent(IN) :: Nmol
    allocate(offset(Nmol), neighbors(Nmol), visited(Nmol),core(Nmol), border(Nmol), adjacency(20*Nmol))
    visited(:) = .false.
    core(:) = .false.
    border(:) = .false.
    neighbors(:) = 0
  end subroutine init_dbscan
  
  subroutine build_graph(r,Nmol,dim, rcl2,sidel,jmin)
    integer, intent(IN) :: Nmol, dim, jmin
    real(kind=4), intent(IN) :: r(dim,Nmol), rcl2, sidel(dim)
    real (kind=4) :: rv(dim), rr2
    integer :: i, j, counter
    offset(1) = 1
    do i=1, Nmol
       counter = 0
       do j=1, Nmol
          if (i.ne.j) then
             rv(1:dim) = r(1:dim,i)-r(1:dim,j)
             ! PBC (unwrap if necessary)
             rv = rv-sidel*nint(rv/sidel)
             rr2= Dot_product(rv,rv)
             if (rr2 < rcl2) then
                adjacency(offset(i)+counter) = j
                counter = counter+1
             endif
          end if
       enddo
       neighbors(i) = counter
       if (counter >= jmin-1) core(i)=.true.
       if (i<=Nmol-1) offset(i+1) = counter+offset(i)
    enddo
  end subroutine build_graph

  subroutine bfs(Nmol, done)
    use clusters, only: cluster
    integer, intent(IN) :: Nmol
    integer, intent(INOUT) :: done(Nmol)
    integer :: i, j, k, l, im, nb, counter, nbo, seen, cluster_id, ic, next(nmol), nexto(100)
    logical :: first
    cluster_id =0
    seen = 1
    do i=1, nmol
       next(:) = 0
       nbo = 1
       next(1) = i
       if (.not.visited(i) .and. core(i)) then
          cluster_id = cluster_id+1
          im = 0
          done(i) = cluster_id 
          l=1
          visited(i) = .true.
          done(seen) = i
          nexto(:) = 0
          do while (next(l)>0)
             ic = next(l)
             counter = 0
             do j = 0, neighbors(ic)-1
                k = adjacency(offset(ic)+j)
                if (.not.visited(k)) then
                   seen = seen+1
                   im = im+1
                   done(k) = cluster_id
                   visited(k) = .true.
                   counter=counter+1
                   nexto(counter) = k
                endif
             enddo
             nb = counter
             next(nbo+1:nbo+1+nb) = nexto(1:nb)
             nbo = nbo+nb
             l=l+1
          end do
!          cluster(cluster_id)%clsize = im
       endif
    end do
   end subroutine bfs
end module dbscan

