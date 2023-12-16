module cpucodes
  use comun, only : r, natms, pi
  implicit none

  Interface
     subroutine cluster_cpu(list,Nmol,dim, rcl2,side2,sidel)
       integer, intent(IN) :: Nmol, dim
       integer, intent(INOUT) :: list(Nmol)
       real(kind=4), intent (IN) :: rcl2, sidel(dim), side2
     end subroutine cluster_cpu
     subroutine rdfcl(dim,members,nm,sidel,side2,g,lsmax,deltar,vol,jmin)
       integer, intent(IN) :: nm, lsmax,jmin, dim
       integer, intent(IN) :: members(nm)
       real, intent(IN) :: sidel(3), side2, deltar, vol
       real, intent(INOUT) :: g(lsmax)
     end subroutine rdfcl
     subroutine rdfclcl(dim,centers,nm,sidel,side2,g,lsmax,deltar,vol)
       integer, intent(IN) :: nm, lsmax, dim
       real, intent(IN) :: centers(3,nm),sidel(3), side2, deltar, vol
       real, intent(INOUT) :: g(lsmax)
     end subroutine rdfclcl
     
     function cldens(dim,members,nm,sidel,side2,center,rclus)
       integer, intent(IN) :: nm, dim
       integer, intent(IN) :: members(nm)
       real, intent(IN) :: sidel(3), side2, center(3), rclus
       real :: avcldens
     end function cldens

     function minj(dim,members,nm,center,sidel,side2)
       integer, intent(IN) :: nm, dim
       integer, intent(IN) :: members(nm)
       real, intent(IN) :: center(3),sidel(3), side2
     end function minj

     subroutine denspcl(dim,members,nm,sidel,side2,rho,lsmax,deltar,center)
       integer, intent(IN) :: nm, lsmax, dim
       integer, intent(IN) :: members(nm)
       real, intent(IN) :: sidel(3), center(3), side2, deltar
       real, intent(INOUT) :: rho(lsmax)
     end subroutine denspcl
  end Interface
  
contains
  subroutine denspcl(dim,members,nm,sidel,side2,rho,lsmax,deltar,center)
    integer, intent(IN) :: nm, lsmax, dim
    integer, intent(IN) :: members(nm)
    real, intent(IN) :: sidel(3), center(3), side2, deltar
    real, intent(INOUT) :: rho(0:lsmax)
    real :: rv(3), rr2, rr, fact
    integer :: k, l, i, j, ind
    ! g(r) is normalized for each configuration, since particle number may change
    do k = 1, nm
       i = members(k)
       rv(1:dim) = r(1:dim,i)-center(1:dim)
       ! PBC (unwrap if necessary)
       rv(1:dim) = rv(1:dim)-sidel(1:dim)*nint(rv(1:dim)/sidel(1:dim))
       rr2= Dot_product(rv(1:dim),rv(1:dim))
       If (rr2.Lt.side2) Then
          rr = sqrt(rr2)
          ind = Nint(rr/deltar)
          rho(ind) = rho(ind)+1.0
       end If
    end do
  end subroutine denspcl


subroutine rdfcl(dim,members,nm,sidel,side2,g,lsmax,deltar,dens,jmin)
    integer, intent(IN) :: nm, lsmax, jmin, dim
    integer, intent(IN) :: members(nm)
    real, intent(IN) :: sidel(3), side2, deltar, dens
    real, intent(INOUT) :: g(lsmax)
    real :: rv(3), rr2, rr, fact
    integer :: k, l, i, j, ind
    ! g(r) is normalized for each configuration, since particle number may change
    fact = 1.0/dens
    do k = 1, nm
       i = members(k)
       if (k .ne. jmin) then 
          j = members(jmin)
          rv(1:dim) = r(1:dim,i)-r(1:dim,j)
          ! PBC (unwrap if necessary)
          rv(1:dim) = rv(1:dim)-sidel(1:dim)*nint(rv(1:dim)/sidel(1:dim))
          rr2= Dot_product(rv(1:dim),rv(1:dim))
          If (rr2.Lt.side2) Then
             rr = sqrt(rr2)
             ind = Nint(rr/deltar)
             g(ind) = g(ind)+fact
          end If
       endif
    end do
  end subroutine rdfcl

  function minj(dim,members,nm,center,sidel,side2)
    integer, intent(IN) :: nm, dim
    integer, intent(IN) :: members(nm)
    real, intent(IN) :: center(3),sidel(3), side2
    real :: rv(3), dist, r2
    integer :: j, i, jmin, minj
    dist = 100000000.0
    do j=1, nm
       i = members(j)
       rv(1:dim) = r(1:dim,i)-center(1:dim)
       rv(1:dim) = rv(1:dim) - sidel(1:dim)*nint(rv(1:dim)/sidel(1:dim))
       r2 = dot_product(rv(1:dim),rv(1:dim))
       if (r2<= dist) then
          dist = r2
          jmin = j
       end if
    end do
    minj = jmin
  end function minj
  
  function cldens(dim,members,nm,sidel,side2,center,rclus)
    integer, intent(IN) :: nm, dim
    integer, intent(IN) :: members(nm)
    real, intent(IN) :: sidel(3), side2, center(3), rclus
    real :: avcldens
    real :: rv(3), vol, rr2, rclus2
    integer :: i, j, count
    if (dim == 3) then
       vol = 4*pi*rclus**3/3.0
    else
       vol = 4*pi*rclus**2
    end if
    rclus2 = rclus*rclus
    count = 0
    do i = 1, nm
       j = members(i)
       rv(1:dim) = r(1:dim,j)-center(1:dim)
       ! PBC (unwrap if necessary)
       rv(1:dim) = rv(1:dim)-sidel(1:dim)*nint(rv(1:dim)/sidel(1:dim))
       rr2= Dot_product(rv(1:dim),rv(1:dim))
       if (rr2 < rclus2) then
          count = count +1
       endif
    end do
    cldens = real(count)/vol
  end function cldens
  
  subroutine rdfclcl(dim,centers,nm,sidel,side2,g,lsmax,deltar,vol)
    integer, intent(IN) :: nm, lsmax, dim
    real, intent(IN) :: centers(3,nm),sidel(3), side2, deltar, vol
    real, intent(INOUT) :: g(lsmax)
    real :: rv(3), rr2, rr, fact
    integer :: k, l, i, j, ind
    ! g(r) is normalized for each configuration, since particle number may change 
    fact = vol/real(nm*nm)
    do i = 1, nm
       do j = i+1, nm
          rv(1:dim) = centers(1:dim,i)-centers(1:dim,j)
          ! PBC (unwrap if necessary)
          rv(1:dim) = rv(1:dim)-sidel*nint(rv(1:dim)/sidel(1:dim))
          rr2= Dot_product(rv(1:dim),rv(1:dim))
          If (rr2.Lt.side2) Then
             rr = sqrt(rr2)
             ind = Nint(rr/deltar)
             g(ind) = g(ind)+fact
          end If
       end do
    end do
  end subroutine rdfclcl

  subroutine cluster_cpu(list,Nmol,dim, rcl2,side2,sidel)
    integer, intent(IN) :: Nmol, dim
    integer, intent(INOUT) :: list(Nmol)
    real(kind=4), intent (IN) :: rcl2, sidel(dim), side2
    real (kind=4) :: rv(dim), rr2
    integer :: i, j, k

    list(:) = [ (i,i=1,Nmol) ] ! Set up the list

    DO i = 1, Nmol - 1 ! Begin outer loop

       IF ( i == list(i) ) THEN

        j = i

        DO ! Begin inner loop

           DO k = i + 1, Nmol ! Begin innermost loop

              IF ( list(k) == k ) THEN
                  rv(1:dim) = r(1:dim,j)-r(1:dim,k)
                  ! PBC (unwrap if necessary)
                  rv = rv-sidel*nint(rv/sidel)
                  rr2= Dot_product(rv,rv)
                  IF ( rr2 <= rcl2 ) list([k,j]) = list([j,k]) ! swap elements
               END IF

            END DO ! End innermost loop

            j = list(j)
            IF ( j == i ) EXIT

         END DO ! End inner loop

      END IF

   END DO ! End outer loop


 end subroutine cluster_cpu
 
end module cpucodes
