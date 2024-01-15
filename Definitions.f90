Module Comun
  use g_types
  Interface
     Function vector_product(a,b)
       use g_types
       Real (myprec), Dimension(:), Intent(IN) :: a, b
       Real (myprec), Dimension(Size(a)) :: vector_product
     End Function vector_product
     Function surface(a,b)
       use g_types
       Real (myprec), Dimension(:), Intent(IN) :: a, b
       Real (myprec) :: surface
     End Function surface
     Subroutine ReadCfg(iunit,io,iconf,istart)
       Integer, Intent(IN) :: iunit,iconf,istart
       Integer, Intent(OUT) :: io
     end Subroutine ReadCfg
     function radius(members, nm, center,sidel)
       use g_types
       integer, intent(IN) :: nm
       integer, intent(IN) :: members(nm)
       real(myprec), intent(IN) :: center(3), sidel(3)
     end function radius
  End Interface
  !
  !  Module of common parameters. Change here the maximum no. of points
  !
  logical UnitR
  Integer, Parameter :: npmx=8192, nmblock=500,nmgrid=50000,&
       & nspmax=3, nmolmax=220000,nmesmax=10000,nkmax=200,ksmax=6, dimsh=1000
  Integer, Dimension(:), Allocatable :: ntype, itype,  imolt, nctr, ntml, imf, iml, nml, tmol,  list, done
  Integer :: Nconf, natms, imcon, keytrj, nsp, itraj, nmolty, nmzero, ndim
  Integer ::  nqmax, lmaxx, lmaxy, lmaxz, nqmin, lsmax, ndr, width, idir
  Integer, save :: nstep=0
  Real (myprec), Dimension(:,:), Allocatable :: vel, force, r, msq, sft
  Real (myprec), Dimension(:), Allocatable ::  bsc, bscat, masa, mat
  Real (myprec) :: q(npmx),Deltar,dr,cq,dq,cr,tuniti
  Real (myprec), Parameter :: pi=3.1415926535898D0,fconv=4.135,Rgas&
       &=8.3144d7, kelvintokcal=0.00198717,aunit=1.0d-8,tunit=1.0d-12
  Real (myprec), Dimension(3,3) :: cell
  Real (myprec), Dimension(3) :: sidel, rlow, rup
  Real (myprec) :: weight, charge, tstep, side, kmax, volumen, sigma,&
       & pwall, pwallp, volumen_r, side2, bav, bav2, b2av, densty,&
       & bsprod, rdenst
  Real, allocatable :: histomix(:,:,:), gcluster(:), gclustav(:),&
       & rhoclus(:), rhoclusav(:), gclcl(:), sqcl(:), sqf(:), sqfp(:,:), sqfcl(:), densprof(:,:)
  Real ::  fk(3), qmin, qmax, qmin2
  Character :: fname88*18, fname55*18,fname72*18,fname70*18,fname77*18&
       &,fname66*18, fname33*18,fname99*18,fname44*18,fname22*18&
       &,fname23*18,fname11*18,atnam*8,npdim(3)*1=(/'x','y','z'/)
  logical :: pdir(3) = .true.
  Character, Dimension(:), Allocatable :: atoms*8
  Integer, Allocatable :: lty(:), nq(:), ipos(:,:), histomixi(:,:,:), densprofi(:,:)

End Module Comun

module dev_def
  use cudafor
  !
  ! Device code
  !
  real, device, allocatable ::r_d(:,:), rclxyzd(:,:)
  integer, device, allocatable :: itype_d(:), ipos_d(:,:), ntype_d(:)&
       &, nq_d(:), gdone(:), histomix_d(:,:,:), densprof_d(:,:)
  real, device, allocatable :: sqf_d(:), sqfp_d(:,:), sqfcl_d(:)
  real, device ::   fk_d(3), sidel_d(3)
  type (cudaDeviceProp) :: prop
  type(cudaEvent) :: startEvent, stopEvent

end module dev_def

Module Fourier
  use g_types
  !
  ! Transfer parameters for initialization of FFT routines
  !
  Use Comun, Only : pi, npmx
  Integer N , NR , NU
  Real (myprec), Dimension (npmx) :: SINe , COSe , FR , FI
  Interface
     Subroutine FFTSIN(F,C,Tf)
       use g_types
       use comun, only : npmx
       Implicit None
       real (myprec), intent(IN) :: C
       real (myprec), dimension(NPMX), intent(INOUT) :: F
       real (myprec), dimension(NPMX), intent(OUT) :: TF
     end Subroutine FFTSIN
     Subroutine FFTCOS(F0,F,Tf)
       use g_types
       use comun, only : npmx
       real (myprec), intent(IN) :: F0
       real (myprec), dimension(NPMX), intent(INOUT) :: F
       real (myprec), dimension(NPMX), intent(OUT) :: TF
     end Subroutine FFTCOS
  end Interface

End Module Fourier



module clusters
  use g_types
  implicit none
  type clus
     sequence
     real (myprec) :: vl(3)
     real (myprec) :: center(3)
     real (myprec) :: ekin
     real (myprec) :: radio
     real (myprec) :: cldens
     real (myprec) :: mass
     Integer :: clsize
     integer, allocatable :: members(:)
  end type clus
  type (clus), Allocatable :: cluster(:)
  Real (myprec), Dimension(:), Allocatable :: sizedist, densav
  Real (myprec), Dimension(:,:), Allocatable ::  rclxyz
  real (myprec) :: rcl, ekclaver, ekinclsav, dcl, ekincl, ekincls
  real (myprec) ::  tadj=0, tgraph=0, tthrus=0, tbfs=0, avradio=0, averdens=0, drho
  integer, allocatable ::   radii(:), densclus(:), contador(:)
  integer :: jmin, minclsize, Nu_clus, ndrho
  real (myprec) :: NTclus
end module clusters




Function vector_product(a,b)
  use g_types
  Real (myprec), Dimension(:), Intent(IN) :: a, b
  Real (myprec), Dimension(Size(a)) :: vector_product
  vector_product(1) = a(2)*b(3)-a(3)*b(2)
  vector_product(2) = -a(1)*b(3)+a(3)*b(1)
  vector_product(3) = a(1)*b(2)-a(2)*b(1)
End Function vector_product

Function surface(a,b)
  use g_types
  Real (myprec), Dimension(:), Intent(IN) :: a, b
  Real (myprec) :: surface
  surface = Abs(a(1)*b(2)-a(2)*b(1))
End Function surface


function radius(members, nm, center,sidel)
  use comun, only : r
  use g_types
  implicit none
  integer, intent(IN) :: nm
  integer, intent(IN) :: members(nm)
  real(myprec), intent(IN) :: center(3), sidel(3)
  real(myprec) :: rr2, rv(3)
  integer :: i
  real(myprec) :: radius
  radius = 0
  do i = 1, nm
     rv(:)= r(:,members(i))-center(:)
     rv = rv-sidel*nint(rv/sidel)
     rr2= Dot_product(rv,rv)
     if (rr2 > radius) radius = rr2
  end do
  radius = sqrt(radius)
end function radius
  
