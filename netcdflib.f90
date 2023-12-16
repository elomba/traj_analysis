!
! This code reads a LAMMPS/AMBER NetCDF file with a series of
! configurations, extracts all data, and then rewrites the
! configurations in to a new NetCDF4 file to illustrate how to define
! a data frame and dump it to a destination file. Specific atom types
! can be selected at run time
!
! Programmed by Enrique Lomba, IQF, Madrid, November 2023
!
! IMPORTANT NOTE !!! LAMMPS netcdf dump format specifies that the
! cell length of non pediodic dimensions is set to 0, irrespective of
! the information set in the lammps script. To cope with this the
! program redifnes the length using abs(max(r_xyz)-min(r_xyz))+10 or
! 2*abs(origin_xyz) depending on the dimension involved. 
!
! to be compiled with
! ifort -I /usr/local/intel_netcdf/include \
!           -L/usr/local/intel_netcdf/lib64 read_nc.f90 -lnetcdff
! export LD_LIBRARY_PATH=/usr/local/intel_netcdf/lib64:$LD_LIBRARY_PATH
!


module myncdf
  ! 
  use g_types
  interface 
     subroutine read_nc_cfg(ncid,ncstart,io,unit)
       integer, intent(in) :: ncid, ncstart
       integer, intent (in), optional :: unit
       integer, intent(out) :: io
     end subroutine read_nc_cfg
     subroutine write_nc_cfg(ncid,ncstart)
       integer, intent(in) :: ncid, ncstart
     end subroutine write_nc_cfg
     subroutine write_nc_cfg_red(ncid,ncstart,ntms)
       integer, intent(in) :: ncid, ncstart,  ntms
     end subroutine write_nc_cfg_red
     SUBROUTINE check(istatus,ioerr)
       INTEGER, INTENT (IN) :: istatus
       INTEGER, INTENT (OUT) :: ioerr
     end SUBROUTINE check
  end interface
  ! Define a new type to hold all details of a given configuration variable
  type config
     character (len=:), allocatable :: varname
     integer :: varid
     integer :: numdims
     integer :: numattrs
     integer :: xtype
     integer, dimension(:), allocatable :: dimlen, dimids
     real(myprec) :: scale
     character (len=:), allocatable :: units
  end type config
  ! Character variables of flexible length
  type globalat
     character (len=:), allocatable :: name
     character (len=:), allocatable :: val
  end type globalat
  type dimens
     integer :: length
     character (len=:), allocatable :: name
  end type dimens
  ! Remap NetCDF F90 variable types
  character (len=*), parameter, dimension(6) :: tipos =(/"NF90_BYTE  ","&
       &NF90_CHAR  ","NF90_SHORT ","NF90_INT   ","NF90_FLOAT ","NF90_D&
       &OUBLE"/)
contains
  ! Routine to return NetCDF error codes
  SUBROUTINE check(istatus,ioerr)
    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: istatus
    INTEGER, INTENT (OUT) :: ioerr
    IF (istatus /= nf90_noerr) THEN
       write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
       ioerr = istatus
    else
       ioerr = nf90_noerr
    END IF
  END SUBROUTINE check
end module myncdf
!
! Define general variables to store configuration vales 
! 
module configuration
  use myncdf
  use netcdf
  type (config), dimension(:), allocatable :: conf
  type (dimens), dimension(:), allocatable :: mydims
  type (globalat), dimension(:), allocatable :: myglobatts
  real(myprec), dimension(:,:,:), allocatable :: r, v
  real(myprec), dimension(3,1) :: org, cell, cell_a
  real(myprec) :: time, scale
  integer, dimension(:,:), allocatable :: idi, ity, imol
  character (len = NF90_MAX_NAME) :: nombre, nom_att
  character (len = 3) :: spatial, cell_spatial
  character (len = 5), dimension(1:3) :: cell_ang
  character (len = 55) :: fname, fnamew, attr
  integer, dimension(nf90_max_var_dims) :: rhDimIds, rdimId
  integer :: vtype, nDims, nVars, nGlobalAtts, unlimDimID, step(1), ntypes, nmconf
  integer :: start(3), count(3), starti(2), counti(2),&
       & cstart(2), ccount(2), csstart(1),cscount(1), ln, ntm&
       &, status, varid, numatts, numdims, natoms, nconf, nlen, nwty
  integer, dimension(:), allocatable :: ndimid, nvarid, atypes, wtypes
  logical :: ex_vel=.false., ex_mol=.false., periodic(3)=.true.
end module configuration

subroutine write_nc_cfg(ncid,ncstart)
  use myncdf
  use netcdf
  use configuration
  implicit none
  integer, intent(in) :: ncid, ncstart
  logical, save :: first=.true.
  integer :: i, j , k, tipo, ioerr
  ! 3 spatial coordinate
  count(1) = 3
  ! read in frame after frame 
  count(3) = 1
  !
  starti(:) = 1
  starti(2) = ncstart
  counti(:) = 1
  start(1) = 1
  start(2) = 1
  start(3) = ncstart
  cstart(1) = 1
  cstart(2) = ncstart
  ccount(1) = 3
  ccount(2) = 1
  if (first) then
     ! Write global attributes
     do i = 1, nglobalatts
        call check(nf90_put_att(ncid,NF90_GLOBAL,myglobatts(i)&
             &%name,myglobatts(i)%val),ioerr)
     end do
     allocate(ndimid(ndims))
     !
     ! define dimensions ( "frame" is defined as unilimted, it will
     ! store configurations as they are written)
     ! nf90_def_dim(ncid,name,lenght_of_dim,dimension_id), dimension_id=output
     !
     do i = 1, ndims
        if (mydims(i)%name == "frame") mydims(i)%length = NF90_UNLIMITED
        call check(nf90_def_dim(ncid,mydims(i)%name,mydims(i)%length&
             &,ndimid(i)),ioerr)
     enddo
     allocate(nvarid(nvars))

     do i=1, nvars
        !
        ! define variables,
        !
        ! nf90_def_var(ncid,name,NF90_TYPE,array_with_dimIDs(:),varid), varid (integer
        !             identifying variable, output)
        !
        call check(nf90_def_var(ncid,conf(i)%varname,conf(i)%xtype,conf(i)&
             &%dimids,nvarid(i)),ioerr)
     enddo
     ! Exite define mode
     call check(nf90_enddef(ncid),ioerr)
  end if
  !
  !   Store variable values, start=starting position, count=number of elements to transfer.
  !
  do i = 1, nvars
     j = nvarid(i)
     do k = 1, conf(j)%numattrs
        call check(nf90_put_att(ncid,j,"units",conf(j)%units),ioerr)
        if (conf(j)%varname == "time") then
           call check(nf90_put_att(ncid,j,"scale_factor",conf(j)%scale),ioerr)
        endif
     end do
     if (conf(j)%varname == "coordinates") then
        start(3) = ncstart
        count(2) = natoms
        call check(nf90_put_var(ncid,j,r,start,count),ioerr)
     else if (conf(j)%varname == "velocities") then
        start(3) = ncstart
        count(2) = natoms
        status =nf90_put_var(ncid,j,v,start,count)        
     else if (conf(j)%varname == "type") then
        counti(1) = conf(j)%dimlen(1)
        status =nf90_put_var(ncid,j,ity,starti,counti)
     else if (conf(j)%varname == "id") then
        counti(1) =  conf(j)%dimlen(1)
        status =nf90_put_var(ncid,j,idi,starti,counti)
     else if (conf(i)%varname == "mol") then
        counti(1) = conf(i)%dimlen(1)
        status =nf90_put_var(ncid,i,imol,starti,counti)
     else if (conf(j)%varname == "cell_origin") then
        ccount(2) = 1
        status =nf90_put_var(ncid,j,org,cstart,ccount)  
     else if (conf(j)%varname == "cell_lengths") then
        ccount(2) = 1
        status =nf90_put_var(ncid,j,cell,cstart,ccount)
     else if (conf(j)%varname == "cell_angles") then
        ccount(2) = 1
        status =nf90_put_var(ncid,j,cell_a,cstart,ccount)
     else if (conf(j)%varname == "cell_spatial") then
        if (first) status = nf90_put_var(ncid,j,cell_spatial)
     else if (conf(j)%varname == "cell_angular") then
        if (first) status = nf90_put_var(ncid,j,cell_ang)
     else if (conf(j)%varname == "spatial") then
        if (first) status = nf90_put_var(ncid,j,spatial(1:3))
     else if (conf(j)%varname == "time") then   
        csstart(1) = ncstart
        cscount(1) = 1
        status = nf90_put_var(ncid,j,step, csstart, cscount)
     end if

  end do
  first = .false.
end subroutine write_nc_cfg

subroutine read_nc_cfg(ncid,ncstart,io,unit)
  use myncdf
  use netcdf
  use configuration
  implicit none
  integer, intent(in) :: ncid, ncstart
  integer, intent (in), optional :: unit
  integer, intent(out ) :: io
  logical, save :: first=.true., typedefined=.false. 
  integer :: i, j , tipo, iunit, k, ioerr
  if (present(unit)) then
     iunit = unit
  else
     iunit = 6
  endif
  ! 3 spatial coordinate
  count(1) = 3
  ! read in frame after frame 
  count(3) = 1
  !
  starti(:) = 1
  starti(2) = ncstart
  counti(:) = 1
  start(1) = 1
  start(2) = 1
  start(3) = ncstart
  cstart(1) = 1
  cstart(2) = ncstart
  ccount(1) = 3
  ccount(2) = 1

  if (first) then
     if (iunit .ne. 6) then
        Write(*,"(//'**** Using NetCDF library version ',A/'*** output to netcdf.log')") trim(nf90_inq_libvers())
        open(iunit,file="netcdf.log")
     endif
     Write(iunit,"(//'**** Using NetCDF library version ',A)") trim(nf90_inq_libvers())
     !
     ! Inquire dimensions and number of vars
     !
     call check(nf90_inquire(ncid, ndims, nvars, nglobalatts,&
          & unlimdimid),ioerr)
     write(iunit,"(// ' No. of dimensions ',i3,' No of vars ', i3,' No. of g&
          &lobal attibutes ', i3, ' Unlimited dimension ', i3)") ndims, nvars, nglobalatts,&
          & unlimdimid
     if (nglobalatts>0) allocate(myglobatts(nglobalatts))
     do j = 1, nglobalatts
        call check(nf90_inq_attname(ncid,NF90_GLOBAL,j, nom_att),ioerr)
        nlen = len(trim(adjustl(nom_att)))
        allocate(character (len=nlen) :: myglobatts(j)%name)
        myglobatts(j)%name = trim(adjustl(nom_att))
        call check(nf90_inquire_attribute(ncid,NF90_GLOBAL,  myglobatts(j)%name,xtype&
             &=tipo, len=nlen),ioerr)
        allocate(character (len=nlen) :: myglobatts(j)%val)
        call check(nf90_get_att(ncid,NF90_GLOBAL, nom_att, attr),ioerr)
        myglobatts(j)%val = trim(adjustl(attr))
        write(iunit,"(' - Global Attribute: ',2A25,' length',i3,' valor: ',A)")nom_att&
             &,tipos(tipo),nlen,  myglobatts(j)%val 
     end do
     allocate(conf(nvars))
     allocate(mydims(ndims))
     write(iunit,"(/' ** Readings dim id, name, length ',/)")
     do j = 1, ndims
        call check(nf90_inquire_dimension(ncid, j, name = nombre, len =&
             & nlen),ioerr)
        mydims(j)%length = nlen
        ln = len( trim(adjustl(nombre)))
        allocate(character(len=ln) :: mydims(j)%name)
        mydims(j)%name = trim(adjustl(nombre))
        write(iunit,"( 'dimension id',i2,' nombre ',A, ' length=',i)")j,mydims(j)%name,mydims(j)%length
        if ( mydims(j)%name == "frame" ) nmconf = mydims(j)%length
     enddo
     write(iunit,"(/' ** Readings vars, type, no. dimensions, no. a&
          &ttrs ',/)")
  endif
  if (ncstart > nmconf) then
     write(*,"(' ** Warning : trying to read past last configuration no. ',i7 )") nmconf
     io = -1
     return
  end if
  do i=1, nvars
     if (first) then
        call check(nf90_inquire_variable(ncid, i, xtype=tipo, name = nombre, ndims =&
             & numdims, natts = numatts),ioerr)
        conf(i)%xtype = tipo
        conf(i)%varid = i
        write(iunit, "('varid=',i2,' ',2A15,' ndims =',i2,' nattrs=', i2)")i&
             &,nombre, tipos(tipo), numdims, numatts
        ln = len( trim(adjustl(nombre)))
        allocate(character(len=ln) :: conf(i)%varname)
        conf(i)%varname(1:ln) =  trim(adjustl(nombre))
        conf(i)%numdims = numdims
        conf(i)%numattrs = numatts
        allocate(conf(i)%dimlen(numdims))
        call check(nf90_inquire_variable(ncid, i, dimids =&
             & rhdimids(1:numdims)),ioerr)
        allocate(conf(i)%dimids(numdims))
        conf(i)%dimids(1:numdims)=rhdimids(1:numdims)
        do j = 1, numdims
           call check(nf90_inquire_dimension(ncid, rhdimids(j), len =&
                & nlen),ioerr)
           conf(i)%dimlen(j) = nlen
           write(iunit,"('  ** dim id =',i3,' length = ',i6)")rhdimids(j)&
                &, nlen
        end do
        do j = 1, numatts
           call check(nf90_inq_attname(ncid,i,j, nom_att),ioerr)
           call check(nf90_inquire_attribute(ncid,i, nom_att,xtype&
                &=tipo, len=nlen),ioerr)
           write(iunit,"('    · Attribute ',2A15,' length',i3)")nom_att&
                &,tipos(tipo),nlen
        enddo
        if (numatts > 0 ) then
           call check(nf90_inq_attname(ncid,i,1, nom_att),ioerr)
           call check(nf90_get_att(ncid,i, nom_att, attr),ioerr)
           ln = len(trim(adjustl(attr)))
           allocate(character(len=ln) :: conf(i)%units )
           conf(i)%units = trim(adjustl(attr))
        endif
        if (numatts == 2) then
           call check(nf90_inq_attname(ncid,i,2, nom_att),ioerr)
           call check(nf90_get_att(ncid,i, nom_att, scale),ioerr)
           conf(i)%scale = scale
        else if  (numatts > 2) then
           write(*,"(' ** Error: program cannot handles more than two attr&
                &ibutes for var ',A)") nombre
           stop
        end if
     end if

     
     
     if (conf(i)%varname == "coordinates") then
        if (first) allocate(r( conf(i)%dimlen(1), conf(i)%dimlen(2),1))
        natoms = conf(i)%dimlen(2)
        nconf = conf(i)%dimlen(3)
        start(3) = ncstart
        count(2) = natoms
        call check(nf90_get_var(ncid,i,r,start,count),ioerr)
        do k=1, 3
           if (cell(k,1) == 0) then
              periodic(k) = .false.
              cell(k,1) = max(abs(maxval(r(k,1:natoms,1))-minval(r(k&
                   &,1:natoms,1)))+10.0,2*abs(org(k,1)))
           end if
        end do
     else if (conf(i)%varname == "velocities") then
        ex_vel = .true.
        if (first) allocate(v( conf(i)%dimlen(1), conf(i)%dimlen(2),1))
        start(3) = ncstart
        count(2) = natoms
        call check(nf90_get_var(ncid,i,v,start,count),ioerr)
     else if (conf(i)%varname == "type") then
        typedefined = .true.
        if (first) allocate(ity( conf(i)%dimlen(1),1))
        counti(1) = conf(i)%dimlen(1)
        call check(nf90_get_var(ncid,i,ity,starti,counti),ioerr)
     else if (conf(i)%varname == "mol") then
        ex_mol=.true.
        if (first) allocate(imol( conf(i)%dimlen(1),1))
        counti(1) = conf(i)%dimlen(1)
        call check(nf90_get_var(ncid,i,imol,starti,counti),ioerr)
     else if (conf(i)%varname == "id") then
        if (first) allocate(idi( conf(i)%dimlen(1),1))
        counti(1) =  conf(i)%dimlen(1)
        ! Reasign id in a contiguos fashion
        do j=1, counti(1)
           idi(j,1) = j
        end do
        call check(nf90_get_var(ncid,i,idi,starti,counti),ioerr)
     else if (conf(i)%varname == "cell_origin") then
        ccount(2) = 1
        call check(nf90_get_var(ncid,i,org,cstart,ccount),ioerr)
     else if (conf(i)%varname == "cell_lengths") then
        ccount(2) = 1
        call check(nf90_get_var(ncid,i,cell,cstart,ccount),ioerr)
     else if (conf(i)%varname == "cell_angles") then
        ccount(2) = 1
        call check(nf90_get_var(ncid,i,cell_a,cstart,ccount),ioerr)
     else if (conf(i)%varname == "cell_spatial") then
        call check(nf90_get_var(ncid,i,cell_spatial),ioerr)
     else if (conf(i)%varname == "cell_angular") then
        call check(nf90_get_var(ncid,i,cell_ang),ioerr)
     else if (conf(i)%varname == "spatial") then
        call check(nf90_get_var(ncid,i,spatial(1:3)),ioerr)
     else if (conf(i)%varname == "time") then   
        csstart(1) = ncstart
        cscount(1) = 1
        call check(nf90_get_var(ncid,i,step, csstart, cscount),ioerr)
     end if
     if (ioerr /= nf90_noerr) then
        ioerr = -1
        return
     end if
  enddo
  if (typedefined) then
     if (first) then
        ntypes = maxval(ity)
        write(*, "(/' ** Number of atoms types =',i3)") ntypes
        allocate(atypes(ntypes))
        atypes(:) = 0
        do i=1, natoms
           atypes(ity(i,1)) = atypes(ity(i,1))+1
        end do
        do i = 1, ntypes
           write(*, "(/'       Number of atoms of types ',i2,'=',i8)") i, atypes(i)
        end do
     endif
  else
     print *, " **** Warning no atom types defined in trajectory file"
  endif
  if (iunit.ne. 6 .and. first) close(iunit)
  first = .false.

end subroutine read_nc_cfg


subroutine write_nc_cfg_red(ncid,ncstart,ntms)
  use myncdf
  use netcdf
  use configuration
  implicit none
  real(myprec), dimension(:,:,:), allocatable, save :: rn, vn
  integer, dimension(:,:), allocatable, save :: idin, ityn, imoln
  integer, intent(in) :: ncid, ncstart,  ntms
  logical, save :: first=.true.
  integer :: i, j , k, tipo(1), index, ioerr
  if (first) then
     allocate(rn(3,ntms,1),ityn(ntms,1),idin(ntms,1))
     if (ex_vel) allocate(vn(3,ntms,1))
     if (ex_mol) allocate(imoln(ntms,1))
     index = 0
     do i = 1, natoms
        if (any(wtypes == ity(i,1))) then
           index = index + 1
           tipo = findloc(wtypes,ity(i,1))
           ityn(index,1) = tipo(1)
           idin(index,1) = index
        end if
     end do
     if (index .ne. ntms) then
        print *, " *** Error in selected no. of atoms "
        stop
     end if
  endif
  index = 0
  do i = 1, natoms
     if (any(wtypes == ity(i,1))) then
        index = index + 1
        tipo(:) = findloc(wtypes,ity(i,1))
        rn(1:3,index,1) = r(1:3,i,1)
        if (ex_vel) vn(1:3,index,1) = v(1:3,i,1)
        if (ex_mol) imoln(index,1) = imol(i,1)
     end if
  end do
  ! 3 spatial coordinate
  count(1) = 3
  ! read in frame after frame 
  count(3) = 1
  !
  starti(:) = 1
  starti(2) = ncstart
  counti(:) = 1
  start(1) = 1
  start(2) = 1
  start(3) = ncstart
  cstart(1) = 1
  cstart(2) = ncstart
  ccount(1) = 3
  ccount(2) = 1
  if (first) then
     ! Write global attributes
     do i = 1, nglobalatts
        call check(nf90_put_att(ncid,NF90_GLOBAL,myglobatts(i)&
             &%name,myglobatts(i)%val),ioerr)
     end do
     allocate(ndimid(ndims))
     !
     ! define dimensions ( "frame" is defined as unlimted, it will
     ! store configurations as they are written)
     ! nf90_def_dim(ncid,name,lenght_of_dim,dimension_id), dimension_id=output
     !
     do i = 1, ndims
        if (mydims(i)%name == "frame") then
           mydims(i)%length = NF90_UNLIMITED
        else if (mydims(i)%name == "atom") then

           mydims(i)%length = ntms
        end if
        call check(nf90_def_dim(ncid,mydims(i)%name,mydims(i)%length&
             &,ndimid(i)),ioerr)
     enddo
     allocate(nvarid(nvars))

     do i=1, nvars
        !
        ! define variables,
        !
        ! nf90_def_var(ncid,name,NF90_TYPE,array_with_dimIDs(:),varid), varid (integer
        !             identifying variable, output)
        !
        call check(nf90_def_var(ncid,conf(i)%varname,conf(i)%xtype,conf(i)&
             &%dimids,nvarid(i)),ioerr)
     enddo
     ! Exite define mode
     call check(nf90_enddef(ncid),ioerr)
  end if
  !
  !   Store variable values, start=starting position, count=number of elements to transfer.
  !
  do i = 1, nvars
     j = nvarid(i)
     do k = 1, conf(j)%numattrs
        call check(nf90_put_att(ncid,j,"units",conf(j)%units),ioerr)
        if (conf(j)%varname == "time") then
           call check(nf90_put_att(ncid,j,"scale_factor",conf(j)%scale),ioerr)
        endif
     end do
     if (conf(j)%varname == "coordinates") then
        start(3) = ncstart
        count(2) = ntms
        call check(nf90_put_var(ncid,j,rn,start,count),ioerr)
     else if (conf(j)%varname == "velocities") then
        start(3) = ncstart
        count(2) = ntms
        call check(nf90_put_var(ncid,j,vn,start,count),ioerr)
     else if (conf(j)%varname == "type") then
        counti(1) = ntms
        call check(nf90_put_var(ncid,j,ityn,starti,counti),ioerr)
     else if (conf(j)%varname == "id") then
        counti(1) =  ntms
        call check(nf90_put_var(ncid,j,idin,starti,counti),ioerr)
     else if (conf(i)%varname == "mol") then
        counti(1) = ntms
        call check(nf90_put_var(ncid,i,imoln,starti,counti),ioerr)
     else if (conf(j)%varname == "cell_origin") then
        ccount(2) = 1
        call check(nf90_put_var(ncid,j,org,cstart,ccount),ioerr)
     else if (conf(j)%varname == "cell_lengths") then
        ccount(2) = 1
        call check(nf90_put_var(ncid,j,cell,cstart,ccount),ioerr)
     else if (conf(j)%varname == "cell_angles") then
        ccount(2) = 1
        call check(nf90_put_var(ncid,j,cell_a,cstart,ccount),ioerr)
     else if (conf(j)%varname == "cell_spatial") then
        if (first) call check( nf90_put_var(ncid,j,cell_spatial),ioerr)
     else if (conf(j)%varname == "cell_angular") then
        if (first) call check( nf90_put_var(ncid,j,cell_ang),ioerr)
     else if (conf(j)%varname == "spatial") then
        if (first) call check( nf90_put_var(ncid,j,spatial(1:3)),ioerr)
     else if (conf(j)%varname == "time") then   
        csstart(1) = ncstart
        cscount(1) = 1
        call check( nf90_put_var(ncid,j,step, csstart, cscount),ioerr)
     end if

  end do
  first = .false.
end subroutine write_nc_cfg_red
