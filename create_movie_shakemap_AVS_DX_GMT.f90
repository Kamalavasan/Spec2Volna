!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!
!---  create a movie of the vertical component of surface displacement or velocity
!---  or a ShakeMap(R) (i.e. map of the maximum absolute value of the two horizontal components
!---  of the velocity vector) in AVS, OpenDX or GMT format
!
!
! AVS UCD format descriptions:
! https://lanl.github.io/LaGriT/pages/docs/read_avs.html
! http://people.sc.fsu.edu/~jburkardt/data/ucd/ucd.html
! http://www.hnware.de/rismo/dokumente/anwenderdoku/formate/avs_ucd.html

  program create_movie_shakemap

  use constants
  use shared_parameters

  implicit none

!-------------------------------------------------------------------------------------------------
! user parameters

! normalizes field display values
  logical, parameter :: NORMALIZE_OUTPUT = .false.

! threshold in percent of the maximum below which we cut the amplitude
  logical, parameter :: APPLY_THRESHOLD = .false.
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD = 1._CUSTOM_REAL / 100._CUSTOM_REAL

! coefficient of power law used for non linear scaling
  logical, parameter :: NONLINEAR_SCALING = .false.
  real(kind=CUSTOM_REAL), parameter :: POWER_SCALING = 0.13_CUSTOM_REAL

! muting source region
  logical, parameter :: MUTE_SOURCE = .false.
  real(kind=CUSTOM_REAL), parameter :: RADIUS_TO_MUTE = 1000._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: X_SOURCE_EXT_MESH = -9023.021484375
  real(kind=CUSTOM_REAL), parameter :: Y_SOURCE_EXT_MESH = 6123.611328125
  real(kind=CUSTOM_REAL), parameter :: Z_SOURCE_EXT_MESH = 17.96331405639648

!-------------------------------------------------------------------------------------------------

  integer :: it,it1,it2,ivalue,nspectot_AVS_max,ispec
  integer :: nframes,iframe,inumber,inorm,iscaling_shake
  ! integer :: ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4

  logical :: plot_shaking_map

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,display
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: vectorx,vectory,vectorz,vectornorm

  double precision :: min_field_current,max_field_current,max_absol

  character(len=MAX_STRING_LEN) :: outputname
  character(len=MAX_STRING_LEN) :: line

  integer :: ipoin

  ! GMT
  ! double precision :: lat,long

  ! for sorting routine
  integer :: npointot,ilocnum,nglob,i,j,ielm,ieoff,ispecloc
  integer, dimension(:), allocatable :: iglob,locval,ireorder
  logical, dimension(:), allocatable :: ifseg,mask_point
  double precision, dimension(:), allocatable :: xp,yp,zp,xp_save,yp_save,zp_save,field_display

  ! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz

  ! copy array for sorting
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: all_data
  real(kind=CUSTOM_REAL):: pt(6,8)
  integer:: i_v, vmeshpts
  ! count_pt
  ! logical::pt_valid
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable:: coords
  integer :: ier

  ! order of points representing the 2D square element
  integer,dimension(NGNOD2D_FOUR_CORNERS_AVS_DX),parameter :: iorder = (/1,3,2,4/)

  integer :: NSPEC_SURFACE_EXT_MESH
  logical :: BROADCAST_AFTER_READ

  ! time measurement
  real :: start, finish

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all movie frames to create a movie'
  print *

  print *
  print *,'reading parameter file'
  print *

  ! initializes
  myrank = 0
  BROADCAST_AFTER_READ = .false.

  ! read the parameter file
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! only one global array for movie data, but stored for all surfaces defined
  ! in file 'surface_from_mesher.h'
  open(unit=IIN,file=trim(OUTPUT_FILES)//'surface_from_mesher.h',status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(OUTPUT_FILES)//'surface_from_mesher.h'
    print *
    print *,'please run xgenerate_databases or xspecfem3D first to create this file, exiting now...'
    stop 'error opening moviedata header file'
  endif
  ! skips first few lines
  do i=1,6
    read(IIN,'(a)') line
  enddo
  ! line with info, e.g. "integer,parameter :: NSPEC_SURFACE_EXT_MESH = 23855"
  read(IIN,'(a)') line
  close(IIN)
  ! gets number from substring after = sign
  i = index(line,'=')
  if (i == 0) stop 'error reading in NSPEC_SURFACE_EXT_MESH from file OUTPUT_FILES/surface_from_mesher.h'

  read(line(i+1:len_trim(line)),*) NSPEC_SURFACE_EXT_MESH

  ! calculates number of total surface points
  if (USE_HIGHRES_FOR_MOVIES) then
     ilocnum = NSPEC_SURFACE_EXT_MESH*NGLLSQUARE
  else
     ilocnum = NSPEC_SURFACE_EXT_MESH*NGNOD2D_FOUR_CORNERS_AVS_DX
  endif
  print *,'  high-resolution: ',USE_HIGHRES_FOR_MOVIES
  print *,'  moviedata element surfaces: ',NSPEC_SURFACE_EXT_MESH
  print *,'  moviedata total elements all: ',ilocnum
  print *

  if (SAVE_DISPLACEMENT) then
    print *,'Displacement will be shown in movie'
  else
    print *,'Velocity will be shown in movie'
  endif
  print *

  if (MUTE_SOURCE) then
    print *,'Muting source region:'
    print *,'  radius = ',RADIUS_TO_MUTE
    print *,'  source location x/y/z = ',X_SOURCE_EXT_MESH,Y_SOURCE_EXT_MESH,Z_SOURCE_EXT_MESH
    print *
  endif

  ! user input

  plot_shaking_map = .false.
  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  ! saving only the last time step 
  it1 = (NSTEP/NTSTEP_BETWEEN_FRAMES)*NTSTEP_BETWEEN_FRAMES
  it2 = it1
  ! frame number for the output file
  inumber = 1
  if (.not. plot_shaking_map) then

    ! limits to maximum of NSTEP
    if (it2 > NSTEP) then
      it2 = NSTEP
    endif

    ! saving only last time step
    print *
    print *,'looping from ',it1,' to ',it2,' every ',NTSTEP_BETWEEN_FRAMES,' time steps'
    ! count number of movie frames
    nframes = 0
    do it = it1,it2
      if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) nframes = nframes + 1
    enddo
  else
    ! only one frame if shaking map
    nframes = 1
    it1 = 1
    it2 = 1
  endif
  print *
  print *,'total number of frames will be ',nframes
  if (nframes == 0) stop 'null number of frames'

  iscaling_shake = 0
  print *, 'norm of displacement will be saved in bathy files'
  inorm = 1
  


  ! define the total number of elements at the surface
  if (USE_HIGHRES_FOR_MOVIES) then
     nspectot_AVS_max = NSPEC_SURFACE_EXT_MESH * (NGLLX-1) * (NGLLY-1)
  else
     nspectot_AVS_max = NSPEC_SURFACE_EXT_MESH
  endif

  ! maximum theoretical number of points at the surface
  npointot = NGNOD2D_FOUR_CORNERS_AVS_DX * nspectot_AVS_max

  ! allocate arrays for sorting routine
  allocate(iglob(npointot), &
           locval(npointot), &
           ifseg(npointot), &
           xp(npointot), &
           yp(npointot), &
           zp(npointot), &
           xp_save(npointot), &
           yp_save(npointot), &
           zp_save(npointot), &
           field_display(npointot), &
           mask_point(npointot), &
           ireorder(npointot),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for sorting routine'

  ! allocates data arrays
  allocate(store_val_x(ilocnum), &
           store_val_y(ilocnum), &
           store_val_z(ilocnum), &
           store_val_ux(ilocnum), &
           store_val_uy(ilocnum), &
           store_val_uz(ilocnum),stat=ier)

  if (ier /= 0) stop 'Error allocating arrays for data arrays'

  allocate(all_data(6,npointot), stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for all data arrays'


  if (USE_HIGHRES_FOR_MOVIES) then
    allocate(x(NGLLX,NGLLY), &
             y(NGLLX,NGLLY), &
             z(NGLLX,NGLLY), &
             display(NGLLX,NGLLY),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for highres'
  endif

  ! user output
  print *
  print *,'there are a total of ',nspectot_AVS_max,' elements at the surface'
  print *
  print *
  if (APPLY_THRESHOLD .and. .not. plot_shaking_map) &
    print *,'Will apply a threshold to amplitude below ',100.*THRESHOLD,' %'
  if (NONLINEAR_SCALING .and. (.not. plot_shaking_map .or. iscaling_shake == 1)) &
    print *,'Will apply a non linear scaling with coef ',POWER_SCALING

  iframe = 0

  ! reading volna mesh points
  call mesh2vec(coords, vmeshpts)
  print *, 'vmeshpts: ', vmeshpts


! loop on all the time steps in the range entered
  do it = it1,it2

    ! check if time step corresponds to a movie frame
    if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0 .or. plot_shaking_map) then

      iframe = iframe + 1

      print *
      if (plot_shaking_map) then
        print *,'reading shaking map snapshot'
      else
        print *,'reading snapshot time step ',it,' out of ',NSTEP
      endif
      print *

      ! read all the elements from the same file
      if (plot_shaking_map) then
        write(outputname,"('/shakingdata')")
      else
        write(outputname,"('/moviedata',i6.6)") it
      endif
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname),status='old', &
            action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'error: ',trim(OUTPUT_FILES)//trim(outputname)
        stop 'error opening moviedata file'
      endif

      read(IOUT) store_val_x
      read(IOUT) store_val_y
      read(IOUT) store_val_z
      read(IOUT) store_val_ux
      read(IOUT) store_val_uy
      read(IOUT) store_val_uz
      close(IOUT)

      ! clear number of elements kept
      ispec = 0

      ! reset point number
      ipoin = 0

      do ispecloc = 1, NSPEC_SURFACE_EXT_MESH

        if (USE_HIGHRES_FOR_MOVIES) then
          ! assign the OpenDX "elements"
          do j = 1,NGLLY
            do i = 1,NGLLX
              ipoin = ipoin + 1

              ! x,y,z coordinates
              xcoord = store_val_x(ipoin)
              ycoord = store_val_y(ipoin)
              zcoord = store_val_z(ipoin)

              ! note:
              ! for shakemaps: ux = norm displacement, uy = norm velocity, uz = norm acceleration
              ! for movies: ux = velocity x-component, uy = velocity y-component, uz = velocity z-component
              vectorx = store_val_ux(ipoin)
              vectory = store_val_uy(ipoin)
              vectorz = store_val_uz(ipoin)
              ! print *, all_data(1,ipoin), all_data(2,ipoin), all_data(3,ipoin)

              x(i,j) = xcoord
              y(i,j) = ycoord
              z(i,j) = zcoord

              ! shakemap
              if (plot_shaking_map) then
                ! chooses norm
                if (inorm == 1) then
                  ! norm displacement
                  display(i,j) = vectorx
                else if (inorm == 2) then
                  ! norm velocity
                  display(i,j) = vectory
                else
                  ! norm acceleration
                  display(i,j) = vectorz
                endif
                !!!! NL NL mute value near source
                if (MUTE_SOURCE) then
                  if ( (sqrt(((x(i,j) - (X_SOURCE_EXT_MESH))**2 + &
                              (y(i,j) - (Y_SOURCE_EXT_MESH))**2 + &
                              (z(i,j) - (Z_SOURCE_EXT_MESH))**2)) < RADIUS_TO_MUTE) ) then
                    display(i,j) = 0.0
                  endif
                endif
              else
                ! movie
                if (inorm == 1) then
                  ! norm of velocity
                  vectornorm = sqrt(vectorz*vectorz + vectory*vectory + vectorx*vectorx)
                  display(i,j) = vectornorm
                else if (inorm == 2) then
                  ! velocity x-component
                  display(i,j) = vectorx
                else if (inorm == 3) then
                  ! velocity y-component
                  display(i,j) = vectory
                else if (inorm == 4) then
                  ! velocity z-component
                  display(i,j) = vectorz
                endif
              endif

            enddo
          enddo

          ! assign the values of the corners of the OpenDX "elements"
          ispec = ispec + 1
          ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)

          do j = 1,NGLLY-1
            do i = 1,NGLLX-1
              ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ielm+(i-1)+(j-1)*(NGLLX-1))
              do ilocnum = 1,NGNOD2D_FOUR_CORNERS_AVS_DX

                if (ilocnum == 1) then
                  xp(ieoff+ilocnum) = dble(x(i,j))
                  yp(ieoff+ilocnum) = dble(y(i,j))
                  zp(ieoff+ilocnum) = dble(z(i,j))
                  field_display(ieoff+ilocnum) = dble(display(i,j))
                else if (ilocnum == 2) then

                  ! accounts for different ordering of square points
                  xp(ieoff+ilocnum) = dble(x(i+1,j+1))
                  yp(ieoff+ilocnum) = dble(y(i+1,j+1))
                  zp(ieoff+ilocnum) = dble(z(i+1,j+1))
                  field_display(ieoff+ilocnum) = dble(display(i+1,j+1))

                else if (ilocnum == 3) then

                  ! accounts for different ordering of square points
                  xp(ieoff+ilocnum) = dble(x(i+1,j))
                  yp(ieoff+ilocnum) = dble(y(i+1,j))
                  zp(ieoff+ilocnum) = dble(z(i+1,j))
                  field_display(ieoff+ilocnum) = dble(display(i+1,j))

                else
                  xp(ieoff+ilocnum) = dble(x(i,j+1))
                  yp(ieoff+ilocnum) = dble(y(i,j+1))
                  zp(ieoff+ilocnum) = dble(z(i,j+1))
                  field_display(ieoff+ilocnum) = dble(display(i,j+1))
                endif

              enddo

            enddo
          enddo

        else
          ! low-resolution (only spectral element corners)
          ispec = ispec + 1
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)

          ! four points for each element
          do i = 1,NGNOD2D_FOUR_CORNERS_AVS_DX

            ! accounts for different ordering of square points
            ilocnum = iorder(i)

            ipoin = ipoin + 1

            xcoord = store_val_x(ipoin)
            ycoord = store_val_y(ipoin)
            zcoord = store_val_z(ipoin)

            ! print *, all_data(1,ipoin), all_data(2,ipoin), all_data(3,ipoin)

            vectorx = store_val_ux(ipoin)
            vectory = store_val_uy(ipoin)
            vectorz = store_val_uz(ipoin)


            xp(ilocnum+ieoff) = dble(xcoord)
            yp(ilocnum+ieoff) = dble(ycoord)
            zp(ilocnum+ieoff) = dble(zcoord)

            ! shakemap
            if (plot_shaking_map) then
              if (inorm == 1) then
                ! norm of displacement
                field_display(ilocnum+ieoff) = dble(vectorx)
              else if (inorm == 2) then
                ! norm of velocity
                field_display(ilocnum+ieoff) = dble(vectory)
              else
                ! norm of acceleration
                field_display(ilocnum+ieoff) = dble(vectorz)
              endif
              !!!! NL NL mute value near source
              if (MUTE_SOURCE) then
                if (sqrt(((dble(xcoord) - (X_SOURCE_EXT_MESH))**2 + &
                          (dble(ycoord) - (Y_SOURCE_EXT_MESH))**2 + &
                          (dble(zcoord) - (Z_SOURCE_EXT_MESH))**2)) < RADIUS_TO_MUTE) then
                  field_display(ilocnum+ieoff) = 0.0
                endif
              endif
            else
              ! movie
              if (inorm == 1) then
                ! norm of velocity
                field_display(ilocnum+ieoff) = sqrt(vectorz**2+vectory**2+vectorx**2)
              else if (inorm == 2) then
                ! velocity x-component
                field_display(ilocnum+ieoff) = vectorx
              else if (inorm == 3) then
                ! velocity y-component
                field_display(ilocnum+ieoff) = vectory
              else
                ! velocity z-component
                field_display(ilocnum+ieoff) = vectorz
              endif
            endif

          enddo
        endif ! USE_HIGHRES_FOR_MOVIES
      enddo ! NSPEC_SURFACE_EXT_MESH

      ! copy coordinate arrays since the sorting routine does not preserve them
      xp_save(:) = xp(:)
      yp_save(:) = yp(:)
      zp_save(:) = zp(:)
      

      ! sort the list based upon coordinates to get rid of multiples
      print *,'sorting list of points'
      call get_global(npointot,xp,yp,zp,iglob,locval,ifseg,nglob, &
           dble(minval(store_val_x(:))),dble(maxval(store_val_x(:))))

      ! print total number of points found
      print *
      print *,'found a total of ',nglob,' points'
      print *,'initial number of points (with multiples) was ',npointot


      !  normalize and scale vector field

      ! compute min and max of data value to normalize
      min_field_current = minval(field_display(:))
      max_field_current = maxval(field_display(:))

      if (plot_shaking_map) then
        ! print minimum and maximum amplitude in current snapshot
        print *
        print *,'minimum amplitude in current snapshot after removal = ',min_field_current
        print *,'maximum amplitude in current snapshot after removal = ',max_field_current
        print *
      else
        ! print minimum and maximum amplitude in current snapshot
        print *
        print *,'minimum amplitude in current snapshot = ',min_field_current
        print *,'maximum amplitude in current snapshot = ',max_field_current
        print *
      endif

      ! apply scaling in all cases for movies
      if (.not. plot_shaking_map) then

        ! normalizes values
        if (NORMALIZE_OUTPUT) then
          ! make sure range is always symmetric and center is in zero
          ! this assumption works only for fields that can be negative
          ! would not work for norm of vector for instance
          ! (we would lose half of the color palette if no negative values)
          max_absol = max(abs(min_field_current),abs(max_field_current))
          min_field_current = - max_absol
          max_field_current = + max_absol

          ! normalize field to [0:1]
          if (abs(max_field_current - min_field_current) > TINYVAL) &
            field_display(:) = (field_display(:) - min_field_current) / (max_field_current - min_field_current)

          ! rescale to [-1,1]
          field_display(:) = 2.*field_display(:) - 1.

          ! apply threshold to normalized field
          if (APPLY_THRESHOLD) &
            where(abs(field_display(:)) <= THRESHOLD) field_display = 0.
        endif

        ! apply non linear scaling to normalized field if needed
        if (NONLINEAR_SCALING) then
          where(field_display(:) >= 0.)
            field_display = field_display ** POWER_SCALING
          elsewhere
            field_display = - abs(field_display) ** POWER_SCALING
          endwhere
        endif

        ! normalizes values
        if (NORMALIZE_OUTPUT) then
          ! map back to [0,1]
          field_display(:) = (field_display(:) + 1.) / 2.

          ! map field to [0:255] for AVS color scale
          field_display(:) = 255. * field_display(:)
        endif

      ! apply scaling only if selected for shaking map
      else if (NONLINEAR_SCALING .and. iscaling_shake == 1) then

        ! normalize field to [0:1]
        if (abs(max_field_current) > TINYVAL) &
          field_display(:) = field_display(:) / max_field_current

        ! apply non linear scaling to normalized field
        field_display = field_display ** POWER_SCALING

        ! map field to [0:255] for AVS color scale
        field_display(:) = 255. * field_display(:)

      endif

      !--- ****** create AVS file using sorted list ******
      !--- writing bathymetry files
      write(outputname,"('/bathy',i4.4,'.txt')") ivalue
      open(unit=12,file=trim(OUTPUT_FILES)//outputname,status='unknown')
      all_data(1,:) = xp_save(:)
      all_data(2,:) = yp_save(:)
      all_data(3,:) = zp_save(:)
      all_data(4,:) = field_display(:)
      call cpu_time(start)
      print *,'starting to write bathymetry file: ', trim(OUTPUT_FILES)//outputname
      do i_v = 1, vmeshpts
        call neighbour_finder(coords(:,i_v), all_data, pt(:,1))
        write(12,*) pt(4,1)
      end do
      call cpu_time(finish)
      print *,'Time to finish writing bathymetry file ',finish-start

      close(12)

    ! end of loop and test on all the time steps for all the movie images
   endif
enddo ! it

  deallocate(all_data)
  deallocate(coords)

  deallocate(store_val_x)
  deallocate(store_val_y)
  deallocate(store_val_z)
  deallocate(store_val_ux)
  deallocate(store_val_uy)
  deallocate(store_val_uz)


  ! deallocate arrays for sorting routine
  deallocate(iglob,locval)
  deallocate(ifseg)
  deallocate(xp,yp,zp)
  deallocate(xp_save,yp_save,zp_save)
  deallocate(field_display)
  deallocate(mask_point)
  deallocate(ireorder)

  if (USE_HIGHRES_FOR_MOVIES) then
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(display)
  endif

  contains

  subroutine mesh2vec(A, meshpts)
    implicit none
    real(kind=CUSTOM_REAL), intent(out), allocatable :: A(:,:)
    integer, intent(out) :: meshpts

    integer, parameter:: MAX_STRING_LENGTH = 256
    integer::FID = 1
    integer::i,ier
     
    integer::npoints, nelements, index = 1, n_tri_meshp = 0, tmp, count = 0
    integer::pt0,pt1,pt2,pt3,pt4,pt5,pt6,pt7
    integer, allocatable :: c0(:), c1(:), c2(:)
    real, allocatable :: x(:), y(:), z(:)
    real, allocatable :: c_x(:,:)
    
    
    character(Len=MAX_STRING_LENGTH)::line
    
    open(unit=fid, file='../volna_ascii.msh', status='old', action='read', iostat=ier)
    if(ier /= 0) then
      print *, 'error in opening file: ', '../volna_ascii.msh'
    endif
    
    !skipping first for lines
    do i=1,4
      read(FID,'(a)') line
    end do
    
     read(FID, *) npoints
    
    allocate(x(npoints), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for x(:)'
    allocate(y(npoints), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for y(:)'
    allocate(z(npoints), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for z(:)'
    
    do i=1,npoints
    read(FID,*) tmp, x(i), y(i), z(i)
    end do
    
    !skipping next two lines
    do i=1,2
    read(FID,'(a)') line
    end do
    
    read(FID, *) nelements
    print*, 'npoints, nelements', npoints, nelements
  
    allocate(c0(nelements), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for c0(:)'
    allocate(c1(nelements), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for c1(:)'
    allocate(c2(nelements), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for c2(:)'
  
    do i=1,nelements
    read(FID,'(a)') line
    read(line, *, iostat=ier) pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7
    if(ier /= 0) then
      count = count + 1
    else
      c0(index) = pt5
      c1(index) = pt6
      c2(index) = pt7
  !    print *, c0(index), c1(index), c2(index)
      index = index + 1
    end if
    end do
  
    n_tri_meshp = index - 1;
  
    allocate(A(6,n_tri_meshp), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for c_x(:)'
    allocate(c_x(3,n_tri_meshp), stat=ier)
    if(ier /= 0) print *, 'error in allocating memory for c_x(:)'
 
    ! stop 
    do i=1,n_tri_meshp
      if (c0(i) < 1 .or. c1(i) < 1 .or. c2(i) < 1) then
        print *,'something went wrong'
        exit
      end if

      A(1,i) = (x(c0(i)) + x(c1(i)) + x(c2(i)))/3.0
      A(2,i) = (y(c0(i)) + y(c1(i)) + y(c2(i)))/3.0
      A(3,i) = (z(c0(i)) + z(c1(i)) + z(c2(i)))/3.0
    end do
    
  
    print *, count , ' Number of lines are not tranigular mesh points'
    meshpts = n_tri_meshp
    

    deallocate(x, y, z)
    deallocate(c0, c1, c2)
    deallocate(c_x)

    end subroutine mesh2vec

  subroutine neighbour_finder(point, mesh, neighbour)
    implicit none
    real(kind=CUSTOM_REAL), intent(in), dimension(:) :: point
    real(kind=CUSTOM_REAL), intent(in), dimension(:,:) :: mesh
    real(kind=CUSTOM_REAL), intent(out), dimension(:) :: neighbour

    real(kind=CUSTOM_REAL) :: distance, distance_min
    integer::i
    distance_min = (point(1) - mesh(1,1))**2 + (point(2) - mesh(2,1))**2 + (point(3) - mesh(3,1))**2
    neighbour = mesh(:,1)
    do i = 2, size(mesh(1,:))
      distance = (point(1) - mesh(1,i))**2 + (point(2) - mesh(2,i))**2 + (point(3) - mesh(3,i))**2
      if(distance < distance_min) then
        neighbour = mesh(:,i)
        distance_min = distance
      end if
    end do
  end subroutine neighbour_finder
end program create_movie_shakemap


