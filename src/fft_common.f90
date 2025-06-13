!! SPDX-License-Identifier: BSD-3-Clause

! This file contains common code shared by all FFT engines

integer, pointer, save :: format => null()  ! input X-pencil or Z-pencil

! Global size of the FFT
integer, pointer, save :: nx_fft => null(), ny_fft => null(), nz_fft => null()

! 2D processor grid
! FIXME this is already available in the module decomp_2d
integer, save, dimension(2) :: dims

! Decomposition objects
TYPE(DECOMP_INFO), pointer, save :: ph => null()  ! physical space
TYPE(DECOMP_INFO), pointer, save :: sp => null()  ! spectral space

! In-place FFT
logical, pointer, save :: inplace => null(), &
                          inplace_r2c => null(), &
                          inplace_c2r => null()

! Skip some c2c transforms
logical, pointer, save :: skip_x_c2c => null()
logical, pointer, save :: skip_y_c2c => null()
logical, pointer, save :: skip_z_c2c => null()

!
! Multiple grid options
!
! Number of FFT grids
integer, save :: n_grid = 0
type(decomp_2d_fft_engine), allocatable, target, save :: fft_engines(:)

public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
          decomp_2d_fft_finalize, decomp_2d_fft_get_size, &
          decomp_2d_fft_get_ph, decomp_2d_fft_get_sp, &
          decomp_2d_fft_get_ngrid, decomp_2d_fft_set_ngrid, &
          decomp_2d_fft_use_grid, decomp_2d_fft_engine, &
          decomp_2d_fft_get_engine, decomp_2d_fft_get_format, &
          decomp_2d_fft_get_inplace, decomp_2d_fft_get_inplace_r2c, &
          decomp_2d_fft_get_inplace_c2r

! Declare generic interfaces to handle different inputs

interface decomp_2d_fft_init
   module procedure fft_init_noarg
   module procedure fft_init_arg
   module procedure fft_init_one_grid
   module procedure fft_init_multiple_grids
end interface

interface decomp_2d_fft_3d
   module procedure fft_3d_c2c
   module procedure fft_3d_r2c
   module procedure fft_3d_c2r
end interface

interface
   module subroutine decomp_2d_fft_log(backend)
      character(len=*), intent(in) :: backend
   end subroutine decomp_2d_fft_log
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise the FFT module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_init_noarg

   implicit none

   call fft_init_arg(PHYSICAL_IN_X)  ! default input is X-pencil data

   return
end subroutine fft_init_noarg

subroutine fft_init_arg(pencil, opt_skip_XYZ_c2c)     ! allow to handle Z-pencil input

   implicit none

   integer, intent(IN) :: pencil
   logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c

   call fft_init_one_grid(pencil, nx_global, ny_global, nz_global, &
                          opt_skip_XYZ_c2c=opt_skip_XYZ_c2c)

   return
end subroutine fft_init_arg

! Initialise the FFT library to perform arbitrary size transforms
subroutine fft_init_one_grid(pencil, nx, ny, nz, &
                             opt_inplace, &
                             opt_inplace_r2c, &
                             opt_inplace_c2r, &
                             opt_skip_XYZ_c2c)

   implicit none

   integer, intent(IN) :: pencil, nx, ny, nz
   logical, intent(in), optional :: opt_inplace, opt_inplace_r2c, opt_inplace_c2r
   logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c

   ! Only one FFT engine will be used
   call decomp_2d_fft_set_ngrid(1)

   ! Initialise the FFT engine
   call decomp_2d_fft_init(pencil, nx, ny, nz, 1, &
                           opt_inplace, &
                           opt_inplace_r2c, &
                           opt_inplace_c2r, &
                           opt_skip_XYZ_c2c=opt_skip_XYZ_c2c)

end subroutine fft_init_one_grid

! Initialise the provided FFT grid
subroutine fft_init_multiple_grids(pencil, nx, ny, nz, igrid, &
                                   opt_inplace, &
                                   opt_inplace_r2c, &
                                   opt_inplace_c2r, &
                                   opt_skip_XYZ_c2c)

   implicit none

   integer, intent(in) :: pencil, nx, ny, nz, igrid
   logical, intent(in), optional :: opt_inplace, opt_inplace_r2c, opt_inplace_c2r
   logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c

   ! Safety check
   if (igrid < 1 .or. igrid > n_grid) then
      call decomp_2d_abort(__FILE__, __LINE__, igrid, "Invalid value for igrid")
   end if

   ! Initialise the engine
   call fft_engines(igrid)%init(pencil, nx, ny, nz, &
                                opt_inplace, &
                                opt_inplace_r2c, &
                                opt_inplace_c2r, &
                                opt_skip_XYZ_c2c)

end subroutine fft_init_multiple_grids

! Initialise the given FFT engine
subroutine decomp_2d_fft_engine_init(engine, pencil, nx, ny, nz, &
                                     opt_inplace, &
                                     opt_inplace_r2c, &
                                     opt_inplace_c2r, &
                                     opt_skip_XYZ_c2c)

   implicit none

   class(decomp_2d_fft_engine), target, intent(inout) :: engine
   integer, intent(in) :: pencil, nx, ny, nz
   logical, intent(in), optional :: opt_inplace, opt_inplace_r2c, opt_inplace_c2r
   logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c

   if (decomp_profiler_fft) call decomp_profiler_start("fft_init")

   ! Safety checks
   if (engine%initialised) then
      call decomp_2d_abort(__FILE__, __LINE__, 4, &
                           'FFT engine should only be initialised once')
   end if
   if (nx <= 0) call decomp_2d_abort(__FILE__, __LINE__, nx, "Invalid value for nx")
   if (ny <= 0) call decomp_2d_abort(__FILE__, __LINE__, ny, "Invalid value for ny")
   if (nz <= 0) call decomp_2d_abort(__FILE__, __LINE__, nz, "Invalid value for nz")

   ! Store the key parameters
   engine%format = pencil
   engine%nx_fft = nx
   engine%ny_fft = ny
   engine%nz_fft = nz

   ! FFT can be inplace
   if (present(opt_inplace)) then
      engine%inplace = opt_inplace
   else
      engine%inplace = DECOMP_2D_FFT_INPLACE
   end if
   if (present(opt_inplace_r2c)) then
      call decomp_2d_abort(__FILE__, __LINE__, -1, 'In-place r2c transform is not yet available')
   end if
   if (present(opt_inplace_c2r)) then
      call decomp_2d_abort(__FILE__, __LINE__, -1, 'In-place c2r transform is not yet available')
   end if

   ! some c2c transforms can be skipped
   if (present(opt_skip_XYZ_c2c)) then
      engine%skip_x_c2c = opt_skip_XYZ_c2c(1)
      engine%skip_y_c2c = opt_skip_XYZ_c2c(2)
      engine%skip_z_c2c = opt_skip_XYZ_c2c(3)
   else
      engine%skip_x_c2c = .false.
      engine%skip_y_c2c = .false.
      engine%skip_z_c2c = .false.
   end if

   ! determine the processor grid in use
   dims = get_decomp_dims()

   ! for c2r/r2c interface:
   ! if in physical space, a real array is of size: nx*ny*nz
   ! in spectral space, the complex array is of size:
   !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
   !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

   if (nx == nx_global .and. ny == ny_global .and. nz == nz_global) then
      engine%ph => decomp_main
   else
      call decomp_info_init(nx, ny, nz, engine%ph_target)
      engine%ph => engine%ph_target
   end if
   if (pencil == PHYSICAL_IN_X) then
      call decomp_info_init(nx / 2 + 1, ny, nz, engine%sp)
   else if (pencil == PHYSICAL_IN_Z) then
      call decomp_info_init(nx, ny, nz / 2 + 1, engine%sp)
   else
      call decomp_2d_abort(__FILE__, __LINE__, pencil, "Invalid value for pencil")
   end if
   !
   ! In case of c2c : assume 2decomp is initialized with complex_pool = .true.
   !
   ! In case of r2c / c2r : assume 2decomp is initialized with complex_pool = .false.
   !                        the line below will make sure complex arrays fit in the memory pool
   !
   if (use_pool) call decomp_pool%new_shape(complex_type, engine%sp)
#ifdef EVEN
   if (use_pool) call decomp_pool%new_shape(complex_type, shp=(/max(engine%sp%x1count * dims(1), engine%sp%y2count * dims(2))/))
#endif

   ! Warning : replace the default engine
   call engine%use_it(opt_force=.true.)

   ! Compute engine-specific stuff
   call init_fft_engine(engine)

   ! Tag the engine as initialised
   engine%initialised = .true.

   ! All the components of the default engine must be updated
   call engine%use_it()

   if (decomp_profiler_fft) call decomp_profiler_end("fft_init")

end subroutine decomp_2d_fft_engine_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final clean up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_finalize

   implicit none

   integer :: igrid

   ! Clean each engine
   do igrid = 1, n_grid
      call fft_engines(igrid)%fin()
   end do

   ! Free memory
   n_grid = 0
   deallocate (fft_engines)

   ! Nullify all the pointers in the module
   nullify (format)
   nullify (nx_fft)
   nullify (ny_fft)
   nullify (nz_fft)
   nullify (ph)
   nullify (sp)
   nullify (inplace)
   nullify (inplace_r2c)
   nullify (inplace_c2r)
   nullify (skip_x_c2c)
   nullify (skip_y_c2c)
   nullify (skip_z_c2c)

   ! Clean engine-specific stuff located in the module
   call finalize_fft_engine()

end subroutine decomp_2d_fft_finalize

! Clean the given FFT engine
subroutine decomp_2d_fft_engine_fin(engine)

   implicit none

   class(decomp_2d_fft_engine), intent(inout) :: engine

   if (decomp_profiler_fft) call decomp_profiler_start("fft_fin")

   if (engine%nx_fft /= nx_global .or. &
       engine%ny_fft /= ny_global .or. &
       engine%nz_fft /= nz_global) then
      call decomp_info_finalize(engine%ph_target)
   end if
   nullify (engine%ph)
   call decomp_info_finalize(engine%sp)

   call finalize_fft_engine(engine)

   engine%initialised = .false.

   if (decomp_profiler_fft) call decomp_profiler_end("fft_fin")

end subroutine decomp_2d_fft_engine_fin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the size, starting/ending index of the distributed array
!  whose global size is (nx/2+1)*ny*nz, for defining data structures
!  in r2c and c2r interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_get_size(istart, iend, isize)

   implicit none
   integer, dimension(3), intent(OUT) :: istart, iend, isize

   if (format == PHYSICAL_IN_X) then
      istart = sp%zst
      iend = sp%zen
      isize = sp%zsz
   else if (format == PHYSICAL_IN_Z) then
      istart = sp%xst
      iend = sp%xen
      isize = sp%xsz
   end if

   return
end subroutine decomp_2d_fft_get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return a pointer to the decomp_info object ph
!
! The caller should not apply decomp_info_finalize on the pointer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function decomp_2d_fft_get_ph()

   implicit none

   type(decomp_info), pointer :: decomp_2d_fft_get_ph

   if (.not. associated(ph)) then
      call decomp_2d_abort(__FILE__, __LINE__, -1, 'FFT library must be initialised first')
   end if
   decomp_2d_fft_get_ph => ph

end function decomp_2d_fft_get_ph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return a pointer to the decomp_info object sp
!
! The caller should not apply decomp_info_finalize on the pointer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function decomp_2d_fft_get_sp()

   implicit none

   type(decomp_info), pointer :: decomp_2d_fft_get_sp

   if (.not. associated(ph)) then
      call decomp_2d_abort(__FILE__, __LINE__, -1, 'FFT library must be initialised first')
   end if
   decomp_2d_fft_get_sp => sp

end function decomp_2d_fft_get_sp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the number of FFT grids
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function decomp_2d_fft_get_ngrid()

   implicit none

   integer :: decomp_2d_fft_get_ngrid

   decomp_2d_fft_get_ngrid = n_grid

end function decomp_2d_fft_get_ngrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the number of FFT grids
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_set_ngrid(ngrd)

   implicit none

   integer, intent(in) :: ngrd

   type(decomp_2d_fft_engine), allocatable, target :: tmp(:)

   ! Safety check
   if (ngrd < 1) then
      call decomp_2d_abort(__FILE__, __LINE__, ngrd, "Invalid value for n_grid")
   end if

   ! Change the number of grids
   if (n_grid > 0 .or. allocated(fft_engines)) then
      ! Save the current grids and deallocate
      call move_alloc(fft_engines, tmp)
      ! Re-allocate
      allocate (fft_engines(ngrd))
      ! Restore
      fft_engines(1:min(ngrd, n_grid)) = tmp(1:min(ngrd, n_grid))
      ! Set the number of FFT grids
      n_grid = ngrd
      ! Free memory and return
      deallocate (tmp)
      return
   end if

   ! Set the number of FFT grids
   n_grid = ngrd
   allocate (fft_engines(n_grid))

end subroutine decomp_2d_fft_set_ngrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make the provided grid the default one
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_use_grid(igrid)

   implicit none

   integer, optional, intent(in) :: igrid

   if (present(igrid)) then

      ! Safety check
      if (igrid < 1 .or. igrid > n_grid) then
         call decomp_2d_abort(__FILE__, __LINE__, igrid, "Invalid value for igrid")
      end if
      call fft_engines(igrid)%use_it()

   else

      ! Safety check
      if (n_grid < 1 .or. .not. allocated(fft_engines)) then
         call decomp_2d_abort(__FILE__, __LINE__, n_grid, &
                              "The FFT module was not initialised")
      end if
      call fft_engines(1)%use_it()

   end if

end subroutine decomp_2d_fft_use_grid

! Associate the pointers in the module to the ones in the engine
subroutine decomp_2d_fft_engine_use_it(engine, opt_force)

   implicit none

   class(decomp_2d_fft_engine), intent(in), target :: engine
   logical, intent(in), optional :: opt_force ! Skip safety check

   ! Safety check
   if (.not. engine%initialised) then
      if (present(opt_force)) then
         if (.not. opt_force) then
            call decomp_2d_abort(__FILE__, __LINE__, -1, "FFT engine is not ready")
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, "FFT engine is not ready")
      end if
   end if

   ! All the engines have this
   format => engine%format
   nx_fft => engine%nx_fft
   ny_fft => engine%ny_fft
   nz_fft => engine%nz_fft
   ph => engine%ph
   sp => engine%sp
   inplace => engine%inplace
   skip_x_c2c => engine%skip_x_c2c
   skip_y_c2c => engine%skip_y_c2c
   skip_z_c2c => engine%skip_z_c2c

   ! Engine-specific stuff
   call use_fft_engine(engine)

end subroutine decomp_2d_fft_engine_use_it

! The external code can obtain direct access to internal FFT engines
!
! Warning : The external code should not clean the provided engine
!           The subroutine decomp_2d_fft_finalize will do it
function decomp_2d_fft_get_engine(igrid)

   implicit none

   type(decomp_2d_fft_engine), pointer :: decomp_2d_fft_get_engine
   integer, intent(in), optional :: igrid

   if (present(igrid)) then

      ! Safety check
      if (igrid < 1 .or. igrid > n_grid) then
         call decomp_2d_abort(__FILE__, __LINE__, igrid, "Invalid value for igrid")
      end if
      decomp_2d_fft_get_engine => fft_engines(igrid)

   else

      ! Safety check
      if (n_grid < 1 .or. .not. allocated(fft_engines)) then
         call decomp_2d_abort(__FILE__, __LINE__, n_grid, &
                              "The FFT module was not initialised")
      end if
      decomp_2d_fft_get_engine => fft_engines(1)

   end if

end function decomp_2d_fft_get_engine

! The external code can check if the engine is physical_in_x or physical_in_z
function decomp_2d_fft_get_format()

   implicit none

   integer :: decomp_2d_fft_get_format

   decomp_2d_fft_get_format = format

end function decomp_2d_fft_get_format

! The external code can check if the FFT is inplace
function decomp_2d_fft_get_inplace()

   implicit none

   logical :: decomp_2d_fft_get_inplace

   decomp_2d_fft_get_inplace = inplace

end function decomp_2d_fft_get_inplace

! The external code can check if the r2c is inplace
function decomp_2d_fft_get_inplace_r2c()

   implicit none

   logical :: decomp_2d_fft_get_inplace_r2c

   decomp_2d_fft_get_inplace_r2c = .false. ! inplace_r2c

end function decomp_2d_fft_get_inplace_r2c

! The external code can check if the c2r is inplace
function decomp_2d_fft_get_inplace_c2r()

   implicit none

   logical :: decomp_2d_fft_get_inplace_c2r

   decomp_2d_fft_get_inplace_c2r = .false. ! inplace_c2r

end function decomp_2d_fft_get_inplace_c2r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Wrappers for calling 3D FFT directly using the engine object
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_engine_fft_c2c(engine, in, out, isign)

   implicit none

   class(decomp_2d_fft_engine), intent(in), target :: engine
   complex(mytype), dimension(:, :, :), intent(INOUT) :: in
   complex(mytype), dimension(:, :, :), intent(OUT) :: out
   integer, intent(IN) :: isign

   call engine%use_it()
   call decomp_2d_fft_3d(in, out, isign)

end subroutine decomp_2d_fft_engine_fft_c2c

subroutine decomp_2d_fft_engine_fft_r2c(engine, in, out)

   implicit none

   class(decomp_2d_fft_engine), intent(in), target :: engine
   real(mytype), dimension(:, :, :), intent(INOUT) :: in
   complex(mytype), dimension(:, :, :), intent(OUT) :: out

   call engine%use_it()
   call decomp_2d_fft_3d(in, out)

end subroutine decomp_2d_fft_engine_fft_r2c

subroutine decomp_2d_fft_engine_fft_c2r(engine, in, out)

   implicit none

   class(decomp_2d_fft_engine), intent(in), target :: engine
   complex(mytype), dimension(:, :, :), intent(INOUT) :: in
   real(mytype), dimension(:, :, :), intent(OUT) :: out

   call engine%use_it()
   call decomp_2d_fft_3d(in, out)

end subroutine decomp_2d_fft_engine_fft_c2r
