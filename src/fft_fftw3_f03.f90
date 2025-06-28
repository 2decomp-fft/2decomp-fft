!! SPDX-License-Identifier: BSD-3-Clause

! This is the FFTW implementation of the FFT library using
! the Fortran 2003 interface introduced in FFTW 3.3-beta1

module decomp_2d_fft

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_profiler
   use, intrinsic :: iso_c_binding
   use m_decomp_pool

   implicit none

   include "fftw3.f03"

   private        ! Make everything private unless declared public

   ! engine-specific global variables
   integer, parameter :: plan_type = FFTW_MEASURE
   integer, parameter :: plan_type_dtt = FFTW_MEASURE + FFTW_UNALIGNED

   ! FFTW plans
   ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
   ! For c2c transforms:
   !     use plan(-1,j) for  forward transform;
   !     use plan( 1,j) for backward transform;
   ! For r2c/c2r transforms:
   !     use plan(0,j) for r2c transforms;
   !     use plan(2,j) for c2r transforms;
   type(C_PTR), contiguous, pointer, save :: plan(:, :) => null()

   integer, parameter, public :: D2D_FFT_BACKEND = D2D_FFT_BACKEND_FFTW3_F03

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

   ! Flag to check if DTT is active
   logical, pointer, save :: with_dtt => null()

   ! array with the DTT setup
   integer, contiguous, pointer, save :: dtt(:) => null()

   ! FFTW plans for DTT
   !    1, 2, 3 : forward transform in X, Y and Z
   !    4, 5, 6 : backward transform in X, Y and Z
   type(C_PTR), contiguous, pointer, save :: dtt_plan(:) => null()

   ! Derived type with all the quantities needed to perform FFT
   type decomp_2d_fft_engine
      type(c_ptr), private :: plan(-1:2, 3)
      integer, private :: format
      logical, private :: initialised = .false.
      integer, private :: nx_fft, ny_fft, nz_fft
      type(decomp_info), pointer, public :: ph => null()
      type(decomp_info), private :: ph_target ! ph => ph_target or ph => decomp_main
      type(decomp_info), public :: sp
      logical, private :: inplace, inplace_r2c, inplace_c2r
      logical, private :: skip_x_c2c, skip_y_c2c, skip_z_c2c
      ! Below is specific to DTT
      logical, public :: with_dtt
      integer, allocatable, public, dimension(:) :: dtt
      type(C_PTR), private :: dtt_plan(6)
   contains
      procedure, public :: init => decomp_2d_fft_engine_init
      procedure, public :: fin => decomp_2d_fft_engine_fin
      procedure, public :: use_it => decomp_2d_fft_engine_use_it
      procedure, private :: dtt_init => decomp_2d_fft_engine_dtt_init
      procedure, private :: dtt_fin => decomp_2d_fft_engine_dtt_fin
      procedure, public :: dtt_input => decomp_2d_fft_engine_dtt_input
      generic, public :: fft => c2c, r2c, c2r
      procedure, private :: c2c => decomp_2d_fft_engine_fft_c2c
      procedure, private :: r2c => decomp_2d_fft_engine_fft_r2c
      procedure, private :: c2r => decomp_2d_fft_engine_fft_c2r
   end type decomp_2d_fft_engine

   !
   ! Multiple grid options
   !
   ! Number of FFT grids
   integer, save :: n_grid = 0
   ! Array to store all the grids
   type(decomp_2d_fft_engine), allocatable, target, save :: fft_engines(:)

   public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
             decomp_2d_dtt_3d_r2r, &
             decomp_2d_fft_finalize, decomp_2d_fft_get_size, &
             decomp_2d_fft_get_ph, decomp_2d_fft_get_sp, &
             decomp_2d_fft_get_dtt_input, &
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

   end subroutine fft_init_noarg

   subroutine fft_init_arg(pencil, opt_skip_XYZ_c2c, opt_DTT)     ! allow to handle Z-pencil input

      implicit none

      integer, intent(IN) :: pencil
      logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c
      integer, dimension(:), intent(in), optional :: opt_DTT

      call fft_init_one_grid(pencil, nx_global, ny_global, nz_global, &
                             opt_skip_XYZ_c2c=opt_skip_XYZ_c2c, &
                             opt_DTT=opt_DTT)

   end subroutine fft_init_arg

   ! Initialise one FFT library to perform arbitrary size transforms
   subroutine fft_init_one_grid(pencil, nx, ny, nz, &
                                opt_inplace, &
                                opt_inplace_r2c, &
                                opt_inplace_c2r, &
                                opt_skip_XYZ_c2c, &
                                opt_DTT)

      implicit none

      integer, intent(IN) :: pencil, nx, ny, nz
      logical, intent(in), optional :: opt_inplace, opt_inplace_r2c, opt_inplace_c2r
      logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c
      integer, dimension(:), intent(in), optional :: opt_DTT

      ! Only one FFT engine will be used
      call decomp_2d_fft_set_ngrid(1)

      ! Initialise the FFT engine
      call decomp_2d_fft_init(pencil, nx, ny, nz, 1, &
                              opt_inplace, &
                              opt_inplace_r2c, &
                              opt_inplace_c2r, &
                              opt_skip_XYZ_c2c=opt_skip_XYZ_c2c, &
                              opt_DTT=opt_DTT)

   end subroutine fft_init_one_grid

   ! Initialise the provided FFT grid
   subroutine fft_init_multiple_grids(pencil, nx, ny, nz, igrid, &
                                      opt_inplace, &
                                      opt_inplace_r2c, &
                                      opt_inplace_c2r, &
                                      opt_skip_XYZ_c2c, &
                                      opt_DTT)

      implicit none

      integer, intent(in) :: pencil, nx, ny, nz, igrid
      logical, intent(in), optional :: opt_inplace, opt_inplace_r2c, opt_inplace_c2r
      logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c
      integer, dimension(:), intent(in), optional :: opt_DTT

      ! Safety check
      if (igrid < 1 .or. igrid > n_grid) then
         call decomp_2d_abort(__FILE__, __LINE__, igrid, "Invalid value for igrid")
      end if

      ! Initialise the engine
      call fft_engines(igrid)%init(pencil, nx, ny, nz, &
                                   opt_inplace, &
                                   opt_inplace_r2c, &
                                   opt_inplace_c2r, &
                                   opt_skip_XYZ_c2c, &
                                   opt_DTT=opt_DTT)

   end subroutine fft_init_multiple_grids

   ! Initialise the given FFT engine
   subroutine decomp_2d_fft_engine_init(engine, pencil, nx, ny, nz, &
                                        opt_inplace, &
                                        opt_inplace_r2c, &
                                        opt_inplace_c2r, &
                                        opt_skip_XYZ_c2c, &
                                        opt_DTT)

      implicit none

      class(decomp_2d_fft_engine), intent(inout), target :: engine
      integer, intent(in) :: pencil, nx, ny, nz
      logical, intent(in), optional :: opt_inplace, opt_inplace_r2c, opt_inplace_c2r
      logical, dimension(3), intent(in), optional :: opt_skip_XYZ_c2c
      integer, dimension(:), intent(in), optional :: opt_DTT

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
      if (engine%inplace .and. present(opt_inplace_r2c)) then
         engine%inplace_r2c = opt_inplace_r2c
      else
         engine%inplace_r2c = .false.!DECOMP_2D_FFT_INPLACE ! this is experimental
      end if
      if (engine%inplace .and. present(opt_inplace_c2r)) then
         engine%inplace_c2r = opt_inplace_c2r
      else
         engine%inplace_c2r = .false.!DECOMP_2D_FFT_INPLACE ! this is experimental
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

      if (nx == nx_global .and. &
          ny == ny_global .and. &
          nz == nz_global) then
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

      ! Prepare the DTT components
      if (present(opt_DTT)) then
         engine%with_dtt = .true.
         call engine%dtt_init(opt_DTT)
      else
         engine%with_dtt = .false.
      end if

      ! Warning : replace the default engine
      call engine%use_it(opt_force=.true.)

      ! Compute the plan for the current engine
      call init_fft_engine()

      ! Tag the engine as initialised
      engine%initialised = .true.

      ! All the components of the default engine must be updated
      call engine%use_it()

      if (decomp_profiler_fft) call decomp_profiler_end("fft_init")

   end subroutine decomp_2d_fft_engine_init

   ! Initialise the DTT components of the provided engine
   subroutine decomp_2d_fft_engine_dtt_init(engine, in_DTT)

      implicit none

      class(decomp_2d_fft_engine), intent(inout), target :: engine
      integer, dimension(:), intent(in) :: in_DTT

      ! Safety check
      if (size(in_DTT) < 3) call decomp_2d_abort(__FILE__, __LINE__, size(in_DTT), "Invalid argument")
      if (minval(in_DTT(1:3)) < 1) call decomp_2d_abort(__FILE__, __LINE__, minval(in_DTT(1:3)), "Invalid argument")
      if (maxval(in_DTT(1:3)) > 8) call decomp_2d_abort(__FILE__, __LINE__, maxval(in_DTT(1:3)), "Invalid argument")

      ! Prepare engine%dtt
      ! Mandatory
      !    1:3 => type of forward transform
      ! Optional, default values in dtt_assign_default
      !    4:6 => ifirst, index where the in-place r2r transform starts
      !    7:9 => ndismiss, number of points skipped
      !    10:12 => ofirst, should match ifirst (in-place r2r transform)
      ! Values defined in dtt_invert
      !    12:15 => type of backward transform
      allocate (engine%dtt(15))
      if (size(in_DTT) == 12) then
         engine%dtt(1:12) = in_DTT
         if (any(engine%dtt(4:6) /= engine%dtt(10:12))) then
            call decomp_2d_abort(__FILE__, __LINE__, 1, "Setup is not compatible with in-place r2r transforms")
         end if
      elseif (size(in_DTT) == 3) then
         engine%dtt(1:3) = in_DTT
         call dtt_assign_default(engine%dtt)
      else
         call decomp_2d_abort(__FILE__, __LINE__, size(in_DTT), "Invalid argument")
      end if
      call dtt_invert(engine%dtt)
      call dtt_for_fftw(engine%dtt(1:3))
      call dtt_for_fftw(engine%dtt(13:15))

      ! Prepare the fftw plans
      engine%dtt_plan = c_null_ptr
      ! in x
      call r2r_1m_x_plan(engine%dtt_plan(1), engine%ph, engine%dtt(1), engine%dtt(7))
      if (engine%dtt(13) /= engine%dtt(1)) call r2r_1m_x_plan(engine%dtt_plan(4), engine%ph, engine%dtt(13), engine%dtt(7))
      ! in y
      call r2r_1m_y_plan(engine%dtt_plan(2), engine%ph, engine%dtt(2), engine%dtt(8))
      if (engine%dtt(14) /= engine%dtt(2)) call r2r_1m_y_plan(engine%dtt_plan(5), engine%ph, engine%dtt(14), engine%dtt(8))
      ! in z
      call r2r_1m_z_plan(engine%dtt_plan(3), engine%ph, engine%dtt(3), engine%dtt(9))
      if (engine%dtt(15) /= engine%dtt(3)) call r2r_1m_z_plan(engine%dtt_plan(6), engine%ph, engine%dtt(15), engine%dtt(9))

   end subroutine decomp_2d_fft_engine_dtt_init

   ! Set default values in the DTT config
   ! FIXME add constants in the module decomp_2d_constants ?
   subroutine dtt_assign_default(arg_dtt)

      implicit none

      integer, intent(inout) :: arg_dtt(15)

      integer :: k

      ! Generic values
      ! ifirst = 1, use the array from the beginning
      arg_dtt(4:6) = 1
      ! ndismiss = 1, skip one point
      arg_dtt(7:9) = 1

      ! Specific cases
      do k = 1, 3
         if (arg_dtt(k) == 0) arg_dtt(k + 6) = 0 ! Periodic, ndismiss = 0
         if (arg_dtt(k) == 1) arg_dtt(k + 6) = 0 ! DCT1, ndismiss = 0
         if (arg_dtt(k) == 5) arg_dtt(k + 3) = 2 ! DST1, ifirst = 2
         if (arg_dtt(k) == 5) arg_dtt(k + 6) = 2 ! DST1, ndismiss = 2
         if (arg_dtt(k) == 7) arg_dtt(k + 3) = 2 ! DST3, ifirst = 2
      end do

      ! ofirst = ifirst
      arg_dtt(10:12) = arg_dtt(4:6)

   end subroutine dtt_assign_default

   ! Set the backward transforms in the DTT config
   ! FIXME add constants in the module decomp_2d_constants ?
   subroutine dtt_invert(arg_dtt)

      implicit none

      integer, intent(inout) :: arg_dtt(15)

      integer :: k

      ! Default : inv(DTT) = DTT
      arg_dtt(13:15) = arg_dtt(1:3)

      ! Except for DCT2, DCT3, DST2, DST3
      do k = 1, 3
         if (arg_dtt(k) == 0) arg_dtt(12 + k) = 9
         if (arg_dtt(k) == 2) arg_dtt(12 + k) = 3 ! inv(DCT2) = DCT3
         if (arg_dtt(k) == 3) arg_dtt(12 + k) = 2 ! inv(DCT3) = DCT2
         if (arg_dtt(k) == 6) arg_dtt(12 + k) = 7 ! inv(DST2) = DST3
         if (arg_dtt(k) == 7) arg_dtt(12 + k) = 6 ! inv(DST3) = DST2
      end do

   end subroutine dtt_invert

   ! Adapt the DTT type to FFTW
   ! FIXME add constants in the module decomp_2d_constants ?
   subroutine dtt_for_fftw(arg_dtt)

      implicit none

      integer, intent(inout) :: arg_dtt(:)

      integer :: k

      do k = 1, size(arg_dtt)
         select case (arg_dtt(k))
         case (0)
            arg_dtt(k) = FFTW_FORWARD
         case (1)
            arg_dtt(k) = FFTW_REDFT00
         case (2)
            arg_dtt(k) = FFTW_REDFT10
         case (3)
            arg_dtt(k) = FFTW_REDFT01
         case (4)
            arg_dtt(k) = FFTW_REDFT11
         case (5)
            arg_dtt(k) = FFTW_RODFT00
         case (6)
            arg_dtt(k) = FFTW_RODFT10
         case (7)
            arg_dtt(k) = FFTW_RODFT01
         case (8)
            arg_dtt(k) = FFTW_RODFT11
         case (9)
            arg_dtt(k) = FFTW_BACKWARD
         case default
            call decomp_2d_abort(__FILE__, __LINE__, arg_dtt(k), "Invalid value")
         end select
      end do

   end subroutine dtt_for_fftw

   ! Adapt the FFTW type to DTT
   ! FIXME add constants in the module decomp_2d_constants ?
   function fftw_for_dtt(arg_dtt)

      implicit none

      integer, dimension(3) :: fftw_for_dtt
      integer, intent(in) :: arg_dtt(3)

      ! Local variables
      integer :: i

      ! Convert FFTW transforms into 2decomp transforms
      do i = 1, 3
         select case (arg_dtt(i))
         case (FFTW_FORWARD)
            fftw_for_dtt(i) = 0
         case (FFTW_REDFT00)
            fftw_for_dtt(i) = 1
         case (FFTW_REDFT10)
            fftw_for_dtt(i) = 2
         case (FFTW_REDFT01)
            fftw_for_dtt(i) = 3
         case (FFTW_REDFT11)
            fftw_for_dtt(i) = 4
         case (FFTW_RODFT00)
            fftw_for_dtt(i) = 5
         case (FFTW_RODFT10)
            fftw_for_dtt(i) = 6
         case (FFTW_RODFT01)
            fftw_for_dtt(i) = 7
         case (FFTW_RODFT11)
            fftw_for_dtt(i) = 8
         case default
            call decomp_2d_abort(__FILE__, __LINE__, arg_dtt(i), "Invalid value")
         end select
      end do

   end function fftw_for_dtt

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
      nullify (plan)
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
      nullify (with_dtt)
      if (associated(dtt)) then
         nullify (dtt)
         nullify (dtt_plan)
      end if

      ! Clean the FFTW library
#ifdef FFTW_omp
      call fftw_cleanup_threads()
#else
      call fftw_cleanup()
#endif


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

      call finalize_fft_engine(engine%plan)

      if (engine%with_dtt) call engine%dtt_fin()

      engine%initialised = .false.

      if (decomp_profiler_fft) call decomp_profiler_end("fft_fin")

   end subroutine decomp_2d_fft_engine_fin

   ! Clean the DTT components of the provided engine
   subroutine decomp_2d_fft_engine_dtt_fin(engine)

      implicit none

      class(decomp_2d_fft_engine), intent(inout) :: engine

      integer :: i

      ! Restore the dtt flag
      engine%with_dtt = .false.

      ! Deallocate the DTT config
      deallocate (engine%dtt)

      ! Clean the fftw plans
      do i = 1, 6
         if (c_associated(engine%dtt_plan(i))) then
#ifdef DOUBLE_PREC
            call fftw_destroy_plan(engine%dtt_plan(i))
#else
            call fftwf_destroy_plan(engine%dtt_plan(i))
#endif
         end if
      end do
      engine%dtt_plan = c_null_ptr

   end subroutine decomp_2d_fft_engine_dtt_fin

   ! Return the type of the transform in each direction
   function decomp_2d_fft_get_dtt_input()

      implicit none

      integer, dimension(3) :: decomp_2d_fft_get_dtt_input

      if (.not. associated(with_dtt)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid operation")

      if (with_dtt) then
         decomp_2d_fft_get_dtt_input = fftw_for_dtt(dtt(1:3))
      else
         call decomp_2d_abort(__FILE__, __LINE__, 1, "No DTT for the current engine")
      end if

   end function decomp_2d_fft_get_dtt_input

   ! Return the type of the transform in each direction
   function decomp_2d_fft_engine_dtt_input(engine)

      implicit none

      class(decomp_2d_fft_engine), intent(in) :: engine
      integer, dimension(3) :: decomp_2d_fft_engine_dtt_input

      if (engine%with_dtt) then
         decomp_2d_fft_engine_dtt_input = fftw_for_dtt(engine%dtt(1:3))
      else
         call decomp_2d_abort(__FILE__, __LINE__, 1, "No DTT for the current engine")
      end if

   end function decomp_2d_fft_engine_dtt_input

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

      type(decomp_2d_fft_engine), allocatable :: tmp(:)

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

      plan => engine%plan
      format => engine%format
      nx_fft => engine%nx_fft
      ny_fft => engine%ny_fft
      nz_fft => engine%nz_fft
      ph => engine%ph
      sp => engine%sp
      inplace => engine%inplace
      inplace_r2c => engine%inplace_r2c
      inplace_c2r => engine%inplace_c2r
      skip_x_c2c => engine%skip_x_c2c
      skip_y_c2c => engine%skip_y_c2c
      skip_z_c2c => engine%skip_z_c2c
      with_dtt => engine%with_dtt

      if (with_dtt) then
         dtt => engine%dtt
         dtt_plan => engine%dtt_plan
      else if (associated(dtt)) then
         nullify (dtt)
         nullify (dtt_plan)
      end if

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

      decomp_2d_fft_get_inplace_r2c = inplace_r2c

   end function decomp_2d_fft_get_inplace_r2c

   ! The external code can check if the c2r is inplace
   function decomp_2d_fft_get_inplace_c2r()

      implicit none

      logical :: decomp_2d_fft_get_inplace_c2r

      decomp_2d_fft_get_inplace_c2r = inplace_c2r

   end function decomp_2d_fft_get_inplace_c2r

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in X direction
   subroutine c2c_1m_x_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

#ifdef DOUBLE_PREC
      complex(C_DOUBLE_COMPLEX), contiguous, pointer :: a1(:, :, :)
      complex(C_DOUBLE_COMPLEX), contiguous, pointer :: a1o(:, :, :)
#else
      complex(C_FLOAT_COMPLEX), contiguous, pointer :: a1(:, :, :)
      complex(C_FLOAT_COMPLEX), contiguous, pointer :: a1o(:, :, :)
#endif
      type(C_PTR) :: a1_p
      integer(C_SIZE_T) :: sz

      sz = decomp%xsz(1) * decomp%xsz(2) * decomp%xsz(3)
      a1_p = fftw_alloc_complex(sz)
      call c_f_pointer(a1_p, a1, [decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)])
      call c_f_pointer(a1_p, a1o, [decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)])

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft(1, decomp%xsz(1), &
                                 decomp%xsz(2) * decomp%xsz(3), a1, decomp%xsz(1), 1, &
                                 decomp%xsz(1), a1o, decomp%xsz(1), 1, decomp%xsz(1), &
                                 isign, plan_type)
#else
      plan1 = fftwf_plan_many_dft(1, decomp%xsz(1), &
                                  decomp%xsz(2) * decomp%xsz(3), a1, decomp%xsz(1), 1, &
                                  decomp%xsz(1), a1o, decomp%xsz(1), 1, decomp%xsz(1), &
                                  isign, plan_type)
#endif

      call fftw_free(a1_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Plan creation failed")

   end subroutine c2c_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in Y direction
   subroutine c2c_1m_y_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

#ifdef DOUBLE_PREC
      complex(C_DOUBLE_COMPLEX), contiguous, pointer :: a1(:, :)
      complex(C_DOUBLE_COMPLEX), contiguous, pointer :: a1o(:, :)
#else
      complex(C_FLOAT_COMPLEX), contiguous, pointer :: a1(:, :)
      complex(C_FLOAT_COMPLEX), contiguous, pointer :: a1o(:, :)
#endif
      type(C_PTR) :: a1_p
      integer(C_SIZE_T) :: sz

      ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
      ! done one Z-plane at a time. So plan for 2D data sets here.
      sz = decomp%ysz(1) * decomp%ysz(2)
      a1_p = fftw_alloc_complex(sz)
      call c_f_pointer(a1_p, a1, [decomp%ysz(1), decomp%ysz(2)])
      call c_f_pointer(a1_p, a1o, [decomp%ysz(1), decomp%ysz(2)])

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft(1, decomp%ysz(2), decomp%ysz(1), &
                                 a1, decomp%ysz(2), decomp%ysz(1), 1, a1o, decomp%ysz(2), &
                                 decomp%ysz(1), 1, isign, plan_type)
#else
      plan1 = fftwf_plan_many_dft(1, decomp%ysz(2), decomp%ysz(1), &
                                  a1, decomp%ysz(2), decomp%ysz(1), 1, a1o, decomp%ysz(2), &
                                  decomp%ysz(1), 1, isign, plan_type)
#endif

      call fftw_free(a1_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 2, "Plan creation failed")

   end subroutine c2c_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in Z direction
   subroutine c2c_1m_z_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

#ifdef DOUBLE_PREC
      complex(C_DOUBLE_COMPLEX), contiguous, pointer :: a1(:, :, :)
      complex(C_DOUBLE_COMPLEX), contiguous, pointer :: a1o(:, :, :)
#else
      complex(C_FLOAT_COMPLEX), contiguous, pointer :: a1(:, :, :)
      complex(C_FLOAT_COMPLEX), contiguous, pointer :: a1o(:, :, :)
#endif
      type(C_PTR) :: a1_p
      integer(C_SIZE_T) :: sz

      sz = decomp%zsz(1) * decomp%zsz(2) * decomp%zsz(3)
      a1_p = fftw_alloc_complex(sz)
      call c_f_pointer(a1_p, a1, [decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)])
      call c_f_pointer(a1_p, a1o, [decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)])

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft(1, decomp%zsz(3), &
                                 decomp%zsz(1) * decomp%zsz(2), a1, decomp%zsz(3), &
                                 decomp%zsz(1) * decomp%zsz(2), 1, a1o, decomp%zsz(3), &
                                 decomp%zsz(1) * decomp%zsz(2), 1, isign, plan_type)
#else
      plan1 = fftwf_plan_many_dft(1, decomp%zsz(3), &
                                  decomp%zsz(1) * decomp%zsz(2), a1, decomp%zsz(3), &
                                  decomp%zsz(1) * decomp%zsz(2), 1, a1o, decomp%zsz(3), &
                                  decomp%zsz(1) * decomp%zsz(2), 1, isign, plan_type)
#endif

      call fftw_free(a1_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 3, "Plan creation failed")

   end subroutine c2c_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D r2c FFTs in X direction
   subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

      real(mytype), contiguous, pointer :: a1(:, :, :)
      complex(mytype), contiguous, pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      a1_p = c_null_ptr
      a2_p = c_null_ptr

      sz = product(decomp_ph%xsz)
      if (inplace_r2c) sz = max(sz, 2 * product(int(decomp_sp%xsz, kind=C_SIZE_T)))
      a1_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp_ph%xsz)
      if (inplace_r2c) then
         call c_f_pointer(a1_p, a2, decomp_sp%xsz)
      else
         sz = product(decomp_sp%xsz)
         a2_p = fftw_alloc_complex(sz)
         call c_f_pointer(a2_p, a2, decomp_sp%xsz)
      end if

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft_r2c(1, decomp_ph%xsz(1), &
                                     decomp_ph%xsz(2) * decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
                                     decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
                                     plan_type)
#else
      plan1 = fftwf_plan_many_dft_r2c(1, decomp_ph%xsz(1), &
                                      decomp_ph%xsz(2) * decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
                                      decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
                                      plan_type)
#endif

      call fftw_free(a1_p)
      if (.not. inplace_r2c) call fftw_free(a2_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 4, "Plan creation failed")

   end subroutine r2c_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
   subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

      complex(mytype), contiguous, pointer :: a1(:, :, :)
      real(mytype), contiguous, pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      a1_p = c_null_ptr
      a2_p = c_null_ptr

      sz = product(decomp_ph%xsz)
      if (inplace_c2r) sz = max(sz, 2 * product(int(decomp_sp%xsz, kind=C_SIZE_T)))
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a2_p, a2, decomp_ph%xsz)
      if (inplace_c2r) then
         call c_f_pointer(a2_p, a1, decomp_sp%xsz)
      else
         sz = product(decomp_sp%xsz)
         a1_p = fftw_alloc_complex(sz)
         call c_f_pointer(a1_p, a1, decomp_sp%xsz)
      end if

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft_c2r(1, decomp_ph%xsz(1), &
                                     decomp_ph%xsz(2) * decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
                                     decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
                                     plan_type)
#else
      plan1 = fftwf_plan_many_dft_c2r(1, decomp_ph%xsz(1), &
                                      decomp_ph%xsz(2) * decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
                                      decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
                                      plan_type)
#endif

      if (.not. inplace_c2r) call fftw_free(a1_p)
      call fftw_free(a2_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 5, "Plan creation failed")

   end subroutine c2r_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
   subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

      real(mytype), contiguous, pointer :: a1(:, :, :)
      complex(mytype), contiguous, pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      a1_p = c_null_ptr
      a2_p = c_null_ptr

      sz = product(decomp_ph%zsz)
      if (inplace_r2c) sz = max(sz, 2 * product(int(decomp_sp%zsz, kind=C_SIZE_T)))
      a1_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp_ph%zsz)
      if (inplace_r2c) then
         call c_f_pointer(a1_p, a2, decomp_sp%zsz)
      else
         sz = product(decomp_sp%zsz)
         a2_p = fftw_alloc_complex(sz)
         call c_f_pointer(a2_p, a2, decomp_sp%zsz)
      end if

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft_r2c(1, decomp_ph%zsz(3), &
                                     decomp_ph%zsz(1) * decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
                                     decomp_ph%zsz(1) * decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
                                     decomp_sp%zsz(1) * decomp_sp%zsz(2), 1, plan_type)
#else
      plan1 = fftwf_plan_many_dft_r2c(1, decomp_ph%zsz(3), &
                                      decomp_ph%zsz(1) * decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
                                      decomp_ph%zsz(1) * decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
                                      decomp_sp%zsz(1) * decomp_sp%zsz(2), 1, plan_type)
#endif

      call fftw_free(a1_p)
      if (.not. inplace_r2c) call fftw_free(a2_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 6, "Plan creation failed")

   end subroutine r2c_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
   subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

      complex(mytype), contiguous, pointer :: a1(:, :, :)
      real(mytype), contiguous, pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      a1_p = c_null_ptr
      a2_p = c_null_ptr

      sz = product(decomp_ph%zsz)
      if (inplace_c2r) sz = max(sz, 2 * product(int(decomp_sp%zsz, kind=C_SIZE_T)))
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a2_p, a2, decomp_ph%zsz)
      if (inplace_c2r) then
         call c_f_pointer(a2_p, a1, decomp_sp%zsz)
      else
         sz = product(decomp_sp%zsz)
         a1_p = fftw_alloc_complex(sz)
         call c_f_pointer(a1_p, a1, decomp_sp%zsz)
      end if

#ifdef DOUBLE_PREC
      plan1 = fftw_plan_many_dft_c2r(1, decomp_ph%zsz(3), &
                                     decomp_ph%zsz(1) * decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
                                     decomp_sp%zsz(1) * decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
                                     decomp_ph%zsz(1) * decomp_ph%zsz(2), 1, plan_type)
#else
      plan1 = fftwf_plan_many_dft_c2r(1, decomp_ph%zsz(3), &
                                      decomp_ph%zsz(1) * decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
                                      decomp_sp%zsz(1) * decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
                                      decomp_ph%zsz(1) * decomp_ph%zsz(2), 1, plan_type)
#endif

      if (.not. inplace_c2r) call fftw_free(a1_p)
      call fftw_free(a2_p)
      if (.not. c_associated(plan1)) call decomp_2d_abort(__FILE__, __LINE__, 7, "Plan creation failed")

   end subroutine c2r_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in X direction, with possibility to dismiss points
   subroutine r2r_1m_x_plan(plan, decomp, dtt, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: dtt ! Type of DTT compatible with fftw3
      integer, intent(in) :: ndismiss ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :), a2(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :), a2(:, :, :)
#endif
      integer(C_INT) :: ntmp(1)
      integer(C_FFTW_R2R_KIND) :: tmp(1)

      call decomp_pool_get(a1, decomp%xsz)
      a2 => a1

      ntmp(1) = decomp%xsz(1) - ndismiss
      tmp(1) = dtt
#ifdef DOUBLE_PREC
      plan = fftw_plan_many_r2r(1, ntmp(1), decomp%xsz(2) * decomp%xsz(3), &
#else
      plan = fftwf_plan_many_r2r(1, ntmp(1), decomp%xsz(2) * decomp%xsz(3), & !&
#endif
                                 a1, decomp%xsz(1), 1, decomp%xsz(1), &
                                 a2, decomp%xsz(1), 1, decomp%xsz(1), &
                                 tmp(1), plan_type_dtt)

      call decomp_pool_free(a1)
      nullify (a2)
      if (.not. c_associated(plan)) call decomp_2d_abort(__FILE__, __LINE__, 17, "Plan creation failed")

   end subroutine r2r_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Y direction, with possibility to dismiss points
   subroutine r2r_1m_y_plan(plan, decomp, dtt, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: dtt ! Type of DTT compatible with fftw3
      integer, intent(in) :: ndismiss ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :), a2(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :), a2(:, :, :)
#endif
      integer(C_INT) :: ntmp(1)
      integer(C_FFTW_R2R_KIND) :: tmp(1)

      call decomp_pool_get(a1, (/decomp%ysz(1), decomp%ysz(2), 1/))
      a2 => a1

      ntmp(1) = decomp%ysz(2) - ndismiss
      tmp(1) = dtt
#ifdef DOUBLE_PREC
      plan = fftw_plan_many_r2r(1, ntmp(1), decomp%ysz(1), &
#else
      plan = fftwf_plan_many_r2r(1, ntmp(1), decomp%ysz(1), & !&
#endif
                                 a1, decomp%ysz(2), decomp%ysz(1), 1, &
                                 a2, decomp%ysz(2), decomp%ysz(1), 1, &
                                 tmp(1), plan_type_dtt)

      call decomp_pool_free(a1)
      nullify (a2)
      if (.not. c_associated(plan)) call decomp_2d_abort(__FILE__, __LINE__, 18, "Plan creation failed")

   end subroutine r2r_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Z direction, with possibility to dismiss points
   subroutine r2r_1m_z_plan(plan, decomp, dtt, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: dtt ! Type of DTT compatible with fftw3
      integer, intent(in) :: ndismiss ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :), a2(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :), a2(:, :, :)
#endif
      integer(C_INT) :: ntmp(1)
      integer(C_FFTW_R2R_KIND) :: tmp(1)

      call decomp_pool_get(a1, decomp%zsz)
      a2 => a1

      ntmp(1) = decomp%zsz(3) - ndismiss
      tmp(1) = dtt
#ifdef DOUBLE_PREC
      plan = fftw_plan_many_r2r(1, ntmp(1), decomp%zsz(1) * decomp%zsz(2), &
#else
      plan = fftwf_plan_many_r2r(1, ntmp(1), decomp%zsz(1) * decomp%zsz(2), & !&
#endif
                                 a1, decomp%zsz(3), decomp%zsz(1) * decomp%zsz(2), 1, &
                                 a2, decomp%zsz(3), decomp%zsz(1) * decomp%zsz(2), 1, &
                                 tmp(1), plan_type_dtt)

      call decomp_pool_free(a1)
      nullify (a2)
      if (.not. c_associated(plan)) call decomp_2d_abort(__FILE__, __LINE__, 19, "Plan creation failed")

   end subroutine r2r_1m_z_plan

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time initialisations for the FFT engine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_fft_engine

#ifdef FFTW_omp
      use omp_lib
#endif

      implicit none

#ifdef FFTW_omp
      integer :: ierr
#endif

      call decomp_2d_fft_log("FFTW (F2003 interface)")

#ifdef FFTW_omp
      ierr = fftw_init_threads()
      if (ierr==0) call decomp_2d_abort(__FILE__, __LINE__, ierr, "fftw_init_threads")
      call fftw_plan_with_nthreads(omp_get_max_threads())
#endif

      if (format == PHYSICAL_IN_X) then

         ! For C2C transforms
         call c2c_1m_x_plan(plan(-1, 1), ph, FFTW_FORWARD)
         call c2c_1m_y_plan(plan(-1, 2), ph, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(-1, 3), ph, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(1, 3), ph, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(1, 2), ph, FFTW_BACKWARD)
         call c2c_1m_x_plan(plan(1, 1), ph, FFTW_BACKWARD)

         ! For R2C/C2R tranforms
         call r2c_1m_x_plan(plan(0, 1), ph, sp)
         call c2c_1m_y_plan(plan(0, 2), sp, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(0, 3), sp, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(2, 3), sp, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(2, 2), sp, FFTW_BACKWARD)
         call c2r_1m_x_plan(plan(2, 1), sp, ph)

      else if (format == PHYSICAL_IN_Z) then

         ! For C2C transforms
         call c2c_1m_z_plan(plan(-1, 3), ph, FFTW_FORWARD)
         call c2c_1m_y_plan(plan(-1, 2), ph, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(-1, 1), ph, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(1, 1), ph, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(1, 2), ph, FFTW_BACKWARD)
         call c2c_1m_z_plan(plan(1, 3), ph, FFTW_BACKWARD)

         ! For R2C/C2R tranforms
         call r2c_1m_z_plan(plan(0, 3), ph, sp)
         call c2c_1m_y_plan(plan(0, 2), sp, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(0, 1), sp, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(2, 1), sp, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(2, 2), sp, FFTW_BACKWARD)
         call c2r_1m_z_plan(plan(2, 3), sp, ph)

      end if

   end subroutine init_fft_engine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time finalisations for the FFT engine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine finalize_fft_engine(local_plan)

      implicit none

      type(c_ptr) :: local_plan(-1:2, 3)

      integer :: i, j

      do j = 1, 3
         do i = -1, 2
#ifdef DOUBLE_PREC
            call fftw_destroy_plan(local_plan(i, j))
#else
            call fftwf_destroy_plan(local_plan(i, j))
#endif
         end do
      end do

   end subroutine finalize_fft_engine

   ! Following routines calculate multiple one-dimensional FFTs to form
   ! the basis of three-dimensional FFTs.

   ! c2c transform, multiple 1D FFTs in x direction
   subroutine c2c_1m_x(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR) :: plan1

      if (skip_x_c2c) return

#ifdef DOUBLE_PREC
      call fftw_execute_dft(plan1, inout, inout)
#else
      call fftwf_execute_dft(plan1, inout, inout)
#endif
   end subroutine c2c_1m_x

   ! c2c transform, multiple 1D FFTs in y direction
   subroutine c2c_1m_y(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR) :: plan1

      integer :: k, s3

      if (skip_y_c2c) return

      s3 = size(inout, 3)

      do k = 1, s3  ! transform on one Z-plane at a time
#ifdef DOUBLE_PREC
         call fftw_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
#else
         call fftwf_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
#endif
      end do
   end subroutine c2c_1m_y

   ! c2c transform, multiple 1D FFTs in z direction
   subroutine c2c_1m_z(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR) :: plan1

      if (skip_z_c2c) return

#ifdef DOUBLE_PREC
      call fftw_execute_dft(plan1, inout, inout)
#else
      call fftwf_execute_dft(plan1, inout, inout)
#endif
   end subroutine c2c_1m_z

   ! r2c transform, multiple 1D FFTs in x direction
   subroutine r2c_1m_x(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      complex(mytype), dimension(:, :, :), intent(INOUT) :: output

      if (skip_x_c2c) call decomp_2d_warning(__FILE__, __LINE__, 1, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      call fftw_execute_dft_r2c(plan(0, 1), input, output)
#else
      call fftwf_execute_dft_r2c(plan(0, 1), input, output)
#endif

      return

   end subroutine r2c_1m_x

   ! r2c transform, multiple 1D FFTs in z direction
   subroutine r2c_1m_z(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      complex(mytype), dimension(:, :, :), intent(INOUT) :: output

      if (skip_z_c2c) call decomp_2d_warning(__FILE__, __LINE__, 2, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      call fftw_execute_dft_r2c(plan(0, 3), input, output)
#else
      call fftwf_execute_dft_r2c(plan(0, 3), input, output)
#endif

      return

   end subroutine r2c_1m_z

   ! c2r transform, multiple 1D FFTs in x direction
   subroutine c2r_1m_x(input, output)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      real(mytype), dimension(:, :, :), intent(INOUT) :: output

      if (skip_x_c2c) call decomp_2d_warning(__FILE__, __LINE__, 3, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      call fftw_execute_dft_c2r(plan(2, 1), input, output)
#else
      call fftwf_execute_dft_c2r(plan(2, 1), input, output)
#endif

      return

   end subroutine c2r_1m_x

   ! c2r transform, multiple 1D FFTs in z direction
   subroutine c2r_1m_z(input, output)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: input
      real(mytype), dimension(:, :, :), intent(INOUT) :: output

      if (skip_z_c2c) call decomp_2d_warning(__FILE__, __LINE__, 4, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      call fftw_execute_dft_c2r(plan(2, 3), input, output)
#else
      call fftwf_execute_dft_c2r(plan(2, 3), input, output)
#endif

      return

   end subroutine c2r_1m_z

   ! r2r transform, multiple 1D DTTs in x direction
   subroutine r2r_1m_x(inr, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: ifirst

      ! Exit if needed
      if (skip_x_c2c) return

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(1)
         ifirst = dtt(4)
      else
         if (c_associated(dtt_plan(4))) then
            plan = dtt_plan(4)
         else
            plan = dtt_plan(1)
         end if
         ifirst = dtt(10)
      end if

      ! Perform the DFT
      call wrapper_r2r(plan, &
                       inr, ifirst, size(inr))

   end subroutine r2r_1m_x

   ! r2r transform, multiple 1D DTTs in y direction
   subroutine r2r_1m_y(inr, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: k, ifirst

      ! Exit if needed
      if (skip_y_c2c) return

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(2)
         ifirst = dtt(5)
      else
         if (c_associated(dtt_plan(5))) then
            plan = dtt_plan(5)
         else
            plan = dtt_plan(2)
         end if
         ifirst = dtt(11)
      end if

      ! Perform the DFT
      do k = 1, size(inr, 3)
         call wrapper_r2r(plan, &
                          inr(:, :, k:k), 1 + size(inr, 1) * (ifirst - 1), size(inr, 1) * size(inr, 2))
      end do

   end subroutine r2r_1m_y

   ! r2r transform, multiple 1D DTTs in z direction
   subroutine r2r_1m_z(inr, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: ifirst

      ! Exit if needed
      if (skip_z_c2c) return

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(3)
         ifirst = dtt(6)
      else
         if (c_associated(dtt_plan(6))) then
            plan = dtt_plan(6)
         else
            plan = dtt_plan(3)
         end if
         ifirst = dtt(12)
      end if

      ! Perform the DFT
      call wrapper_r2r(plan, &
                       inr, 1 + size(inr, 1) * size(inr, 2) * (ifirst - 1), size(inr))

   end subroutine r2r_1m_z

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D FFT - complex to complex
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2c(in, out, isign)

      implicit none

      ! Arguments
      complex(mytype), dimension(:, :, :), intent(INOUT) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out
      integer, intent(IN) :: isign

      ! Local variables
      complex(mytype), contiguous, pointer :: wk2_c2c(:, :, :), wk1(:, :, :)

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2c")

      ! Init pointers
      call decomp_pool_get(wk2_c2c, ph%ysz)
      nullify (wk1)

      if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
          format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

         ! ===== 1D FFTs in X =====
         if (inplace) then
            call c2c_1m_x(in, plan(isign, 1))
         else
            call decomp_pool_get(wk1, ph%xsz)
            wk1(:,:,:) = in(:,:,:)
            call c2c_1m_x(wk1, plan(isign, 1))
         end if

         ! ===== Swap X --> Y; 1D FFTs in Y =====

         if (dims(1) > 1) then
            if (inplace) then
               call transpose_x_to_y(in, wk2_c2c, ph)
            else
               call transpose_x_to_y(wk1, wk2_c2c, ph)
            end if
            call c2c_1m_y(wk2_c2c, plan(isign, 2))
         else
            if (inplace) then
               call c2c_1m_y(in, plan(isign, 2))
            else
               call c2c_1m_y(wk1, plan(isign, 2))
            end if
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_c2c, out, ph)
         else
            if (inplace) then
               call transpose_y_to_z(in, out, ph)
            else
               call transpose_y_to_z(wk1, out, ph)
            end if
         end if
         call c2c_1m_z(out, plan(isign, 3))

      else if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_BACKWARD &
               .OR. &
               format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_FORWARD) then

         ! ===== 1D FFTs in Z =====
         if (inplace) then
            call c2c_1m_z(in, plan(isign, 3))
         else
            call decomp_pool_get(wk1, ph%zsz)
            wk1(:,:,:) = in(:,:,:)
            call c2c_1m_z(wk1, plan(isign, 3))
         end if

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            if (inplace) then
               call transpose_z_to_y(in, wk2_c2c, ph)
            else
               call transpose_z_to_y(wk1, wk2_c2c, ph)
            end if
            call c2c_1m_y(wk2_c2c, plan(isign, 2))
         else  ! out==wk2_c2c if 1D decomposition
            if (inplace) then
               call transpose_z_to_y(in, out, ph)
            else
               call transpose_z_to_y(wk1, out, ph)
            end if
            call c2c_1m_y(out, plan(isign, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_c2c, out, ph)
         end if
         call c2c_1m_x(out, plan(isign, 1))

      end if

      ! Free memory
      call decomp_pool_free(wk2_c2c)
      if (associated(wk1)) call decomp_pool_free(wk1)

      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2c")

   end subroutine fft_3d_c2c

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D forward FFT - real to complex
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_r2c(in_r, out_c)

      implicit none

      ! Arguments
      real(mytype), target, contiguous, dimension(:, :, :), intent(INOUT) :: in_r
      complex(mytype), dimension(:, :, :), intent(OUT) :: out_c

      ! Local variable
      complex(mytype), contiguous, pointer :: wk2_r2c(:, :, :), wk0(:, :, :)

      if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")

      ! Init pointers
      call decomp_pool_get(wk2_r2c, sp%ysz)
      nullify (wk0)

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in X =====
         if (inplace_r2c) then
            call c_f_pointer(c_loc(in_r), wk0, sp%xsz)
         else
            call decomp_pool_get(wk0, sp%xsz)
         end if
         call r2c_1m_x(in_r, wk0)

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            call transpose_x_to_y(wk0, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, plan(0, 2))
         else
            call c2c_1m_y(wk0, plan(0, 2))
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, out_c, sp)
         else
            call transpose_y_to_z(wk0, out_c, sp)
         end if
         call c2c_1m_z(out_c, plan(0, 3))

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in Z =====
         if (inplace_r2c) then
            call c_f_pointer(c_loc(in_r), wk0, sp%zsz)
         else
            call decomp_pool_get(wk0, sp%zsz)
         end if
         call r2c_1m_z(in_r, wk0)

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            call transpose_z_to_y(wk0, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, plan(0, 2))
         else  ! out_c==wk2_r2c if 1D decomposition
            call transpose_z_to_y(wk0, out_c, sp)
            call c2c_1m_y(out_c, plan(0, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_r2c, out_c, sp)
         end if
         call c2c_1m_x(out_c, plan(0, 1))

      end if

      ! Free memory
      call decomp_pool_free(wk2_r2c)
      if (.not.inplace_r2c) call decomp_pool_free(wk0)

      if (decomp_profiler_fft) call decomp_profiler_end("fft_r2c")

      return
   end subroutine fft_3d_r2c

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D inverse FFT - complex to real
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2r(in_c, out_r)

      implicit none

      ! Arguments
      complex(mytype), dimension(:, :, :), intent(INOUT) :: in_c
      real(mytype), target, contiguous, dimension(:, :, :), intent(OUT) :: out_r

      ! Local variables
      complex(mytype), contiguous, pointer :: wk2_r2c(:, :, :), wk0(:, :, :), wk1(:, :, :)

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2r")

      ! Init pointers
      call decomp_pool_get(wk2_r2c, sp%ysz)
      nullify (wk0)
      nullify (wk1)

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in Z =====
         if (inplace) then
            call c2c_1m_z(in_c, plan(2, 3))
         else
            call decomp_pool_get(wk1, sp%zsz)
            wk1(:,:,:) = in_c(:,:,:)
            call c2c_1m_z(wk1, plan(2, 3))
         end if

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (inplace) then
            call transpose_z_to_y(in_c, wk2_r2c, sp)
         else
            call transpose_z_to_y(wk1, wk2_r2c, sp)
         end if
         call c2c_1m_y(wk2_r2c, plan(2, 2))

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (inplace_c2r) then
            call c_f_pointer(c_loc(out_r), wk0, sp%xsz)
         else
            call decomp_pool_get(wk0, sp%xsz)
         end if
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_r2c, wk0, sp)
            call c2r_1m_x(wk0, out_r)
         else
            call c2r_1m_x(wk2_r2c, out_r)
         end if

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in X =====
         if (inplace) then
            call c2c_1m_x(in_c, plan(2, 1))
         else
            call decomp_pool_get(wk1, sp%xsz)
            wk1(:,:,:) = in_c(:,:,:)
            call c2c_1m_x(wk1, plan(2, 1))
         end if

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            if (inplace) then
               call transpose_x_to_y(in_c, wk2_r2c, sp)
            else
               call transpose_x_to_y(wk1, wk2_r2c, sp)
            end if
            call c2c_1m_y(wk2_r2c, plan(2, 2))
         else  ! in_c==wk2_r2c if 1D decomposition
            if (inplace) then
               call c2c_1m_y(in_c, plan(2, 2))
            else
               call c2c_1m_y(wk1, plan(2, 2))
            end if
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (inplace_c2r) then
            call c_f_pointer(c_loc(out_r), wk0, sp%zsz)
         else
            call decomp_pool_get(wk0, sp%zsz)
         end if
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, wk0, sp)
         else
            if (inplace) then
               call transpose_y_to_z(in_c, wk0, sp)
            else
               call transpose_y_to_z(wk1, wk0, sp)
            end if
         end if
         call c2r_1m_z(wk0, out_r)

      end if

      ! Free memory
      call decomp_pool_free(wk2_r2c)
      if (.not.inplace_c2r) call decomp_pool_free(wk0)
      if (associated(wk1)) call decomp_pool_free(wk1)

      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2r")

   end subroutine fft_3d_c2r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Forward 3D DTT - real to real
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_2d_dtt_3d_r2r(in, out_real, isign)

      implicit none

      ! Arguments
      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: in
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: out_real
      integer, intent(in) :: isign

      ! Local variable
      real(mytype), dimension(:, :, :), contiguous, pointer :: wk2ra

      if (decomp_profiler_fft) call decomp_profiler_start("decomp_2d_dtt_3d_r2r")

      ! Safety check
      if (isign /= DECOMP_2D_FFT_FORWARD .and. isign /= DECOMP_2D_FFT_BACKWARD) then
         call decomp_2d_abort(__FILE__, __LINE__, isign, "Invalid value")
      end if

      ! Get a buffer for arrays in y-pencil                                                
      call decomp_pool_get(wk2ra, ph%ysz)

      ! Perform the 3D DTT
      if ((format == PHYSICAL_IN_X .and. isign == DECOMP_2D_FFT_FORWARD) .or. &
          (format == PHYSICAL_IN_Z .and. isign == DECOMP_2D_FFT_BACKWARD)) then

         ! DCT / DST in x
         call r2r_1m_x(in, isign)

         ! Transpose x => y
         call transpose_x_to_y(in, wk2ra, ph)

         ! DCT / DST in y
         call r2r_1m_y(wk2ra, isign)

         ! Transpose y => z
         call transpose_y_to_z(wk2ra, out_real, ph)

         ! DCT / DST in z
         call r2r_1m_z(out_real, isign)

      else

         ! DCT / DST in z
         call r2r_1m_z(in, isign)

         ! Transpose z => y
         call transpose_z_to_y(in, wk2ra, ph)

         ! DCT / DST in y
         call r2r_1m_y(wk2ra, isign)

         ! Transpose y => x
         call transpose_y_to_x(wk2ra, out_real, ph)

         ! DCT / DST in x
         call r2r_1m_x(out_real, isign)

      end if

      call decomp_pool_free(wk2ra)

      if (decomp_profiler_fft) call decomp_profiler_end("decomp_2d_dtt_3d_r2r")

   end subroutine decomp_2d_dtt_3d_r2r

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

   !
   ! Wrappers are needed for skipping points
   !
   subroutine wrapper_r2r(plan, &
                          inr, &
                          ii, isz)

      implicit none

      ! Arguments
      type(c_ptr), intent(in) :: plan
      real(mytype), dimension(:, :, :), target, contiguous, intent(inout) :: inr
      integer, intent(in) :: ii, isz

      ! Local variable
      real(mytype), dimension(:), contiguous, pointer :: inr2

      ! Create 1D pointers mapping the provided 3D array
      call c_f_pointer(c_loc(inr), inr2, (/isz/))

      ! Perform DFT starting at the ifirst location
#ifdef DOUBLE_PREC
      call fftw_execute_r2r(plan, inr2(ii:), inr2(ii:))
#else
      call fftwf_execute_r2r(plan, inr2(ii:), inr2(ii:))
#endif

      ! Release pointers
      nullify (inr2)

   end subroutine wrapper_r2r

end module decomp_2d_fft
