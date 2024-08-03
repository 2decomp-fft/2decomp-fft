!! SPDX-License-Identifier: BSD-3-Clause

! This is the FFTW implementation of the FFT library using
! the Fortran 2003 interface introduced in FFTW 3.3-beta1

module decomp_2d_fft

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_profiler
   use, intrinsic :: iso_c_binding

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

   ! Workspace to store the intermediate Y-pencil data
   ! *** TODO: investigate how to use only one workspace array
   complex(mytype), contiguous, pointer :: wk2_c2c(:, :, :) => null(), &
                                           wk2_r2c(:, :, :) => null(), &
                                           wk13(:, :, :) => null()

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

   ! Flag to check if a 1D DTT is periodic
   logical, save :: dtt_x_dft, dtt_y_dft, dtt_z_dft

   ! array with the DTT setup
   integer, contiguous, pointer, save :: dtt(:) => null()

   ! FFTW plans for DTT
   !    1, 2, 3 : forward transform in X, Y and Z
   !    4, 5, 6 : backward transform in X, Y and Z
   type(C_PTR), contiguous, pointer, save :: dtt_plan(:) => null()

   ! Decomposition objects for DTT
   type(decomp_info), pointer, save :: dtt_decomp_xy => null(), &
                                       dtt_decomp_yz => null(), &
                                       dtt_decomp_sp => null()

   ! Workspace for DTT
   real(mytype), contiguous, pointer :: wk1ra(:, :, :) => null(), &
                                        wk1ia(:, :, :) => null(), &
                                        wk1rb(:, :, :) => null(), &
                                        wk1ib(:, :, :) => null(), &
                                        wk2ra(:, :, :) => null(), &
                                        wk2ia(:, :, :) => null(), &
                                        wk2rb(:, :, :) => null(), &
                                        wk2ib(:, :, :) => null(), &
                                        wk3ra(:, :, :) => null(), &
                                        wk3ia(:, :, :) => null(), &
                                        wk3rb(:, :, :) => null(), &
                                        wk3ib(:, :, :) => null()

   ! Derived type with all the quantities needed to perform FFT
   type decomp_2d_fft_engine
      type(c_ptr), private :: plan(-1:2, 3), wk2_c2c_p, wk13_p
      integer, private :: format
      logical, private :: initialised = .false.
      integer, private :: nx_fft, ny_fft, nz_fft
      type(decomp_info), pointer, public :: ph => null()
      type(decomp_info), private :: ph_target ! ph => ph_target or ph => decomp_main
      type(decomp_info), public :: sp
      complex(mytype), contiguous, pointer, private :: wk2_c2c(:, :, :) => null()
      complex(mytype), contiguous, pointer, private :: wk2_r2c(:, :, :) => null()
      complex(mytype), contiguous, pointer, private :: wk13(:, :, :) => null()
      logical, private :: inplace, inplace_r2c, inplace_c2r
      logical, private :: skip_x_c2c, skip_y_c2c, skip_z_c2c
      ! Below is specific to DTT
      logical, public :: with_dtt
      integer, allocatable, public, dimension(:) :: dtt
      type(C_PTR), private :: dtt_plan(6)
      type(decomp_info), pointer, public :: dtt_decomp_sp => null()
      type(decomp_info), pointer, private :: dtt_decomp_xy => null(), &
                                             dtt_decomp_yz => null()
      type(decomp_info), private :: dtt_decomp_sp_target
      type(c_ptr), private :: wk1ra_p, wk1ia_p, wk1rb_p, wk1ib_p, &
                              wk2ra_p, wk2ia_p, wk2rb_p, wk2ib_p, &
                              wk3ra_p, wk3ia_p, wk3rb_p, wk3ib_p
      real(mytype), contiguous, pointer, private :: wk1ra(:, :, :) => null(), &
                                                    wk1ia(:, :, :) => null(), &
                                                    wk1rb(:, :, :) => null(), &
                                                    wk1ib(:, :, :) => null(), &
                                                    wk2ra(:, :, :) => null(), &
                                                    wk2ia(:, :, :) => null(), &
                                                    wk2rb(:, :, :) => null(), &
                                                    wk2ib(:, :, :) => null(), &
                                                    wk3ra(:, :, :) => null(), &
                                                    wk3ia(:, :, :) => null(), &
                                                    wk3rb(:, :, :) => null(), &
                                                    wk3ib(:, :, :) => null()
   contains
      procedure, public :: init => decomp_2d_fft_engine_init
      procedure, public :: fin => decomp_2d_fft_engine_fin
      procedure, public :: use_it => decomp_2d_fft_engine_use_it
      procedure, private :: dtt_init => decomp_2d_fft_engine_dtt_init
      procedure, private :: dtt_fin => decomp_2d_fft_engine_dtt_fin
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
             decomp_2d_dtt_3d_x2r, decomp_2d_dtt_3d_r2x, &
             decomp_2d_fft_finalize, decomp_2d_fft_get_size, &
             decomp_2d_fft_get_ph, decomp_2d_fft_get_sp, &
             decomp_2d_fft_get_dtt_sp, &
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

      integer(C_SIZE_T) :: sz

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
      ! Allocate the workspace for intermediate y-pencil data
      ! The largest memory block needed is the one for c2c transforms
      !
      sz = product(engine%ph%ysz)
      engine%wk2_c2c_p = fftw_alloc_complex(sz)
      call c_f_pointer(engine%wk2_c2c_p, engine%wk2_c2c, engine%ph%ysz)
      !
      ! A smaller memory block is needed for r2c and c2r transforms
      ! wk2_c2c and wk2_r2c start at the same location
      !
      !    Size of wk2_c2c : ph%ysz(1), ph%ysz(2), ph%ysz(3)
      !    Size of wk2_r2c : sp%ysz(1), sp%ysz(2), sp%ysz(3)
      !
      call c_f_pointer(engine%wk2_c2c_p, engine%wk2_r2c, engine%sp%ysz)
      !
      ! Allocate the workspace for r2c and c2r transforms
      !
      ! wk13 can not be easily fused with wk2_*2c due to statements such as
      ! transpose_y_to_x(wk2_r2c, wk13, sp)
      ! transpose_y_to_z(wk2_r2c, wk13, sp)
      !
      if (.not. (engine%inplace_r2c .and. engine%inplace_c2r)) then
         if (pencil == PHYSICAL_IN_X) then
            sz = product(engine%sp%xsz)
            engine%wk13_p = fftw_alloc_complex(sz)
            call c_f_pointer(engine%wk13_p, engine%wk13, engine%sp%xsz)
         else if (pencil == PHYSICAL_IN_Z) then
            sz = product(engine%sp%zsz)
            engine%wk13_p = fftw_alloc_complex(sz)
            call c_f_pointer(engine%wk13_p, engine%wk13, engine%sp%zsz)
         end if
      end if

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

      integer :: i
      integer(C_SIZE_T) :: sz

      ! Safety check
      if (minval(in_DTT) < 0) call decomp_2d_abort(__FILE__, __LINE__, minval(in_DTT), "Invalid argument")
      if (maxval(in_DTT) > 8) call decomp_2d_abort(__FILE__, __LINE__, maxval(in_DTT), "Invalid argument")

      ! Prepare engine%dtt
      ! Mandatory
      !    1:3 => type of forward transform
      ! Optional, default values in dtt_assign_default
      !    4:6 => ifirst, index where the transform starts in the input
      !    7:9 => ndismiss, number of points skipped
      !    10:12 => ofirst, index where the transform starts in the output
      ! Values defined in dtt_invert
      !    12:15 => type of backward transform
      allocate (engine%dtt(15))
      if (size(in_DTT) == 12) then
         engine%dtt(1:12) = in_DTT
      elseif (size(in_DTT) == 3) then
         engine%dtt(1:3) = in_DTT
         call dtt_assign_default(engine%dtt)
      else
         call decomp_2d_abort(__FILE__, __LINE__, size(in_DTT), "Invalid argument")
      end if
      call dtt_invert(engine%dtt)
      call dtt_for_fftw(engine%dtt(1:3))
      call dtt_for_fftw(engine%dtt(13:15))

      ! Prepare the decomp_info objects
      if (engine%format == PHYSICAL_IN_X) then
         if (engine%dtt(1) == FFTW_FORWARD) then
            engine%dtt_decomp_xy => engine%sp
            engine%dtt_decomp_yz => engine%sp
            engine%dtt_decomp_sp => engine%sp
         else if (engine%dtt(2) == FFTW_FORWARD) then
            engine%dtt_decomp_xy => engine%ph
            call decomp_info_init(engine%nx_fft, engine%ny_fft / 2 + 1, engine%nz_fft, engine%dtt_decomp_sp_target)
            engine%dtt_decomp_yz => engine%dtt_decomp_sp_target
            engine%dtt_decomp_sp => engine%dtt_decomp_sp_target
         else if (engine%dtt(3) == FFTW_FORWARD) then
            engine%dtt_decomp_xy => engine%ph
            engine%dtt_decomp_yz => engine%ph
            call decomp_info_init(engine%nx_fft, engine%ny_fft, engine%nz_fft / 2 + 1, engine%dtt_decomp_sp_target)
            engine%dtt_decomp_sp => engine%dtt_decomp_sp_target
         else
            engine%dtt_decomp_xy => engine%ph
            engine%dtt_decomp_yz => engine%ph
            engine%dtt_decomp_sp => engine%ph
         end if
      else
         if (engine%dtt(3) == FFTW_FORWARD) then
            engine%dtt_decomp_yz => engine%sp
            engine%dtt_decomp_xy => engine%sp
            engine%dtt_decomp_sp => engine%sp
         else if (engine%dtt(2) == FFTW_FORWARD) then
            engine%dtt_decomp_yz => engine%ph
            call decomp_info_init(engine%nx_fft, engine%ny_fft / 2 + 1, engine%nz_fft, engine%dtt_decomp_sp_target)
            engine%dtt_decomp_xy => engine%dtt_decomp_sp_target
            engine%dtt_decomp_sp => engine%dtt_decomp_sp_target
         else if (engine%dtt(1) == FFTW_FORWARD) then
            engine%dtt_decomp_yz => engine%ph
            engine%dtt_decomp_xy => engine%ph
            call decomp_info_init(engine%nx_fft / 2 + 1, engine%ny_fft, engine%nz_fft, engine%dtt_decomp_sp_target)
            engine%dtt_decomp_sp => engine%dtt_decomp_sp_target
         else
            engine%dtt_decomp_yz => engine%ph
            engine%dtt_decomp_xy => engine%ph
            engine%dtt_decomp_sp => engine%ph
         end if
      end if

      ! Working arrays in x
      if (engine%format == PHYSICAL_IN_X) then
         ! Allocate the memory
         sz = product(engine%ph%xsz)
         engine%wk1ra_p = fftw_alloc_real(sz)
         engine%wk1ia_p = fftw_alloc_real(sz)
         ! Associated the pointers
         call c_f_pointer(engine%wk1ra_p, engine%wk1ra, engine%ph%xsz)
         call c_f_pointer(engine%wk1ia_p, engine%wk1ia, engine%ph%xsz)
      else
         sz = product(engine%dtt_decomp_sp%xsz)
         engine%wk1ra_p = fftw_alloc_real(sz)
         engine%wk1ia_p = fftw_alloc_real(sz)
         ! Associated the pointers
         call c_f_pointer(engine%wk1ra_p, engine%wk1ra, engine%dtt_decomp_sp%xsz)
         call c_f_pointer(engine%wk1ia_p, engine%wk1ia, engine%dtt_decomp_sp%xsz)
      end if
      sz = product(engine%dtt_decomp_xy%xsz)
      ! Allocate the memory
      engine%wk1rb_p = fftw_alloc_real(sz)
      engine%wk1ib_p = fftw_alloc_real(sz)
      call c_f_pointer(engine%wk1rb_p, engine%wk1rb, engine%dtt_decomp_xy%xsz)
      call c_f_pointer(engine%wk1ib_p, engine%wk1ib, engine%dtt_decomp_xy%xsz)
      ! Set to zero
      engine%wk1ra = 0._mytype
      engine%wk1ia = 0._mytype
      engine%wk1rb = 0._mytype
      engine%wk1ib = 0._mytype

      ! Working arrays in y
      ! Allocate the memory
      sz = product(engine%dtt_decomp_xy%ysz)
      engine%wk2ra_p = fftw_alloc_real(sz)
      engine%wk2ia_p = fftw_alloc_real(sz)
      sz = product(engine%dtt_decomp_yz%ysz)
      engine%wk2rb_p = fftw_alloc_real(sz)
      engine%wk2ib_p = fftw_alloc_real(sz)
      ! Associated the pointers
      call c_f_pointer(engine%wk2ra_p, engine%wk2ra, engine%dtt_decomp_xy%ysz)
      call c_f_pointer(engine%wk2ia_p, engine%wk2ia, engine%dtt_decomp_xy%ysz)
      call c_f_pointer(engine%wk2rb_p, engine%wk2rb, engine%dtt_decomp_yz%ysz)
      call c_f_pointer(engine%wk2ib_p, engine%wk2ib, engine%dtt_decomp_yz%ysz)
      ! Set to zero
      engine%wk2ra = 0._mytype
      engine%wk2ia = 0._mytype
      engine%wk2rb = 0._mytype
      engine%wk2ib = 0._mytype

      ! Working arrays in z
      ! Allocate the memory
      sz = product(engine%dtt_decomp_yz%zsz)
      engine%wk3ra_p = fftw_alloc_real(sz)
      engine%wk3ia_p = fftw_alloc_real(sz)
      ! Associated the pointers
      call c_f_pointer(engine%wk3ra_p, engine%wk3ra, engine%dtt_decomp_yz%zsz)
      call c_f_pointer(engine%wk3ia_p, engine%wk3ia, engine%dtt_decomp_yz%zsz)
      if (engine%format == PHYSICAL_IN_X) then
         ! Allocate the memory
         sz = product(engine%dtt_decomp_sp%zsz)
         engine%wk3rb_p = fftw_alloc_real(sz)
         engine%wk3ib_p = fftw_alloc_real(sz)
         ! Associated the pointers
         call c_f_pointer(engine%wk3rb_p, engine%wk3rb, engine%dtt_decomp_sp%zsz)
         call c_f_pointer(engine%wk3ib_p, engine%wk3ib, engine%dtt_decomp_sp%zsz)
      else
         ! Allocate the memory
         sz = product(engine%ph%zsz)
         engine%wk3rb_p = fftw_alloc_real(sz)
         engine%wk3ib_p = fftw_alloc_real(sz)
         ! Associated the pointers
         call c_f_pointer(engine%wk3rb_p, engine%wk3rb, engine%ph%zsz)
         call c_f_pointer(engine%wk3ib_p, engine%wk3ib, engine%ph%zsz)
      end if
      ! Set to zero
      engine%wk3ra = 0._mytype
      engine%wk3ia = 0._mytype
      engine%wk3rb = 0._mytype
      engine%wk3ib = 0._mytype

      ! Prepare the fftw plans
      engine%dtt_plan = c_null_ptr
      if (engine%format == PHYSICAL_IN_X) then
         ! in x
         if (engine%dtt(1) == FFTW_FORWARD) then
            call r2rr_1m_x_plan(engine%dtt_plan(1), engine%ph, engine%sp, engine%dtt(7))
            call rr2r_1m_x_plan(engine%dtt_plan(4), engine%sp, engine%ph, engine%dtt(7))
         else
            call r2r_1m_x_plan(engine%dtt_plan(1), engine%ph, engine%dtt(1), engine%dtt(7))
            call r2r_1m_x_plan(engine%dtt_plan(4), engine%ph, engine%dtt(13), engine%dtt(7))
         end if
         ! in y
         if (engine%dtt(2) == FFTW_FORWARD) then
            if (engine%dtt(1) == FFTW_FORWARD) then
               call rr2rr_1m_y_plan(engine%dtt_plan(2), engine%dtt_decomp_xy, engine%dtt(8))
               call rr2rr_1m_y_plan(engine%dtt_plan(5), engine%dtt_decomp_xy, engine%dtt(8)) ! Duplicated plan, this can be improved
            else
               call r2rr_1m_y_plan(engine%dtt_plan(2), engine%dtt_decomp_xy, engine%dtt_decomp_yz, engine%dtt(8))
               call rr2r_1m_y_plan(engine%dtt_plan(5), engine%dtt_decomp_yz, engine%dtt_decomp_xy, engine%dtt(8))
            end if
         else
            call r2r_1m_y_plan(engine%dtt_plan(2), engine%dtt_decomp_xy, engine%dtt(2), engine%dtt(8))
            call r2r_1m_y_plan(engine%dtt_plan(5), engine%dtt_decomp_xy, engine%dtt(14), engine%dtt(8))
         end if
         ! in z
         if (engine%dtt(3) == FFTW_FORWARD) then
            if (engine%dtt(1) == FFTW_FORWARD .or. engine%dtt(2) == FFTW_FORWARD) then
               call rr2rr_1m_z_plan(engine%dtt_plan(3), engine%dtt_decomp_sp, engine%dtt(9))
               call rr2rr_1m_z_plan(engine%dtt_plan(6), engine%dtt_decomp_sp, engine%dtt(9)) ! Duplicated plan, this can be improved
            else
               call r2rr_1m_z_plan(engine%dtt_plan(3), engine%dtt_decomp_yz, engine%dtt_decomp_sp, engine%dtt(9))
               call rr2r_1m_z_plan(engine%dtt_plan(6), engine%dtt_decomp_sp, engine%dtt_decomp_yz, engine%dtt(9))
            end if
         else
            call r2r_1m_z_plan(engine%dtt_plan(3), engine%dtt_decomp_sp, engine%dtt(3), engine%dtt(9))
            call r2r_1m_z_plan(engine%dtt_plan(6), engine%dtt_decomp_sp, engine%dtt(15), engine%dtt(9))
         end if
      else
         ! in z
         if (engine%dtt(3) == FFTW_FORWARD) then
            call r2rr_1m_z_plan(engine%dtt_plan(3), engine%ph, engine%sp, engine%dtt(9))
            call rr2r_1m_z_plan(engine%dtt_plan(6), engine%sp, engine%ph, engine%dtt(9))
         else
            call r2r_1m_z_plan(engine%dtt_plan(3), engine%ph, engine%dtt(3), engine%dtt(9))
            call r2r_1m_z_plan(engine%dtt_plan(6), engine%ph, engine%dtt(15), engine%dtt(9))
         end if
         ! in y
         if (engine%dtt(2) == FFTW_FORWARD) then
            if (engine%dtt(3) == FFTW_FORWARD) then
               call rr2rr_1m_y_plan(engine%dtt_plan(2), engine%dtt_decomp_xy, engine%dtt(8))
               call rr2rr_1m_y_plan(engine%dtt_plan(5), engine%dtt_decomp_xy, engine%dtt(8)) ! Duplicated plan, this can be improved
            else
               call r2rr_1m_y_plan(engine%dtt_plan(2), engine%dtt_decomp_yz, engine%dtt_decomp_xy, engine%dtt(8))
               call rr2r_1m_y_plan(engine%dtt_plan(5), engine%dtt_decomp_xy, engine%dtt_decomp_yz, engine%dtt(8))
            end if
         else
            call r2r_1m_y_plan(engine%dtt_plan(2), engine%dtt_decomp_xy, engine%dtt(2), engine%dtt(8))
            call r2r_1m_y_plan(engine%dtt_plan(5), engine%dtt_decomp_xy, engine%dtt(14), engine%dtt(8))
         end if
         ! in x
         if (engine%dtt(1) == FFTW_FORWARD) then
            if (engine%dtt(2) == FFTW_FORWARD .or. engine%dtt(3) == FFTW_FORWARD) then
               call rr2rr_1m_x_plan(engine%dtt_plan(1), engine%dtt_decomp_sp, engine%dtt(7))
               call rr2rr_1m_x_plan(engine%dtt_plan(4), engine%dtt_decomp_sp, engine%dtt(7)) ! Duplicated plan, this can be improved
            else
               call r2rr_1m_x_plan(engine%dtt_plan(1), engine%dtt_decomp_xy, engine%dtt_decomp_sp, engine%dtt(7))
               call rr2r_1m_x_plan(engine%dtt_plan(4), engine%dtt_decomp_sp, engine%dtt_decomp_xy, engine%dtt(7))
            end if
         else
            call r2r_1m_x_plan(engine%dtt_plan(1), engine%dtt_decomp_sp, engine%dtt(1), engine%dtt(7))
            call r2r_1m_x_plan(engine%dtt_plan(4), engine%dtt_decomp_sp, engine%dtt(13), engine%dtt(7))
         end if
      end if
      do i = 1, 6
         if (.not. c_associated(engine%dtt_plan(i))) call decomp_2d_abort(__FILE__, __LINE__, i, "DTT plan creation failed")
      end do

   end subroutine decomp_2d_fft_engine_dtt_init

   ! Set default values in the DTT config
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
         end select
      end do

   end subroutine dtt_for_fftw

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
      nullify (wk2_c2c)
      nullify (wk2_r2c)
      nullify (wk13)
      nullify (inplace)
      nullify (inplace_r2c)
      nullify (inplace_c2r)
      nullify (skip_x_c2c)
      nullify (skip_y_c2c)
      nullify (skip_z_c2c)
      nullify (with_dtt)
      if (associated(dtt_decomp_xy)) then
         nullify (dtt)
         nullify (dtt_plan)
         nullify (dtt_decomp_xy)
         nullify (dtt_decomp_yz)
         nullify (dtt_decomp_sp)
         nullify (wk1ra)
         nullify (wk1ia)
         nullify (wk1rb)
         nullify (wk1ib)
         nullify (wk2ra)
         nullify (wk2ia)
         nullify (wk2rb)
         nullify (wk2ib)
         nullify (wk3ra)
         nullify (wk3ia)
         nullify (wk3rb)
         nullify (wk3ib)
      end if

      ! Clean the FFTW library
      call fftw_cleanup()

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

      call fftw_free(engine%wk2_c2c_p)
      nullify (engine%wk2_c2c)
      nullify (engine%wk2_r2c)

      if (associated(wk13)) then
         call fftw_free(engine%wk13_p)
         nullify (engine%wk13)
      end if

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

      ! Clean the decomp_info objects
      nullify (engine%dtt_decomp_sp)
      nullify (engine%dtt_decomp_xy)
      nullify (engine%dtt_decomp_yz)
      if (allocated(engine%dtt_decomp_sp_target%x1dist)) &
         call decomp_info_finalize(engine%dtt_decomp_sp_target)

      ! Release the memory
      call fftw_free(engine%wk1ra_p)
      call fftw_free(engine%wk1ia_p)
      call fftw_free(engine%wk1rb_p)
      call fftw_free(engine%wk1ib_p)
      call fftw_free(engine%wk2ra_p)
      call fftw_free(engine%wk2ia_p)
      call fftw_free(engine%wk2rb_p)
      call fftw_free(engine%wk2ib_p)
      call fftw_free(engine%wk3ra_p)
      call fftw_free(engine%wk3ia_p)
      call fftw_free(engine%wk3rb_p)
      call fftw_free(engine%wk3ib_p)

      ! Set pointers to null
      nullify (engine%wk1ra)
      nullify (engine%wk1ia)
      nullify (engine%wk1rb)
      nullify (engine%wk1ib)
      nullify (engine%wk2ra)
      nullify (engine%wk2ia)
      nullify (engine%wk2rb)
      nullify (engine%wk2ib)
      nullify (engine%wk3ra)
      nullify (engine%wk3ia)
      nullify (engine%wk3rb)
      nullify (engine%wk3ib)

      ! Clean the fftw plans
      do i = 1, 6
#ifdef DOUBLE_PREC
         call fftw_destroy_plan(engine%dtt_plan(i))
#else
         call fftwf_destroy_plan(engine%dtt_plan(i))
#endif
      end do
      engine%dtt_plan = c_null_ptr

   end subroutine decomp_2d_fft_engine_dtt_fin

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
   ! Return a pointer to the decomp_info object sp used in DTT
   !
   ! The caller should not apply decomp_info_finalize on the pointer
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function decomp_2d_fft_get_dtt_sp()

      implicit none

      type(decomp_info), pointer :: decomp_2d_fft_get_dtt_sp

      if (.not. associated(dtt_decomp_sp)) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, 'No DTT available')
      end if
      decomp_2d_fft_get_dtt_sp => dtt_decomp_sp

   end function decomp_2d_fft_get_dtt_sp

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
      wk2_c2c => engine%wk2_c2c
      wk2_r2c => engine%wk2_r2c
      wk13 => engine%wk13
      inplace => engine%inplace
      inplace_r2c => engine%inplace_r2c
      inplace_c2r => engine%inplace_c2r
      skip_x_c2c => engine%skip_x_c2c
      skip_y_c2c => engine%skip_y_c2c
      skip_z_c2c => engine%skip_z_c2c
      with_dtt => engine%with_dtt
      dtt_x_dft = .false.
      dtt_y_dft = .false.
      dtt_z_dft = .false.

      if (with_dtt) then
         if (engine%dtt(1) == FFTW_FORWARD) dtt_x_dft = .true.
         if (engine%dtt(2) == FFTW_FORWARD) dtt_y_dft = .true.
         if (engine%dtt(3) == FFTW_FORWARD) dtt_z_dft = .true.
         dtt => engine%dtt
         dtt_plan => engine%dtt_plan
         dtt_decomp_xy => engine%dtt_decomp_xy
         dtt_decomp_yz => engine%dtt_decomp_yz
         dtt_decomp_sp => engine%dtt_decomp_sp
         wk1ra => engine%wk1ra
         wk1ia => engine%wk1ia
         wk1rb => engine%wk1rb
         wk1ib => engine%wk1ib
         wk2ra => engine%wk2ra
         wk2ia => engine%wk2ia
         wk2rb => engine%wk2rb
         wk2ib => engine%wk2ib
         wk3ra => engine%wk3ra
         wk3ia => engine%wk3ia
         wk3rb => engine%wk3rb
         wk3ib => engine%wk3ib
      else if (associated(dtt_decomp_xy)) then
         nullify (dtt)
         nullify (dtt_plan)
         nullify (dtt_decomp_xy)
         nullify (dtt_decomp_yz)
         nullify (dtt_decomp_sp)
         nullify (wk1ra)
         nullify (wk1ia)
         nullify (wk1rb)
         nullify (wk1ib)
         nullify (wk2ra)
         nullify (wk2ia)
         nullify (wk2rb)
         nullify (wk2ib)
         nullify (wk3ra)
         nullify (wk3ia)
         nullify (wk3rb)
         nullify (wk3ib)
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

      return
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

      return
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

      return
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

   end subroutine c2r_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in X direction, with possibility to dismiss points
   subroutine rr2rr_1m_x_plan(plan, decomp, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a4(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a4(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p, a4_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = product(decomp%xsz)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      a3_p = fftw_alloc_real(sz)
      a4_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%xsz)
      call c_f_pointer(a2_p, a2, decomp%xsz)
      call c_f_pointer(a3_p, a3, decomp%xsz)
      call c_f_pointer(a4_p, a4, decomp%xsz)

      dims(1)%n = decomp%xsz(1) - ndismiss
      dims(1)%is = 1
      dims(1)%os = 1
      howmany(1)%n = decomp%xsz(2) * decomp%xsz(3)
      howmany(1)%is = decomp%xsz(1)
      howmany(1)%os = decomp%xsz(1)
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft(1, dims(1), 1, howmany(1), a1, a2, a3, a4, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft(1, dims(1), 1, howmany(1), a1, a2, a3, a4, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      call fftw_free(a4_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)
      nullify (a4)

   end subroutine rr2rr_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Y direction, with possibility to dismiss points
   subroutine rr2rr_1m_y_plan(plan, decomp, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :)
      real(C_DOUBLE), contiguous, pointer :: a4(:, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :)
      real(C_FLOAT), contiguous, pointer :: a4(:, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p, a4_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = decomp%ysz(1) * decomp%ysz(2)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      a3_p = fftw_alloc_real(sz)
      a4_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a2_p, a2, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a3_p, a3, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a4_p, a4, (/decomp%ysz(1), decomp%ysz(2)/))

      dims(1)%n = decomp%ysz(2) - ndismiss
      dims(1)%is = decomp%ysz(1)
      dims(1)%os = decomp%ysz(1)
      howmany(1)%n = decomp%ysz(1)
      howmany(1)%is = 1
      howmany(1)%os = 1
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft(1, dims(1), 1, howmany(1), a1, a2, a3, a4, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft(1, dims(1), 1, howmany(1), a1, a2, a3, a4, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      call fftw_free(a4_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)
      nullify (a4)

   end subroutine rr2rr_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Z direction, with possibility to dismiss points
   subroutine rr2rr_1m_z_plan(plan, decomp, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a4(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a4(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p, a4_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = product(decomp%zsz)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      a3_p = fftw_alloc_real(sz)
      a4_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%zsz)
      call c_f_pointer(a2_p, a2, decomp%zsz)
      call c_f_pointer(a3_p, a3, decomp%zsz)
      call c_f_pointer(a4_p, a4, decomp%zsz)

      dims(1)%n = decomp%zsz(3) - ndismiss
      dims(1)%is = decomp%zsz(1) * decomp%zsz(2)
      dims(1)%os = decomp%zsz(1) * decomp%zsz(2)
      howmany(1)%n = decomp%zsz(1) * decomp%zsz(2)
      howmany(1)%is = 1
      howmany(1)%os = 1
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft(1, dims(1), 1, howmany(1), a1, a2, a3, a4, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft(1, dims(1), 1, howmany(1), a1, a2, a3, a4, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      call fftw_free(a4_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)
      nullify (a4)

   end subroutine rr2rr_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in X direction, with possibility to dismiss points
   subroutine rr2r_1m_x_plan(plan, decomp, decomp2, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp, decomp2
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = product(decomp%xsz)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      sz = product(decomp2%xsz)
      a3_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%xsz)
      call c_f_pointer(a2_p, a2, decomp%xsz)
      call c_f_pointer(a3_p, a3, decomp2%xsz)

      dims(1)%n = decomp2%xsz(1) - ndismiss
      dims(1)%is = 1
      dims(1)%os = 1
      howmany(1)%n = decomp%xsz(2) * decomp%xsz(3)
      howmany(1)%is = decomp%xsz(1)
      howmany(1)%os = decomp2%xsz(1)
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft_c2r(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft_c2r(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)

   end subroutine rr2r_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Y direction, with possibility to dismiss points
   subroutine rr2r_1m_y_plan(plan, decomp, decomp2, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp, decomp2
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = decomp%ysz(1) * decomp%ysz(2)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      sz = decomp2%ysz(1) * decomp2%ysz(2)
      a3_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a2_p, a2, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a3_p, a3, (/decomp2%ysz(1), decomp2%ysz(2)/))

      dims(1)%n = decomp2%ysz(2) - ndismiss
      dims(1)%is = decomp%ysz(1)
      dims(1)%os = decomp2%ysz(1)
      howmany(1)%n = decomp%ysz(1)
      howmany(1)%is = 1
      howmany(1)%os = 1
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft_c2r(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft_c2r(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)

   end subroutine rr2r_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Z direction, with possibility to dismiss points
   subroutine rr2r_1m_z_plan(plan, decomp, decomp2, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp, decomp2
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = product(decomp%zsz)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      sz = product(decomp2%zsz)
      a3_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%zsz)
      call c_f_pointer(a2_p, a2, decomp%zsz)
      call c_f_pointer(a3_p, a3, decomp2%zsz)

      dims(1)%n = decomp2%zsz(3) - ndismiss
      dims(1)%is = decomp%zsz(1) * decomp%zsz(2)
      dims(1)%os = decomp2%zsz(1) * decomp2%zsz(2)
      howmany(1)%n = decomp%zsz(1) * decomp%zsz(2)
      howmany(1)%is = 1
      howmany(1)%os = 1
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft_c2r(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft_c2r(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)

   end subroutine rr2r_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in X direction, with possibility to dismiss points
   subroutine r2rr_1m_x_plan(plan, decomp, decomp2, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp, decomp2
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = product(decomp%xsz)
      a1_p = fftw_alloc_real(sz)
      sz = product(decomp2%xsz)
      a2_p = fftw_alloc_real(sz)
      a3_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%xsz)
      call c_f_pointer(a2_p, a2, decomp2%xsz)
      call c_f_pointer(a3_p, a3, decomp2%xsz)

      dims(1)%n = decomp%xsz(1) - ndismiss
      dims(1)%is = 1
      dims(1)%os = 1
      howmany(1)%n = decomp%xsz(2) * decomp%xsz(3)
      howmany(1)%is = decomp%xsz(1)
      howmany(1)%os = decomp2%xsz(1)
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft_r2c(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft_r2c(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)

   end subroutine r2rr_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Y direction, with possibility to dismiss points
   subroutine r2rr_1m_y_plan(plan, decomp, decomp2, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp, decomp2
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = decomp%ysz(1) * decomp%ysz(2)
      a1_p = fftw_alloc_real(sz)
      sz = decomp2%ysz(1) * decomp2%ysz(2)
      a2_p = fftw_alloc_real(sz)
      a3_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a2_p, a2, (/decomp2%ysz(1), decomp2%ysz(2)/))
      call c_f_pointer(a3_p, a3, (/decomp2%ysz(1), decomp2%ysz(2)/))

      dims(1)%n = decomp%ysz(2) - ndismiss
      dims(1)%is = decomp%ysz(1)
      dims(1)%os = decomp2%ysz(1)
      howmany(1)%n = decomp%ysz(1)
      howmany(1)%is = 1
      howmany(1)%os = 1
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft_r2c(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft_r2c(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)

   end subroutine r2rr_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Z direction, with possibility to dismiss points
   subroutine r2rr_1m_z_plan(plan, decomp, decomp2, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp, decomp2
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a3(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a3(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p, a3_p
      integer(C_SIZE_T) :: sz
#ifdef DOUBLE_PREC
      type(fftw_iodim) :: dims(1), howmany(1)
#else
      type(fftwf_iodim) :: dims(1), howmany(1)
#endif

      sz = product(decomp%zsz)
      a1_p = fftw_alloc_real(sz)
      sz = product(decomp2%zsz)
      a2_p = fftw_alloc_real(sz)
      a3_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%zsz)
      call c_f_pointer(a2_p, a2, decomp2%zsz)
      call c_f_pointer(a3_p, a3, decomp2%zsz)

      dims(1)%n = decomp%zsz(3) - ndismiss
      dims(1)%is = decomp%zsz(1) * decomp%zsz(2)
      dims(1)%os = decomp2%zsz(1) * decomp2%zsz(2)
      howmany(1)%n = decomp%zsz(1) * decomp%zsz(2)
      howmany(1)%is = 1
      howmany(1)%os = 1
#ifdef DOUBLE_PREC
      plan = fftw_plan_guru_split_dft_r2c(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#else
      plan = fftwf_plan_guru_split_dft_r2c(1, dims(1), 1, howmany(1), a1, a2, a3, plan_type_dtt)
#endif

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      call fftw_free(a3_p)
      nullify (a1)
      nullify (a2)
      nullify (a3)

   end subroutine r2rr_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in X direction, with possibility to dismiss points
   subroutine r2r_1m_x_plan(plan, decomp, dtt, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: dtt
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz
      integer(C_INT) :: ntmp(1)
      integer(C_FFTW_R2R_KIND) :: tmp(1)

      sz = product(decomp%xsz)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%xsz)
      call c_f_pointer(a2_p, a2, decomp%xsz)

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

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      nullify (a1)
      nullify (a2)

   end subroutine r2r_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Y direction, with possibility to dismiss points
   subroutine r2r_1m_y_plan(plan, decomp, dtt, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: dtt
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :)
#endif
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz
      integer(C_INT) :: ntmp(1)
      integer(C_FFTW_R2R_KIND) :: tmp(1)

      sz = decomp%ysz(1) * decomp%ysz(2)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, (/decomp%ysz(1), decomp%ysz(2)/))
      call c_f_pointer(a2_p, a2, (/decomp%ysz(1), decomp%ysz(2)/))

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

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      nullify (a1)
      nullify (a2)

   end subroutine r2r_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D DTTs in Z direction, with possibility to dismiss points
   subroutine r2r_1m_z_plan(plan, decomp, dtt, ndismiss)

      implicit none

      type(C_PTR), intent(out) :: plan
      TYPE(DECOMP_INFO), intent(in) :: decomp
      integer, intent(in) :: dtt
      integer, intent(in) :: ndismiss  ! to dismiss n points from the signal

      ! Local variables
#ifdef DOUBLE_PREC
      real(C_DOUBLE), contiguous, pointer :: a1(:, :, :)
      real(C_DOUBLE), contiguous, pointer :: a2(:, :, :)
#else
      real(C_FLOAT), contiguous, pointer :: a1(:, :, :)
      real(C_FLOAT), contiguous, pointer :: a2(:, :, :)
#endif
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz
      integer(C_INT) :: ntmp(1)
      integer(C_FFTW_R2R_KIND) :: tmp(1)

      sz = product(decomp%zsz)
      a1_p = fftw_alloc_real(sz)
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, decomp%zsz)
      call c_f_pointer(a2_p, a2, decomp%zsz)

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

      call fftw_free(a1_p)
      call fftw_free(a2_p)
      nullify (a1)
      nullify (a2)

   end subroutine r2r_1m_z_plan

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time initialisations for the FFT engine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_fft_engine

      implicit none

      call decomp_2d_fft_log("FFTW (F2003 interface)")

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

   ! rr2rr transform, x direction
   subroutine rr2rr_1m_x(inr, ini, outr, outi, isign)

      implicit none

      ! Arguments
      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr, outi
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: ifirst, ofirst

      ! Copy and exit if needed
      if (skip_x_c2c) then
         outr = inr
         outi = ini
         return
      end if

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(1)
         ifirst = dtt(4)
         ofirst = dtt(10)
      else
         plan = dtt_plan(4)
         ifirst = dtt(10)
         ofirst = dtt(4)
      end if

      ! Perform the DFT
      if (isign == DECOMP_2D_FFT_FORWARD) then
         call wrapper_rr2rr(plan, &
                            inr, ini, ifirst, 1, 1, size(inr) - ifirst + 1, &
                            outr, outi, ofirst, 1, 1, size(outr) - ofirst + 1)
      else
         call wrapper_rr2rr(plan, &
                            ini, inr, ifirst, 1, 1, size(inr) - ifirst + 1, &
                            outi, outr, ofirst, 1, 1, size(outr) - ofirst + 1)
      end if

   end subroutine rr2rr_1m_x

   ! rr2rr transform, y direction
   subroutine rr2rr_1m_y(inr, ini, outr, outi, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr, outi
      integer, intent(in) :: isign

      ! Local variables
      integer :: k, ifirst, ofirst
      type(c_ptr) :: plan

      ! Copy and exit if needed
      if (skip_y_c2c) then
         outr = inr
         outi = ini
         return
      end if

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(2)
         ifirst = dtt(5)
         ofirst = dtt(11)
      else
         plan = dtt_plan(5)
         ifirst = dtt(11)
         ofirst = dtt(5)
      end if

      ! Perform the DFT
      if (isign == DECOMP_2D_FFT_FORWARD) then
         do k = 1, size(inr, 3)
            call wrapper_rr2rr(plan, &
                               inr, ini, 1, ifirst, k, size(inr, 1) * (size(inr, 2) - ifirst + 1), &
                               outr, outi, 1, ofirst, k, size(outr, 1) * (size(outr, 2) - ofirst + 1))
         end do
      else
         do k = 1, size(inr, 3)
            call wrapper_rr2rr(plan, &
                               ini, inr, 1, ifirst, k, size(inr, 1) * (size(inr, 2) - ifirst + 1), &
                               outi, outr, 1, ofirst, k, size(outr, 1) * (size(outr, 2) - ofirst + 1))
         end do
      end if

   end subroutine rr2rr_1m_y

   ! rr2rr transform, z direction
   subroutine rr2rr_1m_z(inr, ini, outr, outi, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr, outi
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: ifirst, ofirst

      ! Copy and exit if needed
      if (skip_z_c2c) then
         outr = inr
         outi = ini
         return
      end if

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(3)
         ifirst = dtt(6)
         ofirst = dtt(12)
      else
         plan = dtt_plan(6)
         ifirst = dtt(12)
         ofirst = dtt(6)
      end if

      ! Perform the DFT
      if (isign == DECOMP_2D_FFT_FORWARD) then
         call wrapper_rr2rr(plan, &
                            inr, ini, 1, 1, ifirst, size(inr, 1) * size(inr, 2) * (size(inr, 3) - ifirst + 1), &
                            outr, outi, 1, 1, ofirst, size(outr, 1) * size(outr, 2) * (size(outr, 3) - ofirst + 1))
      else
         call wrapper_rr2rr(plan, &
                            ini, inr, 1, 1, ifirst, size(inr, 1) * size(inr, 2) * (size(inr, 3) - ifirst + 1), &
                            outi, outr, 1, 1, ofirst, size(outr, 1) * size(outr, 2) * (size(outr, 3) - ofirst + 1))
      end if

   end subroutine rr2rr_1m_z

   ! r2r transform, x direction
   subroutine r2r_1m_x(inr, outr, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: ifirst, ofirst

      ! Copy and exit if needed
      if (skip_x_c2c) then
         outr = inr
         return
      end if

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(1)
         ifirst = dtt(4)
         ofirst = dtt(10)
      else
         plan = dtt_plan(4)
         ifirst = dtt(10)
         ofirst = dtt(4)
      end if

      ! Perform the DFT
      call wrapper_r2r(plan, &
                       inr, ifirst, 1, 1, size(inr) - ifirst + 1, &
                       outr, ofirst, 1, 1, size(outr) - ofirst + 1)

   end subroutine r2r_1m_x

   ! r2r transform, y direction
   subroutine r2r_1m_y(inr, outr, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: k, ifirst, ofirst

      ! Copy and exit if needed
      if (skip_y_c2c) then
         outr = inr
         return
      end if

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(2)
         ifirst = dtt(5)
         ofirst = dtt(11)
      else
         plan = dtt_plan(5)
         ifirst = dtt(11)
         ofirst = dtt(5)
      end if

      ! Perform the DFT
      do k = 1, size(inr, 3)
         call wrapper_r2r(plan, &
                          inr, 1, ifirst, k, size(inr, 1) * (size(inr, 2) - ifirst + 1), &
                          outr, 1, ofirst, k, size(outr, 1) * (size(outr, 2) - ofirst + 1))
      end do

   end subroutine r2r_1m_y

   ! r2r transform, z direction
   subroutine r2r_1m_z(inr, outr, isign)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr
      integer, intent(in) :: isign

      ! Local variables
      type(c_ptr) :: plan
      integer :: ifirst, ofirst

      ! Copy and exit if needed
      if (skip_z_c2c) then
         outr = inr
         return
      end if

      ! Get the DTT config
      if (isign == DECOMP_2D_FFT_FORWARD) then
         plan = dtt_plan(3)
         ifirst = dtt(6)
         ofirst = dtt(12)
      else
         plan = dtt_plan(6)
         ifirst = dtt(12)
         ofirst = dtt(6)
      end if

      ! Perform the DFT
      call wrapper_r2r(plan, &
                       inr, 1, 1, ifirst, size(inr, 1) * size(inr, 2) * (size(inr, 3) - ifirst + 1), &
                       outr, 1, 1, ofirst, size(outr, 1) * size(outr, 2) * (size(outr, 3) - ofirst + 1))

   end subroutine r2r_1m_z

   ! rr2r transform, x direction
   subroutine rr2r_1m_x(inr, ini, outr)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr

      ! Perform the DFT
      call wrapper_rr2r(dtt_plan(4), &
                        inr, ini, dtt(10), 1, 1, size(inr) - dtt(10) + 1, &
                        outr, dtt(4), 1, 1, size(outr) - dtt(4) + 1)

   end subroutine rr2r_1m_x

   ! rr2r transform, y direction
   subroutine rr2r_1m_y(inr, ini, outr)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr

      integer :: k

      ! Perform the DFT
      do k = 1, size(inr, 3)
         call wrapper_rr2r(dtt_plan(5), &
                           inr, ini, 1, dtt(11), k, size(inr, 1) * (size(inr, 2) - dtt(11) + 1), &
                           outr, 1, dtt(5), k, size(outr, 1) * (size(outr, 2) - dtt(5) + 1))
      end do

   end subroutine rr2r_1m_y

   ! rr2r transform, z direction
   subroutine rr2r_1m_z(inr, ini, outr)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr

      ! Perform the DFT
      call wrapper_rr2r(dtt_plan(6), &
                        inr, ini, 1, 1, dtt(12), size(inr, 1) * size(inr, 2) * (size(inr, 3) - dtt(12) + 1), &
                        outr, 1, 1, dtt(6), size(outr, 1) * size(outr, 2) * (size(outr, 3) - dtt(6) + 1))

   end subroutine rr2r_1m_z

   ! r2rr transform, x direction
   subroutine r2rr_1m_x(inr, outr, outi)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr, outi

      ! Perform the DFT
      call wrapper_r2rr(dtt_plan(1), &
                        inr, dtt(4), 1, 1, size(inr) - dtt(4) + 1, &
                        outr, outi, dtt(10), 1, 1, size(outr) - dtt(10) + 1)

   end subroutine r2rr_1m_x

   ! r2rr transform, y direction
   subroutine r2rr_1m_y(inr, outr, outi)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr, outi

      integer :: k

      ! Perform the DFT
      do k = 1, size(inr, 3)
         call wrapper_r2rr(dtt_plan(2), &
                           inr, 1, dtt(5), k, size(inr, 1) * (size(inr, 2) - dtt(5) + 1), &
                           outr, outi, 1, dtt(11), k, size(outr, 1) * (size(outr, 2) - dtt(11) + 1))
      end do

   end subroutine r2rr_1m_y

   ! r2rr transform, z direction
   subroutine r2rr_1m_z(inr, outr, outi)

      implicit none

      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: inr
      real(mytype), dimension(:, :, :), contiguous, target, intent(out) :: outr, outi

      ! Perform the DFT
      call wrapper_r2rr(dtt_plan(3), &
                        inr, 1, 1, dtt(6), size(inr, 1) * size(inr, 2) * (size(inr, 3) - dtt(6) + 1), &
                        outr, outi, 1, 1, dtt(12), size(outr, 1) * size(outr, 2) * (size(outr, 3) - dtt(12) + 1))

   end subroutine r2rr_1m_z

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
      complex(mytype), contiguous, pointer :: wk1(:, :, :)
      integer(C_SIZE_T) :: sz
      type(C_PTR) :: wk1_p

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2c")

      ! Initialise to NULL pointer
      nullify (wk1)
      wk1_p = c_null_ptr

      if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
          format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

         ! ===== 1D FFTs in X =====
         if (inplace) then
            call c2c_1m_x(in, plan(isign, 1))
         else
            sz = product(ph%xsz)
            wk1_p = fftw_alloc_complex(sz)
            call c_f_pointer(wk1_p, wk1, ph%xsz)
            wk1 = in
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
            sz = product(ph%zsz)
            wk1_p = fftw_alloc_complex(sz)
            call c_f_pointer(wk1_p, wk1, ph%zsz)
            wk1 = in
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
      if (associated(wk1)) then
         call fftw_free(wk1_p)
         nullify (wk1)
      end if

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
      complex(mytype), contiguous, pointer :: wk0(:, :, :)

      if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")

      ! Initialise to NULL pointer
      nullify (wk0)

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in X =====
         if (inplace_r2c) then
            call c_f_pointer(c_loc(in_r), wk0, sp%xsz)
         else
            wk0 => wk13
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
            wk0 => wk13
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
      if (associated(wk0)) nullify (wk0)

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
      complex(mytype), contiguous, pointer :: wk0(:, :, :), wk1(:, :, :)
      integer(C_SIZE_T) :: sz
      type(C_PTR) :: wk1_p

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2r")

      ! Initialise to NULL pointer
      nullify (wk0)
      nullify (wk1)
      wk1_p = c_null_ptr

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in Z =====
         if (inplace) then
            call c2c_1m_z(in_c, plan(2, 3))
         else
            sz = product(sp%zsz)
            wk1_p = fftw_alloc_complex(sz)
            call c_f_pointer(wk1_p, wk1, sp%zsz)
            wk1 = in_c
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
            wk0 => wk13
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
            sz = product(sp%xsz)
            wk1_p = fftw_alloc_complex(sz)
            call c_f_pointer(wk1_p, wk1, sp%xsz)
            wk1 = in_c
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
            wk0 => wk13
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
      if (associated(wk0)) nullify (wk0)
      if (associated(wk1)) then
         call fftw_free(wk1_p)
         nullify (wk1)
      end if

      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2r")

   end subroutine fft_3d_c2r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Forward 3D DTT - real to real/complex
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_2d_dtt_3d_r2x(in, out_real, out_cplx)

      implicit none

      ! Arguments
      real(mytype), dimension(:, :, :), contiguous, target, intent(inout) :: in
      real(mytype), dimension(:, :, :), contiguous, target, intent(out), optional :: out_real
      complex(mytype), dimension(:, :, :), contiguous, target, intent(out), optional :: out_cplx

      ! Local variables
      logical :: cplx
      integer :: i, j, k

      if (decomp_profiler_fft) call decomp_profiler_start("decomp_2d_dtt_3d_r2x")

      ! Safety check
      if (.not. (present(out_real) .or. present(out_cplx))) then
         call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
      end if

      ! Init
      cplx = .false.

      ! Perform the 3D DTT
      if (format == PHYSICAL_IN_X) then

         ! DFT / DTT in x
         if (dtt_x_dft) then
            call r2rr_1m_x(in, wk1rb, wk1ib)
            cplx = .true.
         else
            call r2r_1m_x(in, wk1rb, DECOMP_2D_FFT_FORWARD)
         end if

         ! Transpose x => y
         call transpose_x_to_y(wk1rb, wk2ra, dtt_decomp_xy)
         if (cplx) then
            call transpose_x_to_y(wk1ib, wk2ia, dtt_decomp_xy)
         end if

         ! DFT / DTT in y
         if (dtt_y_dft) then
            if (cplx) then
               call rr2rr_1m_y(wk2ra, wk2ia, wk2rb, wk2ib, DECOMP_2D_FFT_FORWARD)
            else
               call r2rr_1m_y(wk2ra, wk2rb, wk2ib)
               cplx = .true.
            end if
         else
            call r2r_1m_y(wk2ra, wk2rb, DECOMP_2D_FFT_FORWARD)
            if (cplx) then
               call r2r_1m_y(wk2ia, wk2ib, DECOMP_2D_FFT_FORWARD)
            end if
         end if

         ! Transpose y => z
         call transpose_y_to_z(wk2rb, wk3ra, dtt_decomp_yz)
         if (cplx) then
            call transpose_y_to_z(wk2ib, wk3ia, dtt_decomp_yz)
         end if

         ! DFT / DTT in z
         if (dtt_z_dft) then
            if (cplx) then
               call rr2rr_1m_z(wk3ra, wk3ia, wk3rb, wk3ib, DECOMP_2D_FFT_FORWARD)
            else
               call r2rr_1m_z(wk3ra, wk3rb, wk3ib)
               cplx = .true.
            end if
         else
            if (cplx) then
               call r2r_1m_z(wk3ra, wk3rb, DECOMP_2D_FFT_FORWARD)
               call r2r_1m_z(wk3ia, wk3ib, DECOMP_2D_FFT_FORWARD)
            else
               if (.not. present(out_real)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
               call r2r_1m_z(wk3ra, out_real, DECOMP_2D_FFT_FORWARD)
            end if
         end if

      else

         ! DFT / DTT in z
         if (dtt_z_dft) then
            call r2rr_1m_z(in, wk3ra, wk3ia)
            cplx = .true.
         else
            call r2r_1m_z(in, wk3ra, DECOMP_2D_FFT_FORWARD)
         end if

         ! Transpose z => y
         call transpose_z_to_y(wk3ra, wk2rb, dtt_decomp_yz)
         if (cplx) then
            call transpose_z_to_y(wk3ia, wk2ib, dtt_decomp_yz)
         end if

         ! DFT / DTT in y
         if (dtt_y_dft) then
            if (cplx) then
               call rr2rr_1m_y(wk2rb, wk2ib, wk2ra, wk2ia, DECOMP_2D_FFT_FORWARD)
            else
               call r2rr_1m_y(wk2rb, wk2ra, wk2ia)
               cplx = .true.
            end if
         else
            call r2r_1m_y(wk2rb, wk2ra, DECOMP_2D_FFT_FORWARD)
            if (cplx) then
               call r2r_1m_y(wk2ib, wk2ia, DECOMP_2D_FFT_FORWARD)
            end if
         end if

         ! Transpose y => x
         call transpose_y_to_x(wk2ra, wk1rb, dtt_decomp_xy)
         if (cplx) then
            call transpose_y_to_x(wk2ia, wk1ib, dtt_decomp_xy)
         end if

         ! DFT / DTT in x
         if (dtt_x_dft) then
            if (cplx) then
               call rr2rr_1m_x(wk1rb, wk1ib, wk1ra, wk1ia, DECOMP_2D_FFT_FORWARD)
            else
               call r2rr_1m_x(wk1rb, wk1ra, wk1ia)
               cplx = .true.
            end if
         else
            if (cplx) then
               call r2r_1m_x(wk1rb, wk1ra, DECOMP_2D_FFT_FORWARD)
               call r2r_1m_x(wk1ib, wk1ia, DECOMP_2D_FFT_FORWARD)
            else
               if (.not. present(out_real)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
               call r2r_1m_x(wk1rb, out_real, DECOMP_2D_FFT_FORWARD)
            end if
         end if

      end if

      ! Safety check
      if (cplx .and. (.not. present(out_cplx))) then
         call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
      end if

      ! Recombine the real and imag parts if needed
      if (cplx .and. format == PHYSICAL_IN_X) then
         do k = 1, dtt_decomp_sp%zsz(3)
            do j = 1, dtt_decomp_sp%zsz(2)
               do i = 1, dtt_decomp_sp%zsz(1)
                  out_cplx(i, j, k) = cmplx(wk3rb(i, j, k), wk3ib(i, j, k), kind=mytype)
               end do
            end do
         end do
      else if (cplx) then
         do k = 1, dtt_decomp_sp%xsz(3)
            do j = 1, dtt_decomp_sp%xsz(2)
               do i = 1, dtt_decomp_sp%xsz(1)
                  out_cplx(i, j, k) = cmplx(wk1ra(i, j, k), wk1ia(i, j, k), kind=mytype)
               end do
            end do
         end do
      end if

      if (decomp_profiler_fft) call decomp_profiler_end("decomp_2d_dtt_3d_r2x")

   end subroutine decomp_2d_dtt_3d_r2x

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Backward 3D DTT - real/complex to real
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_2d_dtt_3d_x2r(in_real, in_cplx, out)

      implicit none

      ! Arguments
      real(mytype), dimension(:, :, :), contiguous, target, intent(inout), optional :: in_real
      complex(mytype), dimension(:, :, :), contiguous, target, intent(inout), optional :: in_cplx
      real(mytype), dimension(:, :, :), contiguous, target, intent(out), optional :: out

      ! Local variables
      logical :: cplx
      integer :: i, j, k

      if (decomp_profiler_fft) call decomp_profiler_start("decomp_2d_dtt_3d_x2r")

      ! Init
      if (dtt_x_dft .or. dtt_y_dft .or. dtt_z_dft) then
         cplx = .true.
      else
         cplx = .false.
      end if

      ! Safety check
      if (cplx) then
         if (.not. present(in_cplx)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
      else
         if (.not. present(in_real)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
      end if
      if (.not. present(out)) call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")

      ! Split real and imag parts if needed
      if (cplx) then
         if (format == PHYSICAL_IN_X) then
            do k = 1, dtt_decomp_sp%zsz(3)
               do j = 1, dtt_decomp_sp%zsz(2)
                  do i = 1, dtt_decomp_sp%zsz(1)
                     wk3rb(i, j, k) = real(in_cplx(i, j, k), kind=mytype)
                     wk3ib(i, j, k) = aimag(in_cplx(i, j, k))
                  end do
               end do
            end do
         else
            do k = 1, dtt_decomp_sp%xsz(3)
               do j = 1, dtt_decomp_sp%xsz(2)
                  do i = 1, dtt_decomp_sp%xsz(1)
                     wk1ra(i, j, k) = real(in_cplx(i, j, k), kind=mytype)
                     wk1ia(i, j, k) = aimag(in_cplx(i, j, k))
                  end do
               end do
            end do
         end if
      end if

      ! Perform the 3D DTT
      if (format == PHYSICAL_IN_Z) then

         ! DFT / DTT in x
         if (dtt_x_dft) then
            if (dtt_y_dft .or. dtt_z_dft) then
               call rr2rr_1m_x(wk1ra, wk1ia, wk1rb, wk1ib, DECOMP_2D_FFT_BACKWARD)
            else
               call rr2r_1m_x(wk1ra, wk1ia, wk1rb)
               cplx = .false.
            end if
         else
            if (cplx) then
               call r2r_1m_x(wk1ra, wk1rb, DECOMP_2D_FFT_BACKWARD)
               call r2r_1m_x(wk1ia, wk1ib, DECOMP_2D_FFT_BACKWARD)
            else
               call r2r_1m_x(in_real, wk1rb, DECOMP_2D_FFT_BACKWARD)
            end if
         end if

         ! Transpose x => y
         call transpose_x_to_y(wk1rb, wk2ra, dtt_decomp_xy)
         if (cplx) then
            call transpose_x_to_y(wk1ib, wk2ia, dtt_decomp_xy)
         end if

         ! DFT / DTT in y
         if (dtt_y_dft) then
            if (dtt_z_dft) then
               call rr2rr_1m_y(wk2ra, wk2ia, wk2rb, wk2ib, DECOMP_2D_FFT_BACKWARD)
            else
               call rr2r_1m_y(wk2ra, wk2ia, wk2rb)
               cplx = .false.
            end if
         else
            call r2r_1m_y(wk2ra, wk2rb, DECOMP_2D_FFT_BACKWARD)
            if (cplx) call r2r_1m_y(wk2ia, wk2ib, DECOMP_2D_FFT_BACKWARD)
         end if

         ! Transpose y => z
         call transpose_y_to_z(wk2rb, wk3ra, dtt_decomp_yz)
         if (cplx) then
            call transpose_y_to_z(wk2ib, wk3ia, dtt_decomp_yz)
         end if

         ! DFT / DTT in z
         if (cplx) then
            call rr2r_1m_z(wk3ra, wk3ia, out)
         else
            call r2r_1m_z(wk3ra, out, DECOMP_2D_FFT_BACKWARD)
         end if

      else

         ! DFT / DTT in z
         if (dtt_z_dft) then
            if (dtt_y_dft .or. dtt_x_dft) then
               call rr2rr_1m_z(wk3rb, wk3ib, wk3ra, wk3ia, DECOMP_2D_FFT_BACKWARD)
            else
               call rr2r_1m_z(wk3rb, wk3ib, wk3ra)
               cplx = .false.
            end if
         else
            if (cplx) then
               call r2r_1m_z(wk3rb, wk3ra, DECOMP_2D_FFT_BACKWARD)
               call r2r_1m_z(wk3ib, wk3ia, DECOMP_2D_FFT_BACKWARD)
            else
               call r2r_1m_z(in_real, wk3ra, DECOMP_2D_FFT_BACKWARD)
            end if
         end if

         ! Transpose z => y
         call transpose_z_to_y(wk3ra, wk2rb, dtt_decomp_yz)
         if (cplx) then
            call transpose_z_to_y(wk3ia, wk2ib, dtt_decomp_yz)
         end if

         ! DFT / DTT in y
         if (dtt_y_dft) then
            if (dtt_x_dft) then
               call rr2rr_1m_y(wk2rb, wk2ib, wk2ra, wk2ia, DECOMP_2D_FFT_BACKWARD)
            else
               call rr2r_1m_y(wk2rb, wk2ib, wk2ra)
               cplx = .false.
            end if
         else
            call r2r_1m_y(wk2rb, wk2ra, DECOMP_2D_FFT_BACKWARD)
            if (cplx) call r2r_1m_y(wk2ib, wk2ia, DECOMP_2D_FFT_BACKWARD)
         end if

         ! Transpose y => x
         call transpose_y_to_x(wk2ra, wk1rb, dtt_decomp_xy)
         if (cplx) then
            call transpose_y_to_x(wk2ia, wk1ib, dtt_decomp_xy)
         end if

         ! DFT / DTT in x
         if (cplx) then
            call rr2r_1m_x(wk1rb, wk1ib, out)
         else
            call r2r_1m_x(wk1rb, out, DECOMP_2D_FFT_BACKWARD)
         end if

      end if

      if (decomp_profiler_fft) call decomp_profiler_end("decomp_2d_dtt_3d_x2r")

   end subroutine decomp_2d_dtt_3d_x2r

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
   subroutine wrapper_rr2rr(plan, &
                            inr, &
                            ini, &
                            ii, ij, ik, isz, &
                            outr, &
                            outi, &
                            oi, oj, ok, osz)

      implicit none

      ! Arguments
      type(c_ptr), intent(in) :: plan
      real(mytype), dimension(:, :, :), target, contiguous, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), target, contiguous, intent(out) :: outr, outi
      integer, intent(in) :: ii, ij, ik, isz, oi, oj, ok, osz

      ! Local variables
      real(mytype), dimension(:), contiguous, pointer :: inr2, ini2, outr2, outi2

      ! Create 1D pointers starting at the ifirst / ofirst location
      call c_f_pointer(c_loc(inr(ii, ij, ik)), inr2, (/isz/))
      call c_f_pointer(c_loc(ini(ii, ij, ik)), ini2, (/isz/))
      call c_f_pointer(c_loc(outr(oi, oj, ok)), outr2, (/osz/))
      call c_f_pointer(c_loc(outi(oi, oj, ok)), outi2, (/osz/))

      ! Perform DFT
#ifdef DOUBLE_PREC
      call fftw_execute_split_dft(plan, inr2, ini2, outr2, outi2)
#else
      call fftwf_execute_split_dft(plan, inr2, ini2, outr2, outi2)
#endif

      ! Release pointers
      nullify (inr2)
      nullify (ini2)
      nullify (outr2)
      nullify (outi2)

   end subroutine wrapper_rr2rr
   !
   subroutine wrapper_r2r(plan, &
                          inr, &
                          ii, ij, ik, isz, &
                          outr, &
                          oi, oj, ok, osz)

      implicit none

      ! Arguments
      type(c_ptr), intent(in) :: plan
      real(mytype), dimension(:, :, :), target, contiguous, intent(inout) :: inr
      real(mytype), dimension(:, :, :), target, contiguous, intent(out) :: outr
      integer, intent(in) :: ii, ij, ik, isz, oi, oj, ok, osz

      ! Local variables
      real(mytype), dimension(:), contiguous, pointer :: inr2, outr2

      ! Create 1D pointers starting at the ifirst / ofirst location
      call c_f_pointer(c_loc(inr(ii, ij, ik)), inr2, (/isz/))
      call c_f_pointer(c_loc(outr(oi, oj, ok)), outr2, (/osz/))

      ! Perform DFT
#ifdef DOUBLE_PREC
      call fftw_execute_r2r(plan, inr2, outr2)
#else
      call fftwf_execute_r2r(plan, inr2, outr2)
#endif

      ! Release pointers
      nullify (inr2)
      nullify (outr2)

   end subroutine wrapper_r2r
   !
   subroutine wrapper_rr2r(plan, &
                           inr, &
                           ini, &
                           ii, ij, ik, isz, &
                           outr, &
                           oi, oj, ok, osz)

      implicit none

      ! Arguments
      type(c_ptr), intent(in) :: plan
      real(mytype), dimension(:, :, :), target, contiguous, intent(inout) :: inr, ini
      real(mytype), dimension(:, :, :), target, contiguous, intent(out) :: outr
      integer, intent(in) :: ii, ij, ik, isz, oi, oj, ok, osz

      ! Local variables
      real(mytype), dimension(:), contiguous, pointer :: inr2, ini2, outr2

      ! Create 1D pointers starting at the ifirst / ofirst location
      call c_f_pointer(c_loc(inr(ii, ij, ik)), inr2, (/isz/))
      call c_f_pointer(c_loc(ini(ii, ij, ik)), ini2, (/isz/))
      call c_f_pointer(c_loc(outr(oi, oj, ok)), outr2, (/osz/))

      ! Perform DFT
#ifdef DOUBLE_PREC
      call fftw_execute_split_dft_c2r(plan, inr2, ini2, outr2)
#else
      call fftwf_execute_split_dft_c2r(plan, inr2, ini2, outr2)
#endif

      ! Release pointers
      nullify (inr2)
      nullify (ini2)
      nullify (outr2)

   end subroutine wrapper_rr2r
   !
   subroutine wrapper_r2rr(plan, &
                           inr, &
                           ii, ij, ik, isz, &
                           outr, &
                           outi, &
                           oi, oj, ok, osz)

      implicit none

      ! Arguments
      type(c_ptr), intent(in) :: plan
      real(mytype), dimension(:, :, :), target, contiguous, intent(inout) :: inr
      real(mytype), dimension(:, :, :), target, contiguous, intent(out) :: outr, outi
      integer, intent(in) :: ii, ij, ik, isz, oi, oj, ok, osz

      ! Local variables
      real(mytype), dimension(:), contiguous, pointer :: inr2, outr2, outi2

      ! Create 1D pointers starting at the ifirst / ofirst location
      call c_f_pointer(c_loc(inr(ii, ij, ik)), inr2, (/isz/))
      call c_f_pointer(c_loc(outr(oi, oj, ok)), outr2, (/osz/))
      call c_f_pointer(c_loc(outi(oi, oj, ok)), outi2, (/osz/))

      ! Perform DFT
#ifdef DOUBLE_PREC
      call fftw_execute_split_dft_r2c(plan, inr2, outr2, outi2)
#else
      call fftwf_execute_split_dft_r2c(plan, inr2, outr2, outi2)
#endif

      ! Release pointers
      nullify (inr2)
      nullify (outr2)
      nullify (outi2)

   end subroutine wrapper_r2rr

end module decomp_2d_fft
