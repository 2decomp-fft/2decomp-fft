!! SPDX-License-Identifier: BSD-3-Clause

! This is the FFTW implementation of the FFT library using
! the Fortran 2003 interface introduced in FFTW 3.3-beta1

module decomp_2d_fft

   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d  ! 2D decomposition module
   use, intrinsic :: iso_c_binding

   implicit none

   include "fftw3.f03"

   private        ! Make everything private unless declared public

   ! engine-specific global variables
   integer, save :: plan_type = FFTW_MEASURE

   ! FFTW plans
   ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
   ! For c2c transforms:
   !     use plan(-1,j) for  forward transform;
   !     use plan( 1,j) for backward transform;
   ! For r2c/c2r transforms:
   !     use plan(0,j) for r2c transforms;
   !     use plan(2,j) for c2r transforms;
   type(C_PTR), pointer, save :: plan(:, :)

   integer, parameter, public :: D2D_FFT_BACKEND = D2D_FFT_BACKEND_FFTW3_F03

   integer, pointer, save :: format                 ! input X-pencil or Z-pencil

   ! Global size of the FFT
   integer, pointer, save :: nx_fft, ny_fft, nz_fft

   ! 2D processor grid
   ! FIXME this is already available in the module decomp_2d
   integer, save, dimension(2) :: dims

   ! Decomposition objects
   TYPE(DECOMP_INFO), pointer, save :: ph => null()  ! physical space
   TYPE(DECOMP_INFO), pointer, save :: sp  ! spectral space

   ! Workspace to store the intermediate Y-pencil data
   ! *** TODO: investigate how to use only one workspace array
   complex(mytype), contiguous, pointer :: wk2_c2c(:, :, :), wk2_r2c(:, :, :), wk13(:, :, :)

   ! Derived type with all the quantities needed to perform FFT
   type decomp_2d_fft_engine
      type(c_ptr), private :: plan(-1:2, 3), wk2_c2c_p, wk13_p
      integer, private :: format
      logical, private :: initialised = .false.
      integer, private :: nx_fft, ny_fft, nz_fft
      type(decomp_info), pointer, public :: ph => null()
      type(decomp_info), public :: sp
      complex(mytype), contiguous, pointer, private :: wk2_c2c(:, :, :) => null()
      complex(mytype), contiguous, pointer, private :: wk2_r2c(:, :, :) => null()
      complex(mytype), contiguous, pointer, private :: wk13(:, :, :) => null()
   contains
      procedure, public :: init => decomp_2d_fft_engine_init
      procedure, public :: fin => decomp_2d_fft_engine_fin
      procedure, public :: use_it => decomp_2d_fft_engine_use_it
      generic, public :: fft => c2c, r2c, c2r
      procedure, private :: c2c => decomp_2d_fft_engine_fft_c2c
      procedure, private :: r2c => decomp_2d_fft_engine_fft_r2c
      procedure, private :: c2r => decomp_2d_fft_engine_fft_c2r
   end type decomp_2d_fft_engine

   !
   ! Multigrid options
   !
   ! Number of FFT grids
   integer, save :: n_grid = 0
   ! Array to store all the grids
   type(decomp_2d_fft_engine), allocatable, target, save :: fft_engines(:)

   public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
             decomp_2d_fft_finalize, decomp_2d_fft_get_size, &
             decomp_2d_fft_get_ph, decomp_2d_fft_get_sp, &
             decomp_2d_fft_get_ngrid, decomp_2d_fft_set_ngrid, &
             decomp_2d_fft_use_grid, decomp_2d_fft_engine

   ! Declare generic interfaces to handle different inputs

   interface decomp_2d_fft_init
      module procedure fft_init_noarg
      module procedure fft_init_arg
      module procedure fft_init_general
      module procedure fft_init_multigrid
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

   subroutine fft_init_arg(pencil)     ! allow to handle Z-pencil input

      implicit none

      integer, intent(IN) :: pencil

      call fft_init_general(pencil, nx_global, ny_global, nz_global)

   end subroutine fft_init_arg

   ! Initialise one FFT library to perform arbitrary size transforms
   subroutine fft_init_general(pencil, nx, ny, nz)

      implicit none

      integer, intent(IN) :: pencil, nx, ny, nz

      ! Only one FFT engine will be used
      call decomp_2d_fft_set_ngrid(1)

      ! Initialise the FFT engine
      call decomp_2d_fft_init(pencil, nx, ny, nz, 1)

   end subroutine fft_init_general

   ! Initialise the provided FFT grid
   subroutine fft_init_multigrid(pencil, nx, ny, nz, igrid)

      implicit none

      integer, intent(in) :: pencil, nx, ny, nz, igrid

      ! Safety check
      if (igrid < 1 .or. igrid > n_grid) then
         call decomp_2d_abort(__FILE__, __LINE__, igrid, "Invalid value for igrid")
      end if

      ! Initialise the engine
      call fft_engines(igrid)%init(pencil, nx, ny, nz)

   end subroutine fft_init_multigrid

   ! Initialise the given FFT engine
   subroutine decomp_2d_fft_engine_init(engine, pencil, nx, ny, nz)

      implicit none

      class(decomp_2d_fft_engine), intent(inout), target :: engine
      integer, intent(in) :: pencil, nx, ny, nz

      integer(C_SIZE_T) :: sz

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_start("fft_init")
#endif

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
         if (.not. associated(engine%ph)) allocate (engine%ph)
         call decomp_info_init(nx, ny, nz, engine%ph)
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
      if (pencil == PHYSICAL_IN_X) then
         sz = product(engine%sp%xsz)
         engine%wk13_p = fftw_alloc_complex(sz)
         call c_f_pointer(engine%wk13_p, engine%wk13, engine%sp%xsz)
      else if (pencil == PHYSICAL_IN_Z) then
         sz = product(engine%sp%zsz)
         engine%wk13_p = fftw_alloc_complex(sz)
         call c_f_pointer(engine%wk13_p, engine%wk13, engine%sp%zsz)
      end if

      ! Warning : replace the default engine
      call engine%use_it()

      ! Compute the plan for the current engine
      call init_fft_engine()

      ! Tag the engine as initialised
      engine%initialised = .true.

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_end("fft_init")
#endif

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

      ! Clean the FFTW library
      call fftw_cleanup()

   end subroutine decomp_2d_fft_finalize

   ! Clean the given FFT engine
   subroutine decomp_2d_fft_engine_fin(engine)

      implicit none

      class(decomp_2d_fft_engine), intent(inout) :: engine

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_start("fft_fin")
#endif

      if (engine%nx_fft /= nx_global .or. &
          engine%ny_fft /= ny_global .or. &
          engine%nz_fft /= nz_global) then
         call decomp_info_finalize(engine%ph)
         deallocate (engine%ph)
      end if
      nullify (engine%ph)
      call decomp_info_finalize(engine%sp)

      call fftw_free(engine%wk2_c2c_p)
      nullify (engine%wk2_c2c)
      nullify (engine%wk2_r2c)
      call fftw_free(engine%wk13_p)
      nullify (engine%wk13)

      call finalize_fft_engine(engine%plan)

      engine%initialised = .false.

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_end("fft_fin")
#endif

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

      ! Safety checks
      if (ngrd < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, ngrd, "Invalid value for n_grid")
      end if
      if (n_grid > 0 .or. allocated(fft_engines)) then
         call decomp_2d_abort(__FILE__, __LINE__, n_grid, &
                              "The number of FFT grids was already set")
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
   subroutine decomp_2d_fft_engine_use_it(engine)

      implicit none

      class(decomp_2d_fft_engine), intent(in), target :: engine

      ! Safety checks
      if (.not. engine%initialised) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, "FFT engine is not ready")
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

   end subroutine decomp_2d_fft_engine_use_it

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in X direction
   subroutine c2c_1m_x_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

#ifdef DOUBLE_PREC
      complex(C_DOUBLE_COMPLEX), pointer :: a1(:, :, :)
      complex(C_DOUBLE_COMPLEX), pointer :: a1o(:, :, :)
#else
      complex(C_FLOAT_COMPLEX), pointer :: a1(:, :, :)
      complex(C_FLOAT_COMPLEX), pointer :: a1o(:, :, :)
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
      complex(C_DOUBLE_COMPLEX), pointer :: a1(:, :)
      complex(C_DOUBLE_COMPLEX), pointer :: a1o(:, :)
#else
      complex(C_FLOAT_COMPLEX), pointer :: a1(:, :)
      complex(C_FLOAT_COMPLEX), pointer :: a1o(:, :)
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
      complex(C_DOUBLE_COMPLEX), pointer :: a1(:, :, :)
      complex(C_DOUBLE_COMPLEX), pointer :: a1o(:, :, :)
#else
      complex(C_FLOAT_COMPLEX), pointer :: a1(:, :, :)
      complex(C_FLOAT_COMPLEX), pointer :: a1o(:, :, :)
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

      real(mytype), pointer :: a1(:, :, :)
      complex(mytype), pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      sz = decomp_ph%xsz(1) * decomp_ph%xsz(2) * decomp_ph%xsz(3)
      a1_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, &
                       [decomp_ph%xsz(1), decomp_ph%xsz(2), decomp_ph%xsz(3)])
      sz = decomp_sp%xsz(1) * decomp_sp%xsz(2) * decomp_sp%xsz(3)
      a2_p = fftw_alloc_complex(sz)
      call c_f_pointer(a2_p, a2, &
                       [decomp_sp%xsz(1), decomp_sp%xsz(2), decomp_sp%xsz(3)])

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
      call fftw_free(a2_p)

      return
   end subroutine r2c_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
   subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

      complex(mytype), pointer :: a1(:, :, :)
      real(mytype), pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      sz = decomp_sp%xsz(1) * decomp_sp%xsz(2) * decomp_sp%xsz(3)
      a1_p = fftw_alloc_complex(sz)
      call c_f_pointer(a1_p, a1, &
                       [decomp_sp%xsz(1), decomp_sp%xsz(2), decomp_sp%xsz(3)])
      sz = decomp_ph%xsz(1) * decomp_ph%xsz(2) * decomp_ph%xsz(3)
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a2_p, a2, &
                       [decomp_ph%xsz(1), decomp_ph%xsz(2), decomp_ph%xsz(3)])

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

      call fftw_free(a1_p)
      call fftw_free(a2_p)

      return
   end subroutine c2r_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
   subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

      real(mytype), pointer :: a1(:, :, :)
      complex(mytype), pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      sz = decomp_ph%zsz(1) * decomp_ph%zsz(2) * decomp_ph%zsz(3)
      a1_p = fftw_alloc_real(sz)
      call c_f_pointer(a1_p, a1, &
                       [decomp_ph%zsz(1), decomp_ph%zsz(2), decomp_ph%zsz(3)])
      sz = decomp_sp%zsz(1) * decomp_sp%zsz(2) * decomp_sp%zsz(3)
      a2_p = fftw_alloc_complex(sz)
      call c_f_pointer(a2_p, a2, &
                       [decomp_sp%zsz(1), decomp_sp%zsz(2), decomp_sp%zsz(3)])

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
      call fftw_free(a2_p)

      return
   end subroutine r2c_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
   subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

      implicit none

      type(C_PTR) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

      complex(mytype), pointer :: a1(:, :, :)
      real(mytype), pointer :: a2(:, :, :)
      type(C_PTR) :: a1_p, a2_p
      integer(C_SIZE_T) :: sz

      sz = decomp_sp%zsz(1) * decomp_sp%zsz(2) * decomp_sp%zsz(3)
      a1_p = fftw_alloc_complex(sz)
      call c_f_pointer(a1_p, a1, &
                       [decomp_sp%zsz(1), decomp_sp%zsz(2), decomp_sp%zsz(3)])
      sz = decomp_ph%zsz(1) * decomp_ph%zsz(2) * decomp_ph%zsz(3)
      a2_p = fftw_alloc_real(sz)
      call c_f_pointer(a2_p, a2, &
                       [decomp_ph%zsz(1), decomp_ph%zsz(2), decomp_ph%zsz(3)])

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

      call fftw_free(a1_p)
      call fftw_free(a2_p)

      return
   end subroutine c2r_1m_z_plan

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

      return
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

#ifdef DOUBLE_PREC
      call fftw_execute_dft(plan1, inout, inout)
#else
      call fftwf_execute_dft(plan1, inout, inout)
#endif

      return
   end subroutine c2c_1m_x

   ! c2c transform, multiple 1D FFTs in y direction
   subroutine c2c_1m_y(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR) :: plan1

      integer :: k, s3

      s3 = size(inout, 3)

      do k = 1, s3  ! transform on one Z-plane at a time
#ifdef DOUBLE_PREC
         call fftw_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
#else
         call fftwf_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
#endif
      end do

      return
   end subroutine c2c_1m_y

   ! c2c transform, multiple 1D FFTs in z direction
   subroutine c2c_1m_z(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR) :: plan1

#ifdef DOUBLE_PREC
      call fftw_execute_dft(plan1, inout, inout)
#else
      call fftwf_execute_dft(plan1, inout, inout)
#endif

      return
   end subroutine c2c_1m_z

   ! r2c transform, multiple 1D FFTs in x direction
   subroutine r2c_1m_x(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      complex(mytype), dimension(:, :, :), intent(OUT) :: output

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
      complex(mytype), dimension(:, :, :), intent(OUT) :: output

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
      real(mytype), dimension(:, :, :), intent(OUT) :: output

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
      real(mytype), dimension(:, :, :), intent(OUT) :: output

#ifdef DOUBLE_PREC
      call fftw_execute_dft_c2r(plan(2, 3), input, output)
#else
      call fftwf_execute_dft_c2r(plan(2, 3), input, output)
#endif

      return

   end subroutine c2r_1m_z

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2c(in, out, isign)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out
      integer, intent(IN) :: isign

#ifndef OVERWRITE
      complex(mytype), pointer :: wk1(:, :, :)
      integer(C_SIZE_T) :: sz
      type(C_PTR) :: wk1_p

      wk1_p = c_null_ptr ! Initialise to NULL pointer
#endif

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2c")
#endif

      if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
          format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

         ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
         call c2c_1m_x(in, plan(isign, 1))
#else
         sz = ph%xsz(1) * ph%xsz(2) * ph%xsz(3)
         wk1_p = fftw_alloc_complex(sz)
         call c_f_pointer(wk1_p, wk1, [ph%xsz(1), ph%xsz(2), ph%xsz(3)])
         wk1 = in
         call c2c_1m_x(wk1, plan(isign, 1))
#endif

         ! ===== Swap X --> Y; 1D FFTs in Y =====

         if (dims(1) > 1) then
#ifdef OVERWRITE
            call transpose_x_to_y(in, wk2_c2c, ph)
#else
            call transpose_x_to_y(wk1, wk2_c2c, ph)
#endif
            call c2c_1m_y(wk2_c2c, plan(isign, 2))
         else
#ifdef OVERWRITE
            call c2c_1m_y(in, plan(isign, 2))
#else
            call c2c_1m_y(wk1, plan(isign, 2))
#endif
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_c2c, out, ph)
         else
#ifdef OVERWRITE
            call transpose_y_to_z(in, out, ph)
#else
            call transpose_y_to_z(wk1, out, ph)
#endif
         end if
         call c2c_1m_z(out, plan(isign, 3))

      else if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_BACKWARD &
               .OR. &
               format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_FORWARD) then

         ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
         call c2c_1m_z(in, plan(isign, 3))
#else
         sz = ph%zsz(1) * ph%zsz(2) * ph%zsz(3)
         wk1_p = fftw_alloc_complex(sz)
         call c_f_pointer(wk1_p, wk1, [ph%zsz(1), ph%zsz(2), ph%zsz(3)])
         wk1 = in
         call c2c_1m_z(wk1, plan(isign, 3))
#endif

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
#ifdef OVERWRITE
            call transpose_z_to_y(in, wk2_c2c, ph)
#else
            call transpose_z_to_y(wk1, wk2_c2c, ph)
#endif
            call c2c_1m_y(wk2_c2c, plan(isign, 2))
         else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
            call transpose_z_to_y(in, out, ph)
#else
            call transpose_z_to_y(wk1, out, ph)
#endif
            call c2c_1m_y(out, plan(isign, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_c2c, out, ph)
         end if
         call c2c_1m_x(out, plan(isign, 1))

      end if

#ifndef OVERWRITE
      call fftw_free(wk1_p)
      nullify (wk1)
#endif

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2c")
#endif

      return
   end subroutine fft_3d_c2c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_r2c(in_r, out_c)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT) :: in_r
      complex(mytype), dimension(:, :, :), intent(OUT) :: out_c

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")
#endif

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in X =====
         call r2c_1m_x(in_r, wk13)

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            call transpose_x_to_y(wk13, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, plan(0, 2))
         else
            call c2c_1m_y(wk13, plan(0, 2))
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, out_c, sp)
         else
            call transpose_y_to_z(wk13, out_c, sp)
         end if
         call c2c_1m_z(out_c, plan(0, 3))

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in Z =====
         call r2c_1m_z(in_r, wk13)

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            call transpose_z_to_y(wk13, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, plan(0, 2))
         else  ! out_c==wk2_r2c if 1D decomposition
            call transpose_z_to_y(wk13, out_c, sp)
            call c2c_1m_y(out_c, plan(0, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_r2c, out_c, sp)
         end if
         call c2c_1m_x(out_c, plan(0, 1))

      end if

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_end("fft_r2c")
#endif

      return
   end subroutine fft_3d_r2c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D inverse FFT - complex to real
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2r(in_c, out_r)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: in_c
      real(mytype), dimension(:, :, :), intent(OUT) :: out_r

#ifndef OVERWRITE
      complex(mytype), pointer :: wk1(:, :, :)
      integer(C_SIZE_T) :: sz
      type(C_PTR) :: wk1_p

      wk1_p = c_null_ptr ! Initialise to NULL pointer
#endif

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2r")
#endif

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
         call c2c_1m_z(in_c, plan(2, 3))
#else
         sz = sp%zsz(1) * sp%zsz(2) * sp%zsz(3)
         wk1_p = fftw_alloc_complex(sz)
         call c_f_pointer(wk1_p, wk1, [sp%zsz(1), sp%zsz(2), sp%zsz(3)])
         wk1 = in_c
         call c2c_1m_z(wk1, plan(2, 3))
#endif

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
         call transpose_z_to_y(in_c, wk2_r2c, sp)
#else
         call transpose_z_to_y(wk1, wk2_r2c, sp)
#endif
         call c2c_1m_y(wk2_r2c, plan(2, 2))

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_r2c, wk13, sp)
            call c2r_1m_x(wk13, out_r)
         else
            call c2r_1m_x(wk2_r2c, out_r)
         end if

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
         call c2c_1m_x(in_c, plan(2, 1))
#else
         sz = sp%xsz(1) * sp%xsz(2) * sp%xsz(3)
         wk1_p = fftw_alloc_complex(sz)
         call c_f_pointer(wk1_p, wk1, [sp%xsz(1), sp%xsz(2), sp%xsz(3)])
         wk1 = in_c
         call c2c_1m_x(wk1, plan(2, 1))
#endif

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
#ifdef OVERWRITE
            call transpose_x_to_y(in_c, wk2_r2c, sp)
#else
            call transpose_x_to_y(wk1, wk2_r2c, sp)
#endif
            call c2c_1m_y(wk2_r2c, plan(2, 2))
         else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
            call c2c_1m_y(in_c, plan(2, 2))
#else
            call c2c_1m_y(wk1, plan(2, 2))
#endif
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, wk13, sp)
         else
#ifdef OVERWRITE
            call transpose_y_to_z(in_c, wk13, sp)
#else
            call transpose_y_to_z(wk1, wk13, sp)
#endif
         end if
         call c2r_1m_z(wk13, out_r)

      end if

#ifndef OVERWRITE
      call fftw_free(wk1_p)
      nullify (wk1)
#endif

#ifdef PROFILER
      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2r")
#endif

      return
   end subroutine fft_3d_c2r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Wrappers for calling 3D FFT directly using the engine object
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_2d_fft_engine_fft_c2c(engine, in, out, isign)

      implicit none

      class(decomp_2d_fft_engine), intent(in) :: engine
      complex(mytype), dimension(:, :, :), intent(INOUT) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out
      integer, intent(IN) :: isign

      call engine%use_it()
      call decomp_2d_fft_3d(in, out, isign)

   end subroutine decomp_2d_fft_engine_fft_c2c

   subroutine decomp_2d_fft_engine_fft_r2c(engine, in, out)

      implicit none

      class(decomp_2d_fft_engine), intent(in) :: engine
      real(mytype), dimension(:, :, :), intent(INOUT) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out

      call engine%use_it()
      call decomp_2d_fft_3d(in, out)

   end subroutine decomp_2d_fft_engine_fft_r2c

   subroutine decomp_2d_fft_engine_fft_c2r(engine, in, out)

      implicit none

      class(decomp_2d_fft_engine), intent(in) :: engine
      complex(mytype), dimension(:, :, :), intent(INOUT) :: in
      real(mytype), dimension(:, :, :), intent(OUT) :: out

      call engine%use_it()
      call decomp_2d_fft_3d(in, out)

   end subroutine decomp_2d_fft_engine_fft_c2r

end module decomp_2d_fft
