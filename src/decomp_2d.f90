!! SPDX-License-Identifier: BSD-3-Clause

! This is the main 2D pencil decomposition module

module decomp_2d

   use MPI
   use iso_c_binding, only: c_size_t
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use factor
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_mpi, only: DECOMP_2D_COMM_CART_X => DECOMP_2D_COMM_CART_X
   use decomp_2d_mpi, only: DECOMP_2D_COMM_CART_Y => DECOMP_2D_COMM_CART_Y
   use decomp_2d_mpi, only: DECOMP_2D_COMM_CART_Z => DECOMP_2D_COMM_CART_Z
   use decomp_2d_profiler
#if defined(_GPU)
   use cudafor
   use decomp_2d_cumpi
#if defined(_NCCL)
   use nccl
   use decomp_2d_nccl
#endif
#endif
   use m_info, only: decomp_info => decomp_info ! Expose decomp_info from 2decomp
   use m_info, only: decomp_main => decomp_main ! Expose decomp_main from 2decomp
   use m_decomp_pool
   use m_mem_pool
   use m_halo, only: init_neighbour

   implicit none

   private        ! Make everything private unless declared public

   ! Default parameter opt_global in alloc subroutines
   logical, parameter :: default_opt_global = .false.

   ! some key global variables
   integer, save, public :: nx_global, ny_global, nz_global  ! global size
   logical, save, public, protected :: use_pool

   !
   ! Output for the log can be changed by the external code before calling decomp_2d_init
   !
   !    0 => No log output
   !    1 => Master rank log output to stdout
   !    2 => Master rank log output to the file "decomp_2d_setup.log"
   !    3 => All ranks log output to a dedicated file
   !
   ! The default value is 2 (3 for debug builds)
   !
#ifdef DEBUG
   integer, public, save :: decomp_log = D2D_LOG_TOFILE_FULL
#else
   integer, public, save :: decomp_log = D2D_LOG_TOFILE
#endif

   !
   ! Debug level can be changed by the external code before calling decomp_2d_init
   !
   ! The environment variable "DECOMP_2D_DEBUG" can be used to change the debug level
   !
   ! Debug checks are performed only when the preprocessor variable DEBUG is defined
   !
#ifdef DEBUG
   integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_INFO
#else
   integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_OFF
#endif

   ! staring/ending index and size of data held by current processor
   ! duplicate 'decomp_main', needed by apps to define data structure
   integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
   integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
   integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

   ! These are the buffers used by MPI_ALLTOALL(V) calls
   integer(c_size_t), save :: decomp_buf_size = 0
#if defined(_GPU)
   ! Shared real/complex GPU buffers
   real(mytype), target, device, allocatable, dimension(:) :: work1, work2
   ! Real / complex pointers to GPU buffers
   real(mytype), pointer, contiguous, dimension(:) :: work1_r, work2_r
   complex(mytype), pointer, contiguous, dimension(:) :: work1_c, work2_c
#endif

   ! public user routines
   public :: decomp_2d_init, decomp_2d_finalize, &
             transpose_x_to_y, transpose_y_to_z, &
             transpose_z_to_y, transpose_y_to_x, &
             decomp_info, decomp_main, &
             decomp_info_init, decomp_info_finalize, partition, &
             decomp_info_print, &
             alloc_x, alloc_y, alloc_z, &
             get_decomp_info, &
             get_decomp_dims, &
             d2d_log_is_active, &
             d2d_log_get_unit, &
             d2d_log_close_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! These are routines to perform global data transpositions
   !
   !   Four combinations are available, enough to cover all situations
   !    - transpose_x_to_y (X-pencil --> Y-pencil)
   !    - transpose_y_to_z (Y-pencil --> Z-pencil)
   !    - transpose_z_to_y (Z-pencil --> Y-pencil)
   !    - transpose_y_to_x (Y-pencil --> X-pencil)
   !
   !   Generic interface provided here to support multiple data types
   !    - real and complex types supported through generic interface
   !    - single/double precision supported through pre-processing
   !       * see 'mytype' variable at the beginning
   !    - an optional argument can be supplied to transpose data whose
   !      global size is not the default nx*ny*nz
   !       * as the case in fft r2c/c2r interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface decomp_2d_init
      module procedure decomp_2d_init_ref
   end interface decomp_2d_init

   interface decomp_2d_finalize
      module procedure decomp_2d_finalize_ref
   end interface decomp_2d_finalize

   interface transpose_x_to_y
      module subroutine transpose_x_to_y_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_x_to_y_real_long
      module subroutine transpose_x_to_y_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_x_to_y_real_short
      module subroutine transpose_x_to_y_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_x_to_y_complex_long
      module subroutine transpose_x_to_y_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_x_to_y_complex_short
   end interface transpose_x_to_y

   interface transpose_y_to_z
      module subroutine transpose_y_to_z_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_z_real_long
      module subroutine transpose_y_to_z_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_z_real_short
      module subroutine transpose_y_to_z_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_z_complex_long
      module subroutine transpose_y_to_z_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_z_complex_short
   end interface transpose_y_to_z

   interface transpose_z_to_y
      module subroutine transpose_z_to_y_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_z_to_y_real_long
      module subroutine transpose_z_to_y_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_z_to_y_real_short
      module subroutine transpose_z_to_y_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_z_to_y_complex_long
      module subroutine transpose_z_to_y_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_z_to_y_complex_short
   end interface transpose_z_to_y

   interface transpose_y_to_x
      module subroutine transpose_y_to_x_real_long(src, dst, decomp)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_x_real_long
      module subroutine transpose_y_to_x_real_short(src, dst)
         real(mytype), dimension(:, :, :), intent(IN) :: src
         real(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_x_real_short
      module subroutine transpose_y_to_x_complex_long(src, dst, decomp)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
         TYPE(DECOMP_INFO), intent(IN) :: decomp
      end subroutine transpose_y_to_x_complex_long
      module subroutine transpose_y_to_x_complex_short(src, dst)
         complex(mytype), dimension(:, :, :), intent(IN) :: src
         complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      end subroutine transpose_y_to_x_complex_short
   end interface transpose_y_to_x

   interface alloc_x
      module procedure alloc_x_freal
      module procedure alloc_x_freal_short
      module procedure alloc_x_dreal
      module procedure alloc_x_dreal_short
      module procedure alloc_x_fcplx
      module procedure alloc_x_fcplx_short
      module procedure alloc_x_dcplx
      module procedure alloc_x_dcplx_short
      module procedure alloc_x_ints
      module procedure alloc_x_ints_short
      module procedure alloc_x_logs
      module procedure alloc_x_logs_short
#if defined(_GPU)
      module procedure alloc_x_dev_freal
      module procedure alloc_x_dev_freal_short
      module procedure alloc_x_dev_dreal
      module procedure alloc_x_dev_dreal_short
      module procedure alloc_x_dev_fcplx
      module procedure alloc_x_dev_fcplx_short
      module procedure alloc_x_dev_dcplx
      module procedure alloc_x_dev_dcplx_short
      module procedure alloc_x_dev_ints
      module procedure alloc_x_dev_ints_short
      module procedure alloc_x_dev_logs
      module procedure alloc_x_dev_logs_short
#endif
   end interface alloc_x

   interface alloc_y
      module procedure alloc_y_freal
      module procedure alloc_y_freal_short
      module procedure alloc_y_dreal
      module procedure alloc_y_dreal_short
      module procedure alloc_y_fcplx
      module procedure alloc_y_fcplx_short
      module procedure alloc_y_dcplx
      module procedure alloc_y_dcplx_short
      module procedure alloc_y_ints
      module procedure alloc_y_ints_short
      module procedure alloc_y_logs
      module procedure alloc_y_logs_short
#if defined(_GPU)
      module procedure alloc_y_dev_freal
      module procedure alloc_y_dev_freal_short
      module procedure alloc_y_dev_dreal
      module procedure alloc_y_dev_dreal_short
      module procedure alloc_y_dev_fcplx
      module procedure alloc_y_dev_fcplx_short
      module procedure alloc_y_dev_dcplx
      module procedure alloc_y_dev_dcplx_short
      module procedure alloc_y_dev_ints
      module procedure alloc_y_dev_ints_short
      module procedure alloc_y_dev_logs
      module procedure alloc_y_dev_logs_short
#endif
   end interface alloc_y

   interface alloc_z
      module procedure alloc_z_freal
      module procedure alloc_z_freal_short
      module procedure alloc_z_dreal
      module procedure alloc_z_dreal_short
      module procedure alloc_z_fcplx
      module procedure alloc_z_fcplx_short
      module procedure alloc_z_dcplx
      module procedure alloc_z_dcplx_short
      module procedure alloc_z_ints
      module procedure alloc_z_ints_short
      module procedure alloc_z_logs
      module procedure alloc_z_logs_short
#if defined(_GPU)
      module procedure alloc_z_dev_freal
      module procedure alloc_z_dev_freal_short
      module procedure alloc_z_dev_dreal
      module procedure alloc_z_dev_dreal_short
      module procedure alloc_z_dev_fcplx
      module procedure alloc_z_dev_fcplx_short
      module procedure alloc_z_dev_dcplx
      module procedure alloc_z_dev_dcplx_short
      module procedure alloc_z_dev_ints
      module procedure alloc_z_dev_ints_short
      module procedure alloc_z_dev_logs
      module procedure alloc_z_dev_logs_short
#endif
   end interface alloc_z

   interface

      module function d2d_log_is_active()
         logical :: d2d_log_is_active
      end function d2d_log_is_active

      module function d2d_log_get_unit()
         integer :: d2d_log_get_unit
      end function d2d_log_get_unit

      module subroutine d2d_log_close_unit(io_unit)
         integer, intent(in) :: io_unit
      end subroutine d2d_log_close_unit

      module subroutine d2d_log(given_io_unit)
         integer, intent(in), optional :: given_io_unit
      end subroutine d2d_log

      module subroutine decomp_info_print(d2d, io_unit, d2dname)
         type(decomp_info), intent(in) :: d2d
         integer, intent(in) :: io_unit
         character(len=*), intent(in) :: d2dname
      end subroutine decomp_info_print

   end interface

contains

#include "decomp_2d_init_fin.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Return the default decomposition object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! FIXME avoid a copy and return a pointer to decomp_main
   ! TODO list the external codes using this subroutine
   subroutine get_decomp_info(decomp)

      implicit none

      ! FIXME TYPE(DECOMP_INFO), pointer :: decomp
      TYPE(DECOMP_INFO), intent(OUT) :: decomp

      ! FIXME decomp => decomp_main
      decomp = decomp_main

      return
   end subroutine get_decomp_info

   !
   ! Return the 2D processor grid
   !
   function get_decomp_dims()

      implicit none

      integer, dimension(2) :: get_decomp_dims

      get_decomp_dims = dims

   end function get_decomp_dims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advanced Interface allowing applications to define globle domain of
   ! any size, distribute it, and then transpose data among pencils.
   !  - generate 2D decomposition details as defined in DECOMP_INFO
   !  - the default global data size is nx*ny*nz
   !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
   !  - multiple global sizes can co-exist in one application, each
   !    using its own DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_info_init(nx, ny, nz, decomp)

      use, intrinsic:: iso_c_binding, only: c_f_pointer, c_loc

      implicit none

      ! Arguments
      integer, intent(IN) :: nx, ny, nz
      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      ! Local variables
      integer(c_size_t) :: buf_size
      integer :: errorcode
#if defined(_GPU)
      integer :: status
#endif

      ! verify the global size can actually be distributed as pencils
      if (nx_global < dims(1) .or. ny_global < dims(1) .or. ny_global < dims(2) .or. nz_global < dims(2)) then
         errorcode = 6
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Invalid 2D processor grid. '// &
                              'Make sure that min(nx,ny) >= p_row and '// &
                              'min(ny,nz) >= p_col')
      end if

      ! distribute mesh points
      allocate (decomp%x1dist(0:dims(1) - 1), decomp%y1dist(0:dims(1) - 1), &
                decomp%y2dist(0:dims(2) - 1), decomp%z2dist(0:dims(2) - 1))
      call get_dist(nx, ny, nz, decomp)

      ! generate partition information - starting/ending index etc.
      call partition(nx, ny, nz, (/1, 2, 3/), &
                     decomp%xst, decomp%xen, decomp%xsz)
      call partition(nx, ny, nz, (/2, 1, 3/), &
                     decomp%yst, decomp%yen, decomp%ysz)
      call partition(nx, ny, nz, (/2, 3, 1/), &
                     decomp%zst, decomp%zen, decomp%zsz)

      ! prepare send/receive buffer displacement and count for ALLTOALL(V)
      allocate (decomp%x1cnts(0:dims(1) - 1), decomp%y1cnts(0:dims(1) - 1), &
                decomp%y2cnts(0:dims(2) - 1), decomp%z2cnts(0:dims(2) - 1))
      allocate (decomp%x1disp(0:dims(1) - 1), decomp%y1disp(0:dims(1) - 1), &
                decomp%y2disp(0:dims(2) - 1), decomp%z2disp(0:dims(2) - 1))
      call prepare_buffer(decomp)

      ! Update the shared memory pool
      if (use_pool) then
         if (decomp_pool_ready) then
            call decomp_pool%new_shape(decomp_pool_default_type, decomp)
         else
            call decomp_pool_init(decomp, decomp_pool_default_type)
         end if
      end if

      ! allocate memory for the MPI_ALLTOALL(V) buffers
      ! define the buffers globally for performance reason

      buf_size = max(product(int(decomp%xsz(1:3), c_size_t)), &
                     max(product(int(decomp%ysz(1:3), c_size_t)), &
                         product(int(decomp%zsz(1:3), c_size_t))))

#ifdef EVEN
      ! padded alltoall optimisation may need larger buffer space
      buf_size = max(buf_size, &
                     max(int(decomp%x1count, c_size_t) * int(dims(1), c_size_t), &
                         int(decomp%y2count, c_size_t) * int(dims(2), c_size_t)))
      ! evenly distributed data ?
      if (mod(nx, dims(1)) == 0 .and. mod(ny, dims(1)) == 0 .and. &
          mod(ny, dims(2)) == 0 .and. mod(nz, dims(2)) == 0) then
         decomp%even = .true.
      else
         decomp%even = .false.
      end if
#endif

      ! check if additional memory is required
      if (buf_size > decomp_buf_size) then
         decomp_buf_size = buf_size
#if defined(_GPU)
         if (.not.use_pool) then
            if (associated(work1_r)) nullify (work1_r)
            if (associated(work2_r)) nullify (work2_r)
            if (associated(work1_c)) nullify (work1_c)
            if (associated(work2_c)) nullify (work2_c)
            if (allocated(work1)) deallocate (work1)
            if (allocated(work2)) deallocate (work2)
            allocate (work1(2_c_size_t * buf_size), STAT=status)
            if (status /= 0) then
               errorcode = 2
               call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                    'Out of memory when allocating 2DECOMP workspace')
            end if
            allocate (work2(2_c_size_t * buf_size), STAT=status)
            if (status /= 0) then
               errorcode = 2
               call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                    'Out of memory when allocating 2DECOMP workspace')
            end if
            call c_f_pointer(c_loc(work1), work1_r, [buf_size])
            call c_f_pointer(c_loc(work2), work2_r, [buf_size])
            call c_f_pointer(c_loc(work1), work1_c, [buf_size])
            call c_f_pointer(c_loc(work2), work2_c, [buf_size])
            call decomp_2d_cumpi_init(buf_size, work1, work2)
         end if
#endif
      end if

   end subroutine decomp_info_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Release memory associated with a DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_info_finalize(decomp)

      implicit none

      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      if (allocated(decomp%x1dist)) deallocate (decomp%x1dist)
      if (allocated(decomp%y1dist)) deallocate (decomp%y1dist)
      if (allocated(decomp%y2dist)) deallocate (decomp%y2dist)
      if (allocated(decomp%z2dist)) deallocate (decomp%z2dist)
      if (allocated(decomp%x1cnts)) deallocate (decomp%x1cnts)
      if (allocated(decomp%y1cnts)) deallocate (decomp%y1cnts)
      if (allocated(decomp%y2cnts)) deallocate (decomp%y2cnts)
      if (allocated(decomp%z2cnts)) deallocate (decomp%z2cnts)
      if (allocated(decomp%x1disp)) deallocate (decomp%x1disp)
      if (allocated(decomp%y1disp)) deallocate (decomp%y1disp)
      if (allocated(decomp%y2disp)) deallocate (decomp%y2disp)
      if (allocated(decomp%z2disp)) deallocate (decomp%z2disp)

      return
   end subroutine decomp_info_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Find sub-domain information held by current processor
   !   INPUT:
   !     nx, ny, nz - global data dimension
   !     pdim(3)    - number of processor grid in each dimension,
   !                  valid values: 1 - distibute locally;
   !                                2 - distribute across p_row;
   !                                3 - distribute across p_col
   !   OUTPUT:
   !     lstart(3)  - starting index
   !     lend(3)    - ending index
   !     lsize(3)   - size of the sub-block (redundant)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

      implicit none

      integer, intent(IN) :: nx, ny, nz
      integer, dimension(3), intent(IN) :: pdim
      integer, dimension(3), intent(OUT) :: lstart, lend, lsize

      integer, allocatable, dimension(:) :: st, en, sz
      integer :: i, gsize

      do i = 1, 3

         if (i == 1) then
            gsize = nx
         else if (i == 2) then
            gsize = ny
         else if (i == 3) then
            gsize = nz
         end if

         if (pdim(i) == 1) then        ! all local
            lstart(i) = 1
            lend(i) = gsize
            lsize(i) = gsize
         elseif (pdim(i) == 2) then    ! distribute across dims(1)
            allocate (st(0:dims(1) - 1))
            allocate (en(0:dims(1) - 1))
            allocate (sz(0:dims(1) - 1))
            call distribute(gsize, dims(1), st, en, sz)
            lstart(i) = st(coord(1))
            lend(i) = en(coord(1))
            lsize(i) = sz(coord(1))
            deallocate (st, en, sz)
         elseif (pdim(i) == 3) then    ! distribute across dims(2)
            allocate (st(0:dims(2) - 1))
            allocate (en(0:dims(2) - 1))
            allocate (sz(0:dims(2) - 1))
            call distribute(gsize, dims(2), st, en, sz)
            lstart(i) = st(coord(2))
            lend(i) = en(coord(2))
            lsize(i) = sz(coord(2))
            deallocate (st, en, sz)
         end if

      end do
      return

   end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   - distibutes grid points in one dimension
   !   - handles uneven distribution properly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine distribute(data1, proc, st, en, sz)

      implicit none
      ! data1 -- data size in any dimension to be partitioned
      ! proc  -- number of processors in that dimension
      ! st    -- array of starting index
      ! en    -- array of ending index
      ! sz    -- array of local size  (redundent)
      integer data1, proc, st(0:proc - 1), en(0:proc - 1), sz(0:proc - 1)
      integer i, size1, nl, nu

      size1 = data1 / proc
      nu = data1 - size1 * proc
      nl = proc - nu
      st(0) = 1
      sz(0) = size1
      en(0) = size1
      do i = 1, nl - 1
         st(i) = st(i - 1) + size1
         sz(i) = size1
         en(i) = en(i - 1) + size1
      end do
      size1 = size1 + 1
      do i = nl, proc - 1
         st(i) = en(i - 1) + 1
         sz(i) = size1
         en(i) = en(i - 1) + size1
      end do

      ! Safety checks
      if (en(proc - 1) /= data1) &
         call decomp_2d_abort(__FILE__, __LINE__, en(proc - 1), "Invalid distribution.")
      if (sz(proc - 1) /= (data1 - st(proc - 1) + 1)) &
         call decomp_2d_abort(__FILE__, __LINE__, sz(proc - 1), "Invalid distribution.")

   end subroutine distribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  Define how each dimension is distributed across processors
   !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
   !    such global information is required locally at MPI_ALLTOALLV time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_dist(nx, ny, nz, decomp)

      implicit none

      integer, intent(IN) :: nx, ny, nz
      TYPE(DECOMP_INFO), intent(INOUT) :: decomp
      integer, allocatable, dimension(:) :: st, en

      allocate (st(0:dims(1) - 1))
      allocate (en(0:dims(1) - 1))
      call distribute(nx, dims(1), st, en, decomp%x1dist)
      call distribute(ny, dims(1), st, en, decomp%y1dist)
      deallocate (st, en)

      allocate (st(0:dims(2) - 1))
      allocate (en(0:dims(2) - 1))
      call distribute(ny, dims(2), st, en, decomp%y2dist)
      call distribute(nz, dims(2), st, en, decomp%z2dist)
      deallocate (st, en)

      return
   end subroutine get_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine prepare_buffer(decomp)

      implicit none

      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      integer :: i

      ! MPI_ALLTOALLV buffer information
      ! Avoid loops on arrays starting 0 (and not 1)
      ! with elements that depends on i-1 => Issues with O3 optimiser
      ! Solution explicitly set element (0) to start with

      do i = 0, dims(1) - 1
         decomp%x1cnts(i) = decomp%x1dist(i) * decomp%xsz(2) * decomp%xsz(3)
         decomp%y1cnts(i) = decomp%ysz(1) * decomp%y1dist(i) * decomp%ysz(3)
      end do

      decomp%x1disp(0) = 0  ! displacement is 0-based index
      decomp%y1disp(0) = 0

      do i = 1, dims(1) - 1
         decomp%x1disp(i) = decomp%x1disp(i - 1) + decomp%x1cnts(i - 1)
         decomp%y1disp(i) = decomp%y1disp(i - 1) + decomp%y1cnts(i - 1)
      end do

      do i = 0, dims(2) - 1
         decomp%y2cnts(i) = decomp%ysz(1) * decomp%y2dist(i) * decomp%ysz(3)
         decomp%z2cnts(i) = decomp%zsz(1) * decomp%zsz(2) * decomp%z2dist(i)
      end do

      decomp%y2disp(0) = 0  ! displacement is 0-based index
      decomp%z2disp(0) = 0

      do i = 1, dims(2) - 1
         decomp%y2disp(i) = decomp%y2disp(i - 1) + decomp%y2cnts(i - 1)
         decomp%z2disp(i) = decomp%z2disp(i - 1) + decomp%z2cnts(i - 1)
      end do

      ! MPI_ALLTOALL buffer information
#ifdef EVEN
      ! For evenly distributed data, following is an easier implementation.
      ! But it should be covered by the more general formulation below.
      !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
      !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1)
      !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
      !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

      ! For unevenly distributed data, pad smaller messages. Note the
      ! last blocks along pencils always get assigned more mesh points
      ! for X <=> Y transposes
      decomp%x1count = decomp%x1dist(dims(1) - 1) * &
                       decomp%y1dist(dims(1) - 1) * decomp%xsz(3)
      decomp%y1count = decomp%x1count
      ! for Y <=> Z transposes
      decomp%y2count = decomp%y2dist(dims(2) - 1) * &
                       decomp%z2dist(dims(2) - 1) * decomp%zsz(1)
      decomp%z2count = decomp%y2count
#endif

   end subroutine prepare_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"
#if defined(_GPU)
#include "alloc_dev.f90"
#endif

end module decomp_2d
