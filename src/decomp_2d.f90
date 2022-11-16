!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
! Copyright (C) 2021               the University of Edinburgh (UoE)
!
!=======================================================================

! This is the main 2D pencil decomposition module

module decomp_2d

  use MPI
  use, intrinsic :: iso_fortran_env, only : real32, real64
#if defined(_GPU)
  use cudafor
#if defined(_NCCL)
  use nccl
#endif
#endif

  implicit none

  private        ! Make everything private unless declared public

#ifdef DOUBLE_PREC
  integer, parameter, public :: mytype = KIND(0._real64)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: real2_type = MPI_2DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef SAVE_SINGLE
  integer, parameter, public :: mytype_single = KIND(0._real32)
  integer, parameter, public :: real_type_single = MPI_REAL
#else
  integer, parameter, public :: mytype_single = KIND(0._real64)
  integer, parameter, public :: real_type_single = MPI_DOUBLE_PRECISION
#endif
#else
  integer, parameter, public :: mytype = KIND(0._real32)
  integer, parameter, public :: real_type = MPI_REAL
  integer, parameter, public :: real2_type = MPI_2REAL
  integer, parameter, public :: complex_type = MPI_COMPLEX
  integer, parameter, public :: mytype_single = KIND(0._real32)
  integer, parameter, public :: real_type_single = MPI_REAL
#endif

  integer, save, public :: mytype_bytes

  ! some key global variables
  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: nrank = -1 ! local MPI rank 
  integer, save, public :: nproc = -1 ! total number of processors
  integer, save, public :: decomp_2d_comm = MPI_COMM_NULL ! MPI communicator

  ! parameters for 2D Cartesian topology 
  integer, save, dimension(2) :: dims, coord
  integer, save, public :: DECOMP_2D_COMM_CART_X = MPI_COMM_NULL
  integer, save, public :: DECOMP_2D_COMM_CART_Y = MPI_COMM_NULL
  integer, save, public :: DECOMP_2D_COMM_CART_Z = MPI_COMM_NULL
  integer, save :: DECOMP_2D_COMM_ROW = MPI_COMM_NULL
  integer, save :: DECOMP_2D_COMM_COL = MPI_COMM_NULL

  ! define neighboring blocks (to be used in halo-cell support)
  !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, save, dimension(3,6) :: neighbour 

  ! flags for periodic condition in three dimensions
  logical, save :: periodic_x, periodic_y, periodic_z

  !
  ! Debug level can be changed by the external code before calling decomp_2d_init
  !
  ! The environment variable "DECOMP_2D_DEBUG" can be used to change the debug level
  !
  ! Debug checks are performed only when the preprocessor variable DEBUG is defined
  !
  enum, bind(c)
     enumerator :: D2D_DEBUG_LEVEL_OFF = 0
     enumerator :: D2D_DEBUG_LEVEL_CRITICAL = 1
     enumerator :: D2D_DEBUG_LEVEL_ERROR = 2
     enumerator :: D2D_DEBUG_LEVEL_WARN = 3
     enumerator :: D2D_DEBUG_LEVEL_INFO = 4
     enumerator :: D2D_DEBUG_LEVEL_DEBUG = 5
     enumerator :: D2D_DEBUG_LEVEL_TRACE = 6
  end enum
#ifdef DEBUG
  integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_INFO
#else
  integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_OFF
#endif

#if defined(_GPU)
#if defined(_NCCL)
  integer, save :: row_rank, col_rank
#endif
#endif

#ifdef SHM
  ! derived type to store shared-memory info
  TYPE, public :: SMP_INFO
     integer MPI_COMM          ! SMP associated with this communicator
     integer NODE_ME           ! rank in this communicator
     integer NCPU              ! size of this communicator
     integer SMP_COMM          ! communicator for SMP-node masters
     integer CORE_COMM         ! communicator for cores on SMP-node
     integer SMP_ME            ! SMP-node id starting from 1 ... NSMP
     integer NSMP              ! number of SMP-nodes in this communicator
     integer CORE_ME           ! core id starting from 1 ... NCORE
     integer NCORE             ! number of cores on this SMP-node
     integer MAXCORE           ! maximum no. cores on any SMP-node
     integer N_SND             ! size of SMP shared memory buffer
     integer N_RCV             ! size of SMP shared memory buffer
     integer(8) SND_P          ! SNDBUF address (cray pointer), for real 
     integer(8) RCV_P          ! RCVBUF address (cray pointer), for real
     integer(8) SND_P_c        ! for complex
     integer(8) RCV_P_c        ! for complex
  END TYPE SMP_INFO
#endif

  ! derived type to store decomposition info for a given global data size
  TYPE, public :: DECOMP_INFO
     ! staring/ending index and size of data held by current processor
     integer, dimension(3) :: xst, xen, xsz  ! x-pencil
     integer, dimension(3) :: yst, yen, ysz  ! y-pencil
     integer, dimension(3) :: zst, zen, zsz  ! z-pencil

     ! in addition to local information, processors also need to know 
     ! some global information for global communications to work 

     ! how each dimension is distributed along pencils
     integer, allocatable, dimension(:) :: &
          x1dist, y1dist, y2dist, z2dist

     ! send/receive buffer counts and displacements for MPI_ALLTOALLV
     integer, allocatable, dimension(:) :: &
          x1cnts, y1cnts, y2cnts, z2cnts
     integer, allocatable, dimension(:) :: &
          x1disp, y1disp, y2disp, z2disp

     ! buffer counts for MPI_ALLTOALL: either for evenly distributed data
     ! or for padded-alltoall
     integer :: x1count, y1count, y2count, z2count

     ! evenly distributed data
     logical :: even

#ifdef SHM
     ! For shared-memory implementation

     ! one instance of this derived type for each communicator
     ! shared moemory info, such as which MPI rank belongs to which node
     TYPE(SMP_INFO) :: ROW_INFO, COL_INFO

     ! shared send/recv buffers for ALLTOALLV
     integer, allocatable, dimension(:) :: x1cnts_s, y1cnts_s, &
          y2cnts_s, z2cnts_s
     integer, allocatable, dimension(:) :: x1disp_s, y1disp_s, &
          y2disp_s, z2disp_s
     ! A copy of original buffer displacement (will be overwriten)
     integer, allocatable, dimension(:) :: x1disp_o, y1disp_o, &
          y2disp_o, z2disp_o
#endif
  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), target, save, public :: decomp_main
  ! FIXME The extra decomp_info objects should be defined in the external code, not here
  !       Currently keeping them to avoid breaking external codes
  TYPE(DECOMP_INFO), save, public :: phG,ph1,ph2,ph3,ph4

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

  ! These are the buffers used by MPI_ALLTOALL(V) calls
  real(mytype),    allocatable, dimension(:) :: work1_r, work2_r
  complex(mytype), allocatable, dimension(:) :: work1_c, work2_c

#if defined(_GPU)
  real(mytype), allocatable, dimension(:), device :: work1_r_d, work2_r_d
  complex(mytype), allocatable, dimension(:), device :: work1_c_d, work2_c_d

#if defined(_NCCL)
  integer col_comm_size, row_comm_size
  integer, allocatable, dimension(:) :: local_to_global_col, local_to_global_row
  type(ncclUniqueId) :: nccl_uid_2decomp
  type(ncclComm) :: nccl_comm_2decomp
  integer(kind=cuda_stream_kind) :: cuda_stream_2decomp
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To define smaller arrays using every several mesh points
  integer, save, dimension(3), public :: xszS,yszS,zszS,xstS,ystS,zstS,xenS,yenS,zenS
  integer, save, dimension(3), public :: xszV,yszV,zszV,xstV,ystV,zstV,xenV,yenV,zenV
  integer, save, dimension(3), public :: xszP,yszP,zszP,xstP,ystP,zstP,xenP,yenP,zenP
  logical, save :: coarse_mesh_starts_from_1
  integer, save :: iskipS, jskipS, kskipS
  integer, save :: iskipV, jskipV, kskipV
  integer, save :: iskipP, jskipP, kskipP

  !
  ! Profiler section
  !
  ! Integer to select the profiling tool
  !    0 => no profiling, default
  !    1 => Caliper (https://github.com/LLNL/Caliper)
  !
  enum, bind(c)
    enumerator :: decomp_profiler_none = 0
    enumerator :: decomp_profiler_caliper = 1
  end enum
  integer(kind(decomp_profiler_none)), save, public :: decomp_profiler = decomp_profiler_none
  ! Default : profile everything
  logical, save, public :: decomp_profiler_transpose = .true.
  logical, save, public :: decomp_profiler_io = .true.
  logical, save, public :: decomp_profiler_fft = .true.
  logical, save, public :: decomp_profiler_d2d = .true.

  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       decomp_info_init, decomp_info_finalize, partition, &
       decomp_info_print, decomp_profiler_prep, &
       decomp_profiler_start, decomp_profiler_end, &
       init_coarser_mesh_statS,fine_to_coarseS,&
       init_coarser_mesh_statV,fine_to_coarseV,&
       init_coarser_mesh_statP,fine_to_coarseP,&
       alloc_x, alloc_y, alloc_z, &
       update_halo, decomp_2d_abort, &
       decomp_2d_warning, get_decomp_info, &
       decomp_mpi_comm_free, get_decomp_dims 

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
     module subroutine decomp_2d_init_ref(nx,ny,nz,p_row,p_col,periodic_bc,comm)
       integer, intent(IN) :: nx,ny,nz
       integer, intent(INOUT) :: p_row,p_col
       logical, dimension(3), intent(IN), optional :: periodic_bc
       integer, intent(in), optional :: comm
     end subroutine decomp_2d_init_ref
  end interface decomp_2d_init
  interface decomp_2d_finalize
     module subroutine decomp_2d_finalize_ref()
     end subroutine decomp_2d_finalize_ref
  end interface decomp_2d_finalize
  
  interface transpose_x_to_y
     module procedure transpose_x_to_y_real
     module procedure transpose_x_to_y_real_short
     module procedure transpose_x_to_y_complex
     module procedure transpose_x_to_y_complex_short
  end interface transpose_x_to_y

  interface transpose_y_to_z
     module procedure transpose_y_to_z_real
     module procedure transpose_y_to_z_real_short
     module procedure transpose_y_to_z_complex
     module procedure transpose_y_to_z_complex_short
  end interface transpose_y_to_z

  interface transpose_z_to_y
     module procedure transpose_z_to_y_real
     module procedure transpose_z_to_y_real_short
     module procedure transpose_z_to_y_complex
     module procedure transpose_z_to_y_complex_short
  end interface transpose_z_to_y

  interface transpose_y_to_x
     module procedure transpose_y_to_x_real
     module procedure transpose_y_to_x_real_short
     module procedure transpose_y_to_x_complex
     module procedure transpose_y_to_x_complex_short
  end interface transpose_y_to_x

  interface update_halo
     module procedure update_halo_real
     module procedure update_halo_real_short
     module procedure update_halo_complex
     module procedure update_halo_complex_short
  end interface update_halo

  interface alloc_x
     module procedure alloc_x_real
     module procedure alloc_x_real_short
     module procedure alloc_x_complex
     module procedure alloc_x_complex_short
  end interface alloc_x

  interface alloc_y
     module procedure alloc_y_real
     module procedure alloc_y_real_short
     module procedure alloc_y_complex
     module procedure alloc_y_complex_short
  end interface alloc_y

  interface alloc_z
     module procedure alloc_z_real
     module procedure alloc_z_real_short
     module procedure alloc_z_complex
     module procedure alloc_z_complex_short
  end interface alloc_z

  interface decomp_2d_abort
     module procedure decomp_2d_abort_basic
     module procedure decomp_2d_abort_file_line
#if defined(_GPU) && defined(_NCCL)
     module procedure decomp_2d_abort_nccl_basic
     module procedure decomp_2d_abort_nccl_file_line
#endif
  end interface decomp_2d_abort

  interface decomp_2d_warning
     module procedure decomp_2d_warning_basic
     module procedure decomp_2d_warning_file_line
  end interface decomp_2d_warning

  interface

     module subroutine d2d_listing(given_io_unit)
        integer, intent(in), optional :: given_io_unit
     end subroutine d2d_listing

     module subroutine decomp_info_print(d2d, io_unit, d2dname)
        type(decomp_info), intent(in) :: d2d
        integer, intent(in) :: io_unit
        character(len=*), intent(in) :: d2dname
     end subroutine decomp_info_print

  end interface

   ! Generic interface to initialize the profiler
   interface decomp_profiler_init
      module subroutine decomp_profiler_init_noarg
      end subroutine decomp_profiler_init_noarg
   end interface decomp_profiler_init

   ! Generic interface to finalize the profiler
   interface decomp_profiler_fin
      module subroutine decomp_profiler_fin_noarg
      end subroutine decomp_profiler_fin_noarg
   end interface decomp_profiler_fin

   ! Generic interface for the profiler to log setup
   interface decomp_profiler_log
      module subroutine decomp_profiler_log_int(io_unit)
         integer, intent(in) :: io_unit
      end subroutine decomp_profiler_log_int
   end interface decomp_profiler_log

   ! Generic interface to prepare the profiler before init.
   interface decomp_profiler_prep
      module subroutine decomp_profiler_prep_bool(profiler_setup)
         logical, dimension(4), intent(in), optional :: profiler_setup
      end subroutine decomp_profiler_prep_bool
   end interface decomp_profiler_prep

   ! Generic interface for the profiler to start a given timer
   interface decomp_profiler_start
      module subroutine decomp_profiler_start_char(timer_name)
         character(len=*), intent(in) :: timer_name
      end subroutine decomp_profiler_start_char
   end interface decomp_profiler_start

   ! Generic interface for the profiler to end a given timer
   interface decomp_profiler_end
      module subroutine decomp_profiler_end_char(timer_name)
         character(len=*), intent(in) :: timer_name
      end subroutine decomp_profiler_end_char
   end interface decomp_profiler_end

   !---------------------------------------------------------------------
   ! Advanced Interface allowing applications to define globle domain of
   ! any size, distribute it, and then transpose data among pencils.
   !  - generate 2D decomposition details as defined in DECOMP_INFO
   !  - the default global data size is nx*ny*nz
   !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
   !  - multiple global sizes can co-exist in one application, each
   !    using its own DECOMP_INFO object
   !---------------------------------------------------------------------
   interface decomp_info_init
      module subroutine decomp_info_init_impl(nx,ny,nz,decomp)
        integer, intent(IN) :: nx,ny,nz
        TYPE(DECOMP_INFO), intent(INOUT) :: decomp
      end subroutine decomp_info_init_impl
   end interface decomp_info_init

   !----------------------------------------------------------------------
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
   !----------------------------------------------------------------------
   interface partition
      module subroutine partition_impl(nx, ny, nz, pdim, lstart, lend, lsize)
        integer, intent(IN) :: nx, ny, nz
        integer, dimension(3), intent(IN) :: pdim
        integer, dimension(3), intent(OUT) :: lstart, lend, lsize
      end subroutine partition_impl
   end interface partition
   
contains

  !
  ! Small wrapper to free a MPI communicator
  !
  subroutine decomp_mpi_comm_free(mpi_comm)

    implicit none

    integer, intent(inout) :: mpi_comm
    integer :: ierror

    ! Return if no MPI comm to free
    if (mpi_comm == MPI_COMM_NULL) return

    ! Free the provided MPI communicator
    call MPI_COMM_FREE(mpi_comm, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    mpi_comm = MPI_COMM_NULL

  end subroutine decomp_mpi_comm_free

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
  ! Release memory associated with a DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_finalize(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    if (allocated(decomp%x1dist)) deallocate(decomp%x1dist)
    if (allocated(decomp%y1dist)) deallocate(decomp%y1dist)
    if (allocated(decomp%y2dist)) deallocate(decomp%y2dist)
    if (allocated(decomp%z2dist)) deallocate(decomp%z2dist)
    if (allocated(decomp%x1cnts)) deallocate(decomp%x1cnts)
    if (allocated(decomp%y1cnts)) deallocate(decomp%y1cnts)
    if (allocated(decomp%y2cnts)) deallocate(decomp%y2cnts)
    if (allocated(decomp%z2cnts)) deallocate(decomp%z2cnts)
    if (allocated(decomp%x1disp)) deallocate(decomp%x1disp)
    if (allocated(decomp%y1disp)) deallocate(decomp%y1disp)
    if (allocated(decomp%y2disp)) deallocate(decomp%y2disp)
    if (allocated(decomp%z2disp)) deallocate(decomp%z2disp)

#ifdef SHM
    if (allocated(decomp%x1disp_o)) deallocate(decomp%x1disp_o)
    if (allocated(decomp%y1disp_o)) deallocate(decomp%y1disp_o)
    if (allocated(decomp%y2disp_o)) deallocate(decomp%y2disp_o)
    if (allocated(decomp%z2disp_o)) deallocate(decomp%z2disp_o)
    if (allocated(decomp%x1cnts_s)) deallocate(decomp%x1cnts_s)
    if (allocated(decomp%y1cnts_s)) deallocate(decomp%y1cnts_s)
    if (allocated(decomp%y2cnts_s)) deallocate(decomp%y2cnts_s)
    if (allocated(decomp%z2cnts_s)) deallocate(decomp%z2cnts_s)
    if (allocated(decomp%x1disp_s)) deallocate(decomp%x1disp_s)
    if (allocated(decomp%y1disp_s)) deallocate(decomp%y1disp_s)
    if (allocated(decomp%y2disp_s)) deallocate(decomp%y2disp_s)
    if (allocated(decomp%z2disp_s)) deallocate(decomp%z2disp_s)
#endif

    return
  end subroutine decomp_info_finalize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for statistic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statS(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipS = i_skip
    jskipS = j_skip
    kskipS = k_skip

    skip(1)=iskipS
    skip(2)=jskipS
    skip(3)=kskipS

    do i=1,3
       if (from1) then
          xstS(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstS(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = xend(i)/skip(i)
       end if
       xszS(i) = xenS(i)-xstS(i)+1
    end do

    do i=1,3
       if (from1) then
          ystS(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystS(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = yend(i)/skip(i)
       end if
       yszS(i) = yenS(i)-ystS(i)+1
    end do

    do i=1,3
       if (from1) then
          zstS(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstS(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = zend(i)/skip(i)
       end if
       zszS(i) = zenS(i)-zstS(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for visualization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statV(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipV = i_skip
    jskipV = j_skip
    kskipV = k_skip

    skip(1)=iskipV
    skip(2)=jskipV
    skip(3)=kskipV

    do i=1,3
       if (from1) then
          xstV(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstV(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = xend(i)/skip(i)
       end if
       xszV(i) = xenV(i)-xstV(i)+1
    end do

    do i=1,3
       if (from1) then
          ystV(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystV(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = yend(i)/skip(i)
       end if
       yszV(i) = yenV(i)-ystV(i)+1
    end do

    do i=1,3
       if (from1) then
          zstV(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstV(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = zend(i)/skip(i)
       end if
       zszV(i) = zenV(i)-zstV(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for probe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statP(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipP = i_skip
    jskipP = j_skip
    kskipP = k_skip

    skip(1)=iskipP
    skip(2)=jskipP
    skip(3)=kskipP

    do i=1,3
       if (from1) then
          xstP(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstP(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = xend(i)/skip(i)
       end if
       xszP(i) = xenP(i)-xstP(i)+1
    end do

    do i=1,3
       if (from1) then
          ystP(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystP(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = yend(i)/skip(i)
       end if
       yszP(i) = yenP(i)-ystP(i)+1
    end do

    do i=1,3
       if (from1) then
          zstP(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstP(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = zend(i)/skip(i)
       end if
       zszP(i) = zenP(i)-zstP(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statP

  ! Copy data from a fine-resolution array to a coarse one for statistic
  subroutine fine_to_coarseS(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystS(1):yenS(1),ystS(2):yenS(2),ystS(3):yenS(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstS(1):zenS(1),zstS(2):zenS(2),zstS(3):zenS(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseS

  ! Copy data from a fine-resolution array to a coarse one for visualization
  subroutine fine_to_coarseV(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystV(1):yenV(1),ystV(2):yenV(2),ystV(3):yenV(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstV(1):zenV(1),zstV(2):zenV(2),zstV(3):zenV(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseV

  ! Copy data from a fine-resolution array to a coarse one for probe
  subroutine fine_to_coarseP(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstP(3),xenP(3)
             do j=xstP(2),xenP(2)
                do i=xstP(1),xenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=xstP(3),xenP(3)
             do j=xstP(2),xenP(2)
                do i=xstP(1),xenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystP(1):yenP(1),ystP(2):yenP(2),ystP(3):yenP(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystP(3),yenP(3)
             do j=ystP(2),yenP(2)
                do i=ystP(1),yenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=ystP(3),yenP(3)
             do j=ystP(2),yenP(2)
                do i=ystP(1),yenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstP(1):zenP(1),zstP(2):zenP(2),zstP(3):zenP(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstP(3),zenP(3)
             do j=zstP(2),zenP(2)
                do i=zstP(1),zenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=zstP(3),zenP(3)
             do j=zstP(2),zenP(2)
                do i=zstP(1),zenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseP


#ifdef SHM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Generate shared-memory information 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init_shm(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! a copy of old displacement array (will be overwritten by shm code)
    allocate(decomp%x1disp_o(0:dims(1)-1),decomp%y1disp_o(0:dims(1)-1), &
         decomp%y2disp_o(0:dims(2)-1),decomp%z2disp_o(0:dims(2)-1))
    decomp%x1disp_o = decomp%x1disp
    decomp%y1disp_o = decomp%y1disp
    decomp%y2disp_o = decomp%y2disp
    decomp%z2disp_o = decomp%z2disp

    call prepare_shared_buffer(decomp%ROW_INFO,DECOMP_2D_COMM_ROW,decomp)
    call prepare_shared_buffer(decomp%COL_INFO,DECOMP_2D_COMM_COL,decomp)

    return
  end subroutine decomp_info_init_shm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For shared-memory implementation, prepare send/recv shared buffer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_shared_buffer(C,MPI_COMM,decomp)

    implicit none

    TYPE(SMP_INFO) :: C
    INTEGER :: MPI_COMM
    TYPE(DECOMP_INFO) :: decomp

    INTEGER, ALLOCATABLE :: KTBL(:,:),NARY(:,:),KTBLALL(:,:)
    INTEGER MYSMP, MYCORE, COLOR

    integer :: ierror

    C%MPI_COMM = MPI_COMM
    CALL MPI_COMM_SIZE(MPI_COMM,C%NCPU,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
    CALL MPI_COMM_RANK(MPI_COMM,C%NODE_ME,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_RANK")
    C%SMP_COMM  = MPI_COMM_NULL
    C%CORE_COMM = MPI_COMM_NULL
    C%SMP_ME= 0
    C%NCORE = 0
    C%CORE_ME = 0
    C%MAXCORE = 0
    C%NSMP  = 0
    C%N_SND = 0
    C%N_RCV = 0
    C%SND_P = 0
    C%RCV_P = 0
    C%SND_P_c = 0
    C%RCV_P_c = 0

    ! get smp-node map for this communicator and set up smp communicators
    CALL GET_SMP_MAP(C%MPI_COMM, C%NSMP, MYSMP, &
         C%NCORE, MYCORE, C%MAXCORE)
    C%SMP_ME = MYSMP + 1
    C%CORE_ME = MYCORE + 1
    ! - set up inter/intra smp-node communicators
    COLOR = MYCORE
    IF (COLOR.GT.0) COLOR = MPI_UNDEFINED
    CALL MPI_Comm_split(C%MPI_COMM, COLOR, MYSMP, C%SMP_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SPLIT")
    CALL MPI_Comm_split(C%MPI_COMM, MYSMP, MYCORE, C%CORE_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SPLIT")
    ! - allocate work space
    ALLOCATE(KTBL(C%MAXCORE,C%NSMP),NARY(C%NCPU,C%NCORE))
    ALLOCATE(KTBLALL(C%MAXCORE,C%NSMP))
    ! - set up smp-node/core to node_me lookup table
    KTBL = 0
    KTBL(C%CORE_ME,C%SMP_ME) = C%NODE_ME + 1
    CALL MPI_ALLREDUCE(KTBL,KTBLALL,C%NSMP*C%MAXCORE,MPI_INTEGER, &
         MPI_SUM,MPI_COMM,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")
    KTBL=KTBLALL
    !  IF (SUM(KTBL) /= C%NCPU*(C%NCPU+1)/2) &
    !       CALL MPI_ABORT(...

    ! compute offsets in shared SNDBUF and RCVBUF
    CALL MAPSET_SMPSHM(C, KTBL, NARY, decomp)

    DEALLOCATE(KTBL,NARY)

    return
  end subroutine prepare_shared_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Use Ian Bush's FreeIPC to generate shared-memory information
  !  - system independent solution
  !  - replacing David Tanqueray's implementation in alloc_shm.c
  !    (old C code renamed to get_smp_map2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_smp_map(comm, nnodes, my_node, ncores, my_core, maxcor)

    use FIPC_module

    implicit none

    integer, intent(IN) :: comm
    integer, intent(OUT) :: nnodes, my_node, ncores, my_core, maxcor

    integer :: intra_comm, extra_comm
    integer :: ierror

    call FIPC_init(comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "FIPC_init")

    ! intra_comm: communicator for processes on this shared memory node
    ! extra_comm: communicator for all rank 0 on each shared memory node
    call FIPC_ctxt_intra_comm(FIPC_ctxt_world, intra_comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "FIPC_ctxt_intra_comm")
    call FIPC_ctxt_extra_comm(FIPC_ctxt_world, extra_comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "FIPC_ctxt_extra_comm")

    call MPI_COMM_SIZE(intra_comm,  ncores, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
    call MPI_COMM_RANK(intra_comm, my_core, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")

    ! only rank 0 on each shared memory node member of extra_comm
    ! for others extra_comm = MPI_COMM_NULL
    if (extra_comm /= MPI_COMM_NULL) then
       call MPI_COMM_SIZE(extra_comm,  nnodes, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
       call MPI_COMM_RANK(extra_comm, my_node, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    end if

    ! other ranks share the same information as their leaders
    call MPI_BCAST( nnodes, 1, MPI_INTEGER, 0, intra_comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
    call MPI_BCAST(my_node, 1, MPI_INTEGER, 0, intra_comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

    ! maxcor
    call MPI_ALLREDUCE(ncores, maxcor, 1, MPI_INTEGER, MPI_MAX, &
         decomp_2d_comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")

    call FIPC_finalize(ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "FIPC_finalize")

    return

  end subroutine get_smp_map


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up smp-node based shared memory maps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MAPSET_SMPSHM(C, KTBL, NARY, decomp)

    IMPLICIT NONE

    TYPE (SMP_INFO) C
    INTEGER KTBL(C%MAXCORE,C%NSMP)
    INTEGER NARY(C%NCPU,C%NCORE)
    TYPE (DECOMP_INFO) :: decomp

    INTEGER i, j, k, l, N, PTR, BSIZ, ierror, status, seed
    character*16 s

    BSIZ = C%N_SND

    ! a - SNDBUF
    IF (C%MPI_COMM==DECOMP_2D_COMM_COL) THEN
       ALLOCATE(decomp%x1cnts_s(C%NSMP),decomp%x1disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%x1cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLGATHER")
       PTR = 0
       DO i=1,C%NSMP
          decomp%x1disp_s(i) = PTR
          N = 0
          DO j=1,C%MAXCORE
             k = KTBL(j,i)
             IF (k > 0) then
                DO l=1,C%NCORE
                   IF (l == C%CORE_ME) decomp%x1disp_o(k-1) = PTR
                   N = N + NARY(k,l)
                   PTR = PTR + NARY(k,l)
                END DO
             END IF
          END DO
          decomp%x1cnts_s(i) = N
       END DO
       decomp%x1disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR

    ELSE IF (C%MPI_COMM==DECOMP_2D_COMM_ROW) THEN
       ALLOCATE(decomp%y2cnts_s(C%NSMP),decomp%y2disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%y2cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLGATHER")
       PTR = 0
       DO i=1,C%NSMP
          decomp%y2disp_s(i) = PTR
          N = 0
          DO j=1,C%MAXCORE
             k = KTBL(j,i)
             IF (k > 0) then
                DO l=1,C%NCORE
                   IF (l == C%CORE_ME) decomp%y2disp_o(k-1) = PTR
                   N = N + NARY(k,l)
                   PTR = PTR + NARY(k,l)
                END DO
             END IF
          END DO
          decomp%y2cnts_s(i) = N
       END DO
       decomp%y2disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
    END IF

    ! b - RCVBUF

    IF (C%MPI_COMM==DECOMP_2D_COMM_COL) THEN
       ALLOCATE(decomp%y1cnts_s(C%NSMP),decomp%y1disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%y1cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLGATHER")
       PTR = 0
       DO i=1,C%NSMP
          decomp%y1disp_s(i) = PTR
          N=0
          DO j=1,C%NCORE
             DO l=1,C%MAXCORE
                k = KTBL(l,i)
                IF (k > 0) then
                   IF (j == C%CORE_ME) decomp%y1disp_o(k-1) = PTR
                   N = N + NARY(k,j)
                   PTR = PTR + NARY(k,j)
                END IF
             END DO
          END DO
          decomp%y1cnts_s(i) = N
       END DO
       decomp%y1disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR

    ELSE IF (C%MPI_COMM==DECOMP_2D_COMM_ROW) THEN
       ALLOCATE(decomp%z2cnts_s(C%NSMP),decomp%z2disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%z2cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLGATHER")
       PTR = 0
       DO i=1,C%NSMP
          decomp%z2disp_s(i) = PTR
          N=0
          DO j=1,C%NCORE
             DO l=1,C%MAXCORE
                k = KTBL(l,i)
                IF (k > 0) then
                   IF (j == C%CORE_ME) decomp%z2disp_o(k-1) = PTR
                   N = N + NARY(k,j)
                   PTR = PTR + NARY(k,j)
                END IF
             END DO
          END DO
          decomp%z2cnts_s(i) = N
       END DO
       decomp%z2disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR

    END IF

    ! check buffer size and (re)-allocate buffer space if necessary
    IF (BSIZ > C%N_SND) then
       IF (C%SND_P /= 0) CALL DEALLOC_SHM(C%SND_P, C%CORE_COMM)
       ! make sure each rank has unique keys to get shared memory
       !IF (C%MPI_COMM==DECOMP_2D_COMM_COL) THEN
       !   seed = nrank+nproc*0+1 ! has to be non-zero
       !ELSE IF (C%MPI_COMM==DECOMP_2D_COMM_ROW) THEN
       !   seed = nrank+nproc*1+1
       !END IF
       status = 1
       !CALL ALLOC_SHM(C%SND_P, BSIZ, real_type, C%CORE_COMM, status, &
       !     seed)
       CALL ALLOC_SHM(C%SND_P, BSIZ, real_type, C%CORE_COMM, status)
       C%N_SND = BSIZ

       IF (C%RCV_P /= 0) CALL DEALLOC_SHM(C%RCV_P, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%RCV_P, BSIZ, real_type, C%CORE_COMM, status)
       C%N_RCV = BSIZ

       IF (C%SND_P_c /= 0) CALL DEALLOC_SHM(C%SND_P_c, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%SND_P_c, BSIZ, complex_type, C%CORE_COMM, status)
       C%N_SND = BSIZ

       IF (C%RCV_P_c /= 0) CALL DEALLOC_SHM(C%RCV_P_c, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%RCV_P_c, BSIZ, complex_type, C%CORE_COMM, status)
       C%N_RCV = BSIZ


    END IF

    RETURN
  END SUBROUTINE MAPSET_SMPSHM

#endif


#ifdef OCC
  ! For non-blocking communication code, progress the comminication stack
  subroutine transpose_test(handle)

    implicit none

    integer :: handle, ierror

    call NBC_TEST(handle,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "NBC_TEST")

    return
  end subroutine transpose_test
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transposition routines 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "transpose_x_to_y.f90"
#include "transpose_y_to_z.f90"
#include "transpose_z_to_y.f90"
#include "transpose_y_to_x.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error handling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_abort_basic(errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    integer :: ierror

    if (nrank==0) then
       write(*,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
       write(error_unit,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(error_unit,*) 'ERROR MESSAGE: ' // msg
    end if
    call MPI_ABORT(decomp_2d_comm,errorcode,ierror)

  end subroutine decomp_2d_abort_basic

  subroutine decomp_2d_abort_file_line(file, line, errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode, line
    character(len=*), intent(IN) :: msg, file

    integer :: ierror

    if (nrank==0) then
       write(*,*) '2DECOMP&FFT ERROR'
       write(*,*) '  errorcode:     ', errorcode
       write(*,*) '  error in file  ' // file
       write(*,*) '           line  ', line
       write(*,*) '  error message: ' // msg
       write(error_unit,*) '2DECOMP&FFT ERROR'
       write(error_unit,*) '  errorcode:     ', errorcode
       write(error_unit,*) '  error in file  ' // file
       write(error_unit,*) '           line  ', line
       write(error_unit,*) '  error message: ' // msg
    end if
    call MPI_ABORT(decomp_2d_comm,errorcode,ierror)

  end subroutine decomp_2d_abort_file_line

#if defined(_GPU) && defined(_NCCL)
  !
  ! This is based on the file "nccl.h" in nvhpc 22.1
  !
  function _ncclresult_to_integer(errorcode)

    implicit none

    type(ncclresult), intent(IN) :: errorcode
    integer :: _ncclresult_to_integer

    if (errorcode == ncclSuccess) then
        _ncclresult_to_integer = 0
    elseif (errorcode == ncclUnhandledCudaError) then
        _ncclresult_to_integer = 1
    elseif (errorcode == ncclSystemError) then
        _ncclresult_to_integer = 2
    elseif (errorcode == ncclInternalError) then
        _ncclresult_to_integer = 3
    elseif (errorcode == ncclInvalidArgument) then
        _ncclresult_to_integer = 4
    elseif (errorcode == ncclInvalidUsage) then
        _ncclresult_to_integer = 5
    elseif (errorcode == ncclNumResults) then
        _ncclresult_to_integer = 6
    else
      _ncclresult_to_integer = -1
      call decomp_2d_warning(__FILE__, __LINE__, _ncclresult_to_integer, &
                             "NCCL error handling needs some update")
    end if

  end function _ncclresult_to_integer

  !
  ! Small wrapper for basic NCCL errors
  !
  subroutine decomp_2d_abort_nccl_basic(errorcode, msg)

    implicit none

    type(ncclresult), intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    call decomp_2d_abort(_ncclresult_to_integer(errorcode), &
                         msg // " " // ncclGetErrorString(errorcode))

  end subroutine decomp_2d_abort_nccl_basic

  !
  ! Small wrapper for NCCL errors
  !
  subroutine decomp_2d_abort_nccl_file_line(file, line, errorcode, msg)

    implicit none

    type(ncclresult), intent(IN) :: errorcode
    integer, intent(in) :: line
    character(len=*), intent(IN) :: msg, file

    call decomp_2d_abort(file, &
                         line, &
                         _ncclresult_to_integer(errorcode), &
                         msg // " " // ncclGetErrorString(errorcode))

  end subroutine decomp_2d_abort_nccl_file_line
#endif

  subroutine decomp_2d_warning_basic(errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    if (nrank==0) then
       write(*,*) '2DECOMP&FFT WARNING - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
       write(error_unit,*) '2DECOMP&FFT WARNING - errorcode: ', errorcode
       write(error_unit,*) 'ERROR MESSAGE: ' // msg
    end if

  end subroutine decomp_2d_warning_basic

  subroutine decomp_2d_warning_file_line(file, line, errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode, line
    character(len=*), intent(IN) :: msg, file

    if (nrank==0) then
       write(*,*) '2DECOMP&FFT WARNING'
       write(*,*) '  errorcode:     ', errorcode
       write(*,*) '  error in file  ' // file
       write(*,*) '           line  ', line
       write(*,*) '  error message: ' // msg
       write(error_unit,*) '2DECOMP&FFT WARNING'
       write(error_unit,*) '  errorcode:     ', errorcode
       write(error_unit,*) '  error in file  ' // file
       write(error_unit,*) '           line  ', line
       write(error_unit,*) '  error message: ' // msg
    end if

  end subroutine decomp_2d_warning_file_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"
  
end module decomp_2d

