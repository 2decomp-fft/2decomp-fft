!!!=======================================================================
!!! This is part of the 2DECOMP&FFT library
!!! 
!!! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
!!! decomposition. It also implements a highly scalable distributed
!!! three-dimensional Fast Fourier Transform (FFT).
!!!
!!! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!!! Copyright (C) 2022-              the Xcompact3d developers
!!!
!!!=======================================================================

submodule (decomp_2d) decomp_2d_init_fin

  use mpi
  use factor

  implicit none

  integer, save :: decomp_buf_size = 0

contains

  !======================================================================
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
  !======================================================================
  module subroutine decomp_2d_init_ref(nx,ny,nz,p_row,p_col,periodic_bc,comm)

    use iso_fortran_env, only : output_unit

    implicit none

    integer, intent(IN) :: nx,ny,nz
    integer, intent(INOUT) :: p_row,p_col
    logical, dimension(3), intent(IN), optional :: periodic_bc
    integer, intent(in), optional :: comm

    integer :: errorcode, ierror, row, col, iounit
#if defined(_GPU) && defined(_NCCL)
    integer :: cuda_stat
    type(ncclResult) :: nccl_stat
#endif
    logical, dimension(2) :: periodic
#ifdef DEBUG
    character(len=7) fname ! Sufficient for up to O(1M) ranks
#endif

#ifdef PROFILER
    ! Prepare the profiler if it was not already prepared
    if (decomp_profiler.eq.decomp_profiler_none) call decomp_profiler_prep()
    ! Start the profiler
    call decomp_profiler_init()
    ! Start the timer for decomp_2d_init
    if (decomp_profiler_d2d) call decomp_profiler_start("decomp_2d_init")
#endif

    ! Use the provided MPI communicator if present
    if (present(comm)) then
       decomp_2d_comm = comm
    else
       decomp_2d_comm = MPI_COMM_WORLD
    endif

    ! If the external code has not set nrank and nproc
    if (nrank == -1) then
       call MPI_COMM_RANK(decomp_2d_comm, nrank, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, &
                                             __LINE__, &
                                             ierror, &
                                             "MPI_COMM_RANK")
    endif
    if (nproc == -1) then
       call MPI_COMM_SIZE(decomp_2d_comm, nproc, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, &
                                             __LINE__, &
                                             ierror, &
                                             "MPI_COMM_SIZE")
    endif
#ifdef DEBUG
    ! Check if a modification of the debug level is needed
    call decomp_2d_debug()
#endif

    nx_global = nx
    ny_global = ny
    nz_global = nz

    if (present(periodic_bc)) then
       periodic_x = periodic_bc(1)
       periodic_y = periodic_bc(2)
       periodic_z = periodic_bc(3)
    else
       periodic_x = .false.
       periodic_y = .false.
       periodic_z = .false.
    end if

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, row, col)
       p_row = row
       p_col = col
    else
       if (nproc /= p_row*p_col) then
          errorcode = 1
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Invalid 2D processor grid - nproc /= p_row*p_col')
       else
          row = p_row
          col = p_col
       end if
    end if

    ! Create 2D Catersian topology
    ! Note that in order to support periodic B.C. in the halo-cell code,
    ! need to create multiple topology objects: DECOMP_2D_COMM_CART_?,
    ! corresponding to three pencil orientations. They contain almost
    ! identical topological information but allow different combinations
    ! of periodic conditions.
    dims(1) = row
    dims(2) = col
    periodic(1) = periodic_y
    periodic(2) = periodic_z
    call MPI_CART_CREATE(decomp_2d_comm,2,dims,periodic, &
         .false., &  ! do not reorder rank
         DECOMP_2D_COMM_CART_X, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
    periodic(1) = periodic_x
    periodic(2) = periodic_z
    call MPI_CART_CREATE(decomp_2d_comm,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Y, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
    periodic(1) = periodic_x
    periodic(2) = periodic_y
    call MPI_CART_CREATE(decomp_2d_comm,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Z, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")

    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_COORDS")

    ! derive communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
         DECOMP_2D_COMM_COL,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SUB")
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
         DECOMP_2D_COMM_ROW,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SUB")

    ! gather information for halo-cell support code
    call init_neighbour

    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,decomp_main)

    ! make a copy of the decomposition information associated with the
    ! default global size in these global variables so applications can
    ! use them to create data structures 
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz

    ! determine the number of bytes per float number
    ! do not use 'mytype' which is compiler dependent
    ! also possible to use inquire(iolength=...) 
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")

#ifdef EVEN
    if (nrank==0) write(*,*) 'Padded ALLTOALL optimisation on'
#endif 

#if defined(_GPU)
#if defined(_NCCL)
    call MPI_COMM_RANK(DECOMP_2D_COMM_COL,col_rank,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,row_rank,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    call MPI_COMM_SIZE(DECOMP_2D_COMM_COL,col_comm_size,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
    call MPI_COMM_SIZE(DECOMP_2D_COMM_ROW,row_comm_size,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")

    allocate(local_to_global_col(col_comm_size), local_to_global_row(row_comm_size))
    
    local_to_global_col(:) = 0
    local_to_global_row(:) = 0
    local_to_global_col(col_rank+1) = nrank
    local_to_global_row(row_rank+1) = nrank
    
    call mpi_allreduce(MPI_IN_PLACE,local_to_global_col,col_comm_size,MPI_INTEGER,MPI_SUM,DECOMP_2D_COMM_COL,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")
    call mpi_allreduce(MPI_IN_PLACE,local_to_global_row,row_comm_size,MPI_INTEGER,MPI_SUM,DECOMP_2D_COMM_ROW,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")

    if (nrank .eq. 0) then
       nccl_stat = ncclGetUniqueId(nccl_uid_2decomp)
       if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclGetUniqueId")
    end if
    call MPI_Bcast(nccl_uid_2decomp, int(sizeof(ncclUniqueId)), MPI_BYTE, 0, decomp_2d_comm, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

    nccl_stat = ncclCommInitRank(nccl_comm_2decomp, nproc, nccl_uid_2decomp, nrank)
    if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclCommInitRank")
    cuda_stat = cudaStreamCreate(cuda_stream_2decomp)
    if (cuda_stat /= 0) call decomp_2d_abort(__FILE__, __LINE__, cuda_stat, "cudaStreamCreate")
#endif
#endif

    !
    ! Select the IO unit for decomp_2d setup
    !
#ifdef DEBUG
    write(fname, "(I0)") nrank ! Adapt to magnitude of nrank
    open(newunit=iounit, file='decomp_2d_setup_'//trim(fname)//'.log', iostat=ierror)
#else
    if (nrank == 0) then
       open(newunit=iounit, file="decomp_2d_setup.log", iostat=ierror)
    else
       iounit = output_unit
       ierror = 0
    endif
#endif
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "Could not open log file")
    !
    ! Print the decomp_2d setup
    !
    call d2d_listing(iounit)
    !
    ! Close the IO unit if it was not stdout
    !
    if (iounit /= output_unit) then
       close(iounit, iostat=ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "Could not close log file")
    endif

#ifdef PROFILER
    ! Stop the timer for decomp_2d_init
    if (decomp_profiler_d2d) call decomp_profiler_end("decomp_2d_init")
#endif

    return
  end subroutine decomp_2d_init_ref

  !======================================================================
  ! Routine to be called by applications to clean things up
  !======================================================================
  module subroutine decomp_2d_finalize_ref
    
    implicit none

#if defined(_GPU) && defined(_NCCL)
    type(ncclResult) :: nccl_stat
#endif
 
#ifdef PROFILER
    if (decomp_profiler_d2d) call decomp_profiler_start("decomp_2d_fin")
#endif

    call decomp_mpi_comm_free(DECOMP_2D_COMM_ROW)
    call decomp_mpi_comm_free(DECOMP_2D_COMM_COL)
    call decomp_mpi_comm_free(DECOMP_2D_COMM_CART_X)
    call decomp_mpi_comm_free(DECOMP_2D_COMM_CART_Y)
    call decomp_mpi_comm_free(DECOMP_2D_COMM_CART_Z)

    call decomp_info_finalize(decomp_main)

    decomp_buf_size = 0
    deallocate(work1_r, work2_r, work1_c, work2_c)
#if defined(_GPU)
    deallocate(work1_r_d, work2_r_d, work1_c_d, work2_c_d)

#if defined(_NCCL)
    nccl_stat = ncclCommDestroy(nccl_comm_2decomp)
    if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclCommDestroy")
#endif
#endif

    nrank = -1
    nproc = -1

#ifdef PROFILER
    if (decomp_profiler_d2d) call decomp_profiler_end("decomp_2d_fin")
    ! Finalize the profiler
    call decomp_profiler_fin()
#endif

    return
  end subroutine decomp_2d_finalize_ref

#ifdef DEBUG
  !
  ! Try to read the environment variable DECOMP_2D_DEBUG to change the debug level
  !
  ! The expected value is an integer below 9999
  !
  subroutine decomp_2d_debug

    implicit none

    integer :: ierror
    character(len=4) :: val
    character(len=*), parameter :: varname = "DECOMP_2D_DEBUG"

    ! Read the variable
    call get_environment_variable(varname, value=val, status=ierror)

    ! Return if no variable, or no support for env. variable
    if (ierror >= 1) return

    ! Minor error, print warning and return
    if (ierror /= 0) then
       call decomp_2d_warning(__FILE__, &
                              __LINE__, &
                              ierror, &
                              "Error when reading DECOMP_2D_DEBUG : "//val)
       return
    endif

    ! Conversion to integer if possible
    read(val, '(i4)', iostat=ierror) decomp_debug
    if (ierror /= 0) call decomp_2d_warning(__FILE__, &
                                            __LINE__, &
                                            ierror, &
                                            "Error when reading DECOMP_2D_DEBUG : "//val)

  end subroutine decomp_2d_debug
#endif

  !---------------------------------------------------------------------
  ! Auto-tuning algorithm to select the best 2D processor grid
  !---------------------------------------------------------------------
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)

    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, allocatable, dimension(:) :: factors
    integer :: nfact, i, col, i_best

    if (nrank==0) write(*,*) 'In auto-tuning mode......'

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    i_best=nfact/2+1
    col=factors(i_best)

    best_p_col = col
    best_p_row=iproc/col
    if (nrank==0) print *,'p_row x p_col', best_p_row, best_p_col
    if ((best_p_col==1).and.(nrank==0)) then
       print *,'WARNING: current 2D DECOMP set-up might not work'
    endif
    
    deallocate(factors)

    return
  end subroutine best_2d_grid

  !---------------------------------------------------------------------
  ! To support halo-cell exchange:
  !   find the MPI ranks of neighbouring pencils
  ! XXX: belongs in a halo module
  !---------------------------------------------------------------------
  subroutine init_neighbour

    integer :: ierror

    ! For X-pencil
    neighbour(1,1) = MPI_PROC_NULL               ! east
    neighbour(1,2) = MPI_PROC_NULL               ! west
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
         neighbour(1,4), neighbour(1,3), ierror) ! north & south
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
         neighbour(1,6), neighbour(1,5), ierror) ! top & bottom
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

    ! For Y-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
         neighbour(2,2), neighbour(2,1), ierror) ! east & west
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    neighbour(2,3) = MPI_PROC_NULL               ! north
    neighbour(2,4) = MPI_PROC_NULL               ! south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
         neighbour(2,6), neighbour(2,5), ierror) ! top & bottom
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

    ! For Z-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, &
         neighbour(3,2), neighbour(3,1), ierror) ! east & west
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, &
         neighbour(3,4), neighbour(3,3), ierror) ! north & south
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
    neighbour(3,5) = MPI_PROC_NULL               ! top
    neighbour(3,6) = MPI_PROC_NULL               ! bottom
    return
  end subroutine init_neighbour

  !---------------------------------------------------------------------
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
  !---------------------------------------------------------------------
  module subroutine decomp_info_init_impl(nx,ny,nz,decomp)

    implicit none

    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: buf_size, status, errorcode

    ! verify the global size can actually be distributed as pencils
    if (nx_global<dims(1) .or. ny_global<dims(1) .or. ny_global<dims(2) .or. nz_global<dims(2)) then
       errorcode = 6
       call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
            'Invalid 2D processor grid. ' // &
            'Make sure that min(nx,ny) >= p_row and ' // &
            'min(ny,nz) >= p_col')
    end if

    if (mod(nx,dims(1))==0 .and. mod(ny,dims(1))==0 .and. &
         mod(ny,dims(2))==0 .and. mod(nz,dims(2))==0) then
       decomp%even = .true.
    else
       decomp%even = .false.
    end if

    ! distribute mesh points
    allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
         decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1))
    call get_dist(nx,ny,nz,decomp)

    ! generate partition information - starting/ending index etc.
    call partition(nx, ny, nz, (/ 1,2,3 /), &
         decomp%xst, decomp%xen, decomp%xsz)
    call partition(nx, ny, nz, (/ 2,1,3 /), &
         decomp%yst, decomp%yen, decomp%ysz)
    call partition(nx, ny, nz, (/ 2,3,1 /), &
         decomp%zst, decomp%zen, decomp%zsz)

    ! prepare send/receive buffer displacement and count for ALLTOALL(V)
    allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
         decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1))
    allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
         decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1))
    call prepare_buffer(decomp)

#ifdef SHM
    ! prepare shared-memory information if required
    call decomp_info_init_shm(decomp)
#endif

    ! allocate memory for the MPI_ALLTOALL(V) buffers
    ! define the buffers globally for performance reason

    buf_size = max(decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3), &
         max(decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3), &
         decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)) )
#ifdef EVEN
    ! padded alltoall optimisation may need larger buffer space
    buf_size = max(buf_size, &
         max(decomp%x1count*dims(1),decomp%y2count*dims(2)) ) 
#endif

    ! check if additional memory is required
    ! *** TODO: consider how to share the real/complex buffers 
    if (buf_size > decomp_buf_size) then
       decomp_buf_size = buf_size
#if defined(_GPU)
       if (allocated(work1_r_d)) deallocate(work1_r_d)
       if (allocated(work2_r_d)) deallocate(work2_r_d)
       if (allocated(work1_c_d)) deallocate(work1_c_d)
       if (allocated(work2_c_d)) deallocate(work2_c_d)
       allocate(work1_r_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work1_c_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_r_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_c_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
#endif
       if (allocated(work1_r)) deallocate(work1_r)
       if (allocated(work2_r)) deallocate(work2_r)
       if (allocated(work1_c)) deallocate(work1_c)
       if (allocated(work2_c)) deallocate(work2_c)
       allocate(work1_r(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_r(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work1_c(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_c(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
    end if

    return
  end subroutine decomp_info_init_impl

  !----------------------------------------------------------------------
  ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
  !----------------------------------------------------------------------
  subroutine prepare_buffer(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: i

    !LG : AJOUTS "bidons" pour eviter un plantage en -O3 avec gcc9.3
    !       * la fonction sortait des valeurs 'aleatoires'
    !         et le calcul plantait dans MPI_ALLTOALLV
    !       * pas de plantage en O2
    
    if (nrank==0) then
       open(newunit=i, file='temp.dat', form='unformatted')
       write(i) decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist, &
                decomp%xsz,decomp%ysz,decomp%zsz
       close(i, status='delete')
    endif

    ! MPI_ALLTOALLV buffer information

    do i=0, dims(1)-1
       decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
       if (i==0) then
          decomp%x1disp(i) = 0  ! displacement is 0-based index
          decomp%y1disp(i) = 0
       else
          decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
          decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
       end if
    end do

    do i=0, dims(2)-1
       decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
       decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
       if (i==0) then
          decomp%y2disp(i) = 0  ! displacement is 0-based index
          decomp%z2disp(i) = 0
       else
          decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
          decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
       end if
    end do

    ! MPI_ALLTOALL buffer information

    ! For evenly distributed data, following is an easier implementation.
    ! But it should be covered by the more general formulation below.
    !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
    !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1) 
    !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
    !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

    ! For unevenly distributed data, pad smaller messages. Note the 
    ! last blocks along pencils always get assigned more mesh points
    ! for X <=> Y transposes
    decomp%x1count = decomp%x1dist(dims(1)-1) * &
         decomp%y1dist(dims(1)-1) * decomp%xsz(3)
    decomp%y1count = decomp%x1count
    ! for Y <=> Z transposes
    decomp%y2count = decomp%y2dist(dims(2)-1) * &
         decomp%z2dist(dims(2)-1) * decomp%zsz(1)
    decomp%z2count = decomp%y2count

    return
  end subroutine prepare_buffer

  !----------------------------------------------------------------------
  !  Define how each dimension is distributed across processors
  !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
  !    such global information is required locally at MPI_ALLTOALLV time
  !----------------------------------------------------------------------
  subroutine get_dist(nx,ny,nz,decomp)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
    integer, allocatable, dimension(:) :: st,en

    allocate(st(0:dims(1)-1))
    allocate(en(0:dims(1)-1))
    call distribute(nx,dims(1),st,en,decomp%x1dist)
    call distribute(ny,dims(1),st,en,decomp%y1dist)
    deallocate(st,en)

    allocate(st(0:dims(2)-1))
    allocate(en(0:dims(2)-1))
    call distribute(ny,dims(2),st,en,decomp%y2dist)
    call distribute(nz,dims(2),st,en,decomp%z2dist)
    deallocate(st,en)

    return
  end subroutine get_dist

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
  subroutine partition_impl(nx, ny, nz, pdim, lstart, lend, lsize)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    integer, dimension(3), intent(IN) :: pdim
    integer, dimension(3), intent(OUT) :: lstart, lend, lsize

    integer, allocatable, dimension(:) :: st,en,sz
    integer :: i, gsize

    do i = 1, 3

       if (i==1) then
          gsize = nx
       else if (i==2) then
          gsize = ny
       else if (i==3) then
          gsize = nz
       end if

       if (pdim(i) == 1) then        ! all local
          lstart(i) = 1
          lend(i)   = gsize
          lsize(i)  = gsize
       elseif (pdim(i) == 2) then    ! distribute across dims(1)
          allocate(st(0:dims(1)-1))
          allocate(en(0:dims(1)-1))
          allocate(sz(0:dims(1)-1))
          call distribute(gsize,dims(1),st,en,sz)
          lstart(i) = st(coord(1))
          lend(i)   = en(coord(1))
          lsize(i)  = sz(coord(1))
          deallocate(st,en,sz)
       elseif (pdim(i) == 3) then    ! distribute across dims(2)
          allocate(st(0:dims(2)-1))
          allocate(en(0:dims(2)-1))
          allocate(sz(0:dims(2)-1))
          call distribute(gsize,dims(2),st,en,sz)
          lstart(i) = st(coord(2))
          lend(i)   = en(coord(2))
          lsize(i)  = sz(coord(2))
          deallocate(st,en,sz)
       end if

    end do
    return   

  end subroutine partition_impl

  !----------------------------------------------------------------------
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
  !---------------------------------------------------------------------- 
  subroutine distribute(data1,proc,st,en,sz)

    implicit none
    ! data1 -- data size in any dimension to be partitioned
    ! proc  -- number of processors in that dimension
    ! st    -- array of starting index
    ! en    -- array of ending index
    ! sz    -- array of local size  (redundent)
    integer data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size1,nl,nu

    size1=data1/proc
    nu = data1 - size1 * proc
    nl = proc - nu
    st(0) = 1
    sz(0) = size1
    en(0) = size1
    do i=1,nl-1
       st(i) = st(i-1) + size1
       sz(i) = size1
       en(i) = en(i-1) + size1
    end do
    size1 = size1 + 1
    do i=nl,proc-1
       st(i) = en(i-1) + 1
       sz(i) = size1
       en(i) = en(i-1) + size1
    end do
    en(proc-1)= data1 
    sz(proc-1)= data1-st(proc-1)+1

    return
  end subroutine distribute

end submodule decomp_2d_init_fin

