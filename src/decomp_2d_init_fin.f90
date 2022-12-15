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

  !======================================================================
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
  !======================================================================
  subroutine decomp_2d_init_ref(nx, ny, nz, p_row, p_col, periodic_bc, comm)

     use mpi
     use iso_fortran_env, only: output_unit

     implicit none

     integer, intent(IN) :: nx, ny, nz
     integer, intent(INOUT) :: p_row, p_col
     logical, dimension(3), intent(IN), optional :: periodic_bc
     integer, intent(in), optional :: comm

     integer :: errorcode, ierror, row, col, iounit
     logical, dimension(2) :: periodic
#if defined(_GPU) && defined(_NCCL)
     integer :: cuda_stat
     type(ncclResult) :: nccl_stat
#endif

#ifdef PROFILER
     ! Prepare the profiler if it was not already prepared
     if (decomp_profiler == decomp_profiler_none) call decomp_profiler_prep()
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
     end if

     ! Safety check
     if (MPI_SUCCESS /= 0) call decomp_2d_abort(__FILE__, &
                                                __LINE__, &
                                                MPI_SUCCESS, &
                                                "MPI error check is broken")

     ! If the external code has not set nrank and nproc
     if (nrank == -1) then
        call MPI_COMM_RANK(decomp_2d_comm, nrank, ierror)
        if (ierror /= 0) call decomp_2d_abort(__FILE__, &
                                              __LINE__, &
                                              ierror, &
                                              "MPI_COMM_RANK")
     end if
     if (nproc == -1) then
        call MPI_COMM_SIZE(decomp_2d_comm, nproc, ierror)
        if (ierror /= 0) call decomp_2d_abort(__FILE__, &
                                              __LINE__, &
                                              ierror, &
                                              "MPI_COMM_SIZE")
     end if
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

     if (p_row == 0 .and. p_col == 0) then
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
     call MPI_CART_CREATE(decomp_2d_comm, 2, dims, periodic, &
                          .false., &  ! do not reorder rank
                          DECOMP_2D_COMM_CART_X, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
     periodic(1) = periodic_x
     periodic(2) = periodic_z
     call MPI_CART_CREATE(decomp_2d_comm, 2, dims, periodic, &
                          .false., DECOMP_2D_COMM_CART_Y, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
     periodic(1) = periodic_x
     periodic(2) = periodic_y
     call MPI_CART_CREATE(decomp_2d_comm, 2, dims, periodic, &
                          .false., DECOMP_2D_COMM_CART_Z, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")

     call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X, nrank, 2, coord, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_COORDS")

     ! derive communicators defining sub-groups for ALLTOALL(V)
     call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), &
                       DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SUB")
     call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), &
                       DECOMP_2D_COMM_ROW, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SUB")

     ! gather information for halo-cell support code
     call init_neighbour

     ! actually generate all 2D decomposition information
     call decomp_info_init(nx, ny, nz, decomp_main)

     ! make a copy of the decomposition information associated with the
     ! default global size in these global variables so applications can
     ! use them to create data structures
     xstart = decomp_main%xst
     ystart = decomp_main%yst
     zstart = decomp_main%zst
     xend = decomp_main%xen
     yend = decomp_main%yen
     zend = decomp_main%zen
     xsize = decomp_main%xsz
     ysize = decomp_main%ysz
     zsize = decomp_main%zsz

     ! determine the number of bytes per float number
     ! do not use 'mytype' which is compiler dependent
     ! also possible to use inquire(iolength=...)
     call MPI_TYPE_SIZE(real_type, mytype_bytes, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")

#ifdef EVEN
     if (nrank == 0) write (*, *) 'Padded ALLTOALL optimisation on'
#endif

#if defined(_GPU)
#if defined(_NCCL)
     call MPI_COMM_RANK(DECOMP_2D_COMM_COL, col_rank, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
     call MPI_COMM_RANK(DECOMP_2D_COMM_ROW, row_rank, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
     call MPI_COMM_SIZE(DECOMP_2D_COMM_COL, col_comm_size, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
     call MPI_COMM_SIZE(DECOMP_2D_COMM_ROW, row_comm_size, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")

     allocate (local_to_global_col(col_comm_size), local_to_global_row(row_comm_size))

     local_to_global_col(:) = 0
     local_to_global_row(:) = 0
     local_to_global_col(col_rank + 1) = nrank
     local_to_global_row(row_rank + 1) = nrank

     call mpi_allreduce(MPI_IN_PLACE, local_to_global_col, col_comm_size, MPI_INTEGER, MPI_SUM, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")
     call mpi_allreduce(MPI_IN_PLACE, local_to_global_row, row_comm_size, MPI_INTEGER, MPI_SUM, DECOMP_2D_COMM_ROW, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")

     if (nrank == 0) then
        nccl_stat = ncclGetUniqueId(nccl_uid_2decomp)
     end if
     call MPI_Bcast(nccl_uid_2decomp, int(sizeof(ncclUniqueId)), MPI_BYTE, 0, decomp_2d_comm, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     nccl_stat = ncclCommInitRank(nccl_comm_2decomp, nproc, nccl_uid_2decomp, nrank)
     cuda_stat = cudaStreamCreate(cuda_stream_2decomp)
#endif
#endif

     !
     ! Get the IO unit for decomp_2d setup
     !
     iounit = d2d_listing_get_unit()
     !
     ! Print the decomp_2d setup
     !
     call d2d_listing(iounit)

#ifdef PROFILER
     ! Stop the timer for decomp_2d_init
     if (decomp_profiler_d2d) call decomp_profiler_end("decomp_2d_init")
#endif

     return
  end subroutine decomp_2d_init_ref

  !======================================================================
  ! Routine to be called by applications to clean things up
  !======================================================================
  subroutine decomp_2d_finalize_ref

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
     deallocate (work1_r, work2_r, work1_c, work2_c)
#if defined(_GPU)
     deallocate (work1_r_d, work2_r_d, work1_c_d, work2_c_d)

#if defined(_NCCL)
     nccl_stat = ncclCommDestroy(nccl_comm_2decomp)
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
     end if

     ! Conversion to integer if possible
     read (val, '(i4)', iostat=ierror) decomp_debug
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

     if (nrank == 0) write (*, *) 'In auto-tuning mode......'

     i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors
     allocate (factors(i))
     call findfactor(iproc, factors, nfact)
     if (nrank == 0) write (*, *) 'factors: ', (factors(i), i=1, nfact)

     i_best = nfact/2 + 1
     col = factors(i_best)

     best_p_col = col
     best_p_row = iproc/col
     if (nrank == 0) print *, 'p_row x p_col', best_p_row, best_p_col
     if ((best_p_col == 1) .and. (nrank == 0)) then
        print *, 'WARNING: current 2D DECOMP set-up might not work'
     end if

     deallocate (factors)

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
     neighbour(1, 1) = MPI_PROC_NULL               ! east
     neighbour(1, 2) = MPI_PROC_NULL               ! west
     call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
                         neighbour(1, 4), neighbour(1, 3), ierror) ! north & south
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
     call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
                         neighbour(1, 6), neighbour(1, 5), ierror) ! top & bottom
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

     ! For Y-pencil
     call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
                         neighbour(2, 2), neighbour(2, 1), ierror) ! east & west
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
     neighbour(2, 3) = MPI_PROC_NULL               ! north
     neighbour(2, 4) = MPI_PROC_NULL               ! south
     call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
                         neighbour(2, 6), neighbour(2, 5), ierror) ! top & bottom
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

     ! For Z-pencil
     call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, &
                         neighbour(3, 2), neighbour(3, 1), ierror) ! east & west
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
     call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, &
                         neighbour(3, 4), neighbour(3, 3), ierror) ! north & south
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
     neighbour(3, 5) = MPI_PROC_NULL               ! top
     neighbour(3, 6) = MPI_PROC_NULL               ! bottom
     return
  end subroutine init_neighbour
