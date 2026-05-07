!! SPDX-License-Identifier: BSD-3-Clause

  !======================================================================
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !     periodic_bc  - optional, periodicity flag for halo operations
  !     comm         - optional, MPI communicator, default MPI_COMM_WORLD
  !     complex_pool - optional, flag to use a complex memory pool
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
  !======================================================================
  subroutine decomp_2d_init_ref(nx, ny, nz, p_row, p_col, &
                                periodic_bc, &
                                comm, &
                                complex_pool)

     use mpi
     use iso_fortran_env, only: output_unit

     implicit none

     integer, intent(IN) :: nx, ny, nz
     integer, intent(INOUT) :: p_row, p_col
     logical, dimension(3), intent(IN), optional :: periodic_bc
     integer, intent(in), optional :: comm
     logical, intent(in), optional :: complex_pool

     integer :: errorcode, ierror, row, col, iounit
     logical, dimension(2) :: periodic

     ! Prepare the profiler if it was not already prepared
     if (decomp_profiler == decomp_profiler_none) call decomp_profiler_prep()
     ! Start the profiler
     call decomp_profiler_init()
     ! Start the timer for decomp_2d_init
     if (decomp_profiler_d2d) call decomp_profiler_start("decomp_2d_init")

     call decomp_2d_mpi_init(comm)

     if (nx <= 0) call decomp_2d_abort(__FILE__, __LINE__, nx, "Invalid value for nx")
     if (ny <= 0) call decomp_2d_abort(__FILE__, __LINE__, ny, "Invalid value for ny")
     if (nz <= 0) call decomp_2d_abort(__FILE__, __LINE__, nz, "Invalid value for nz")

     ! Check if the memory pool is available
#if defined(_GPU) || defined(_OPENMP_GPU)
     use_pool = .false.
#else
     use_pool = .true.
#endif

#ifdef DEBUG
     ! Check if a modification of the debug level is needed
     decomp_debug = d2d_get_env_var("DECOMP_2D_DEBUG", decomp_debug)
#endif

     ! Set the default node repartition if needed
     if (decomp_partition_default(1) == DECOMP_PARTITION_UNDEF) then
        decomp_partition_default(1) = d2d_get_env_var("DECOMP_PARTITION_ROW", DECOMP_PARTITION_LAST)
     end if
     if (decomp_partition_default(2) == DECOMP_PARTITION_UNDEF) then
        decomp_partition_default(2) = d2d_get_env_var("DECOMP_PARTITION_COL", DECOMP_PARTITION_LAST)
     end if

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

     d2d_t_y2z_pack = 0.d0
     d2d_t_y2z_mpi = 0.d0
     d2d_t_y2z_unpack = 0.d0
     d2d_t_z2y_pack = 0.d0
     d2d_t_z2y_mpi = 0.d0
     d2d_t_z2y_unpack = 0.d0
     d2d_n_y2z_mpi_calls = 0
     d2d_n_z2y_mpi_calls = 0
     d2d_n_y2z_use_device_ptr = 0
     d2d_n_z2y_use_device_ptr = 0

     if (p_row <= 0 .or. p_col <= 0) then
        ! determine the best 2D processor grid
        call best_2d_grid(nproc, row, col)
        p_row = row
        p_col = col
     else
        if (nproc /= p_row * p_col) then
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
     call decomp_info_init(nx, ny, nz, decomp_main, decomp_partition_default)

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

#if defined(_GPU) && defined(_NCCL)
     call decomp_2d_nccl_init(DECOMP_2D_COMM_COL, DECOMP_2D_COMM_ROW)
#endif

     !
     ! Extend the main memory pool to store complex numbers
     !
     if (use_pool .and. present(complex_pool)) then
        if (complex_pool) then
           decomp_pool_default_type = complex_type
           call mem_pool_set_default_type(complex_type)
           call decomp_pool%new_shape(complex_type, decomp_main)
        end if
     end if

     !
     ! Get the IO unit for decomp_2d setup
     !
     iounit = d2d_log_get_unit()
     !
     ! Print the decomp_2d setup
     !
     call d2d_log(iounit)

     ! Stop the timer for decomp_2d_init
     if (decomp_profiler_d2d) call decomp_profiler_end("decomp_2d_init")

     return
  end subroutine decomp_2d_init_ref

  !======================================================================
  ! Routine to be called by applications to clean things up
  !======================================================================
  subroutine decomp_2d_finalize_ref

     implicit none

     if (decomp_profiler_d2d) call decomp_profiler_start("decomp_2d_fin")

     if (use_pool) call decomp_pool_fin()

     call decomp_2d_mpi_comm_free(DECOMP_2D_COMM_ROW)
     call decomp_2d_mpi_comm_free(DECOMP_2D_COMM_COL)
     call decomp_2d_mpi_comm_free(DECOMP_2D_COMM_CART_X)
     call decomp_2d_mpi_comm_free(DECOMP_2D_COMM_CART_Y)
     call decomp_2d_mpi_comm_free(DECOMP_2D_COMM_CART_Z)

     decomp_partition_default = DECOMP_PARTITION_UNDEF

     call decomp_info_finalize(decomp_main)

     decomp_buf_size = 0
#if defined(_GPU)
     if (associated(work1_r)) nullify (work1_r)
     if (associated(work2_r)) nullify (work2_r)
     if (associated(work1_c)) nullify (work1_c)
     if (associated(work2_c)) nullify (work2_c)
     if (allocated(work1)) deallocate (work1)
     if (allocated(work2)) deallocate (work2)
     call decomp_2d_cumpi_fin()
#if defined(_NCCL)
     call decomp_2d_nccl_fin()
#endif
#elif defined(_OPENMP_GPU)
     if (associated(work1_r)) nullify (work1_r)
     if (associated(work2_r)) nullify (work2_r)
     if (associated(work1_c)) nullify (work1_c)
     if (associated(work2_c)) nullify (work2_c)
     if (work_omp_mapped) then
        !$omp target exit data map(delete:work1(1:size(work1)), work2(1:size(work2)))
        work_omp_mapped = .false.
     end if
     if (allocated(work1)) deallocate (work1)
     if (allocated(work2)) deallocate (work2)
#endif

     call decomp_2d_mpi_fin()

     if (decomp_profiler_d2d) call decomp_profiler_end("decomp_2d_fin")
     ! Finalize the profiler
     call decomp_profiler_fin()

     return
  end subroutine decomp_2d_finalize_ref

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

     i_best = nfact / 2 + 1
     col = factors(i_best)

     best_p_col = col
     best_p_row = iproc / col
     if (nrank == 0) print *, 'p_row x p_col', best_p_row, best_p_col
     if ((best_p_col == 1) .and. (nrank == 0)) then
        print *, 'WARNING: current 2D DECOMP set-up might not work'
     end if

     deallocate (factors)

     return
  end subroutine best_2d_grid
