!! SPDX-License-Identifier: BSD-3-Clause
program node_order

   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_testing
#if defined(_GPU)
   use cudafor
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: ierror, io_unit
   type(decomp_info) :: decomp1, decomp2

   ! Init
   call MPI_INIT(ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_INIT")
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   call decomp_2d_testing_init(p_row, p_col, nx, ny, nz)

   !
   !
   ! Before calling decomp_2d_init, the external code can change the variable
   !    decomp_partition_default(2)
   !
   ! The first element is used for the node repartition in row
   ! The second element is used for the node repartition in column
   !
   ! 3D array in X : Y / Z are distributed in row / col
   ! 3D array in Y : X / Z are distributed in row / col
   ! 3D array in Z : X / Y are distributed in row / col
   !
   ! Valid values are :
   !    DECOMP_PARTITION_FIRST => the first CPUs will have extra nodes
   !    DECOMP_PARTITION_LAST => the last CPUs will have the extra nodes
   !
   ! The default value is DECOMP_PARTITION_LAST
   !
   ! In the present example, the external code is not using the default repartition
   !
   decomp_partition_default(:) = DECOMP_PARTITION_FIRST
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   ! Create decomp_info objects with custom repartition
   call decomp_info_init(nx, ny, nz, decomp1, &
                         (/DECOMP_PARTITION_FIRST, DECOMP_PARTITION_LAST/))
   call decomp_info_init(nx, ny, nz, decomp2, &
                         (/DECOMP_PARTITION_LAST, DECOMP_PARTITION_FIRST/))

   ! Print the setup of the decomp_info objects
   io_unit = d2d_log_get_unit()
   if (d2d_log_is_active()) &
      write (io_unit, *) '==========================================================='
   call decomp_info_print(decomp1, io_unit, "first_last")
   if (d2d_log_is_active()) &
      write (io_unit, *) '==========================================================='
   call decomp_info_print(decomp2, io_unit, "last_first")
   if (d2d_log_is_active()) &
      write (io_unit, *) '==========================================================='
   call d2d_log_close_unit(io_unit)

   call decomp_info_finalize(decomp1)
   call decomp_info_finalize(decomp2)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program node_order
