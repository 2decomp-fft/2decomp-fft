!! SPDX-License-Identifier: BSD-3-Clause
program io_bench

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_io_adios
   use decomp_2d_io_object_adios
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI
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

   real(mytype), allocatable, dimension(:, :, :) :: u1
   type(d2d_io_adios), save :: io

   double precision :: t1, t2
   integer :: ierror

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   call decomp_2d_testing_init(p_row, p_col, nx, ny, nz)

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   call decomp_2d_io_init()
   call decomp_2d_register_var("io.dat", 1, real_type)

   call alloc_x(u1, .true.)
   call random_number(u1)

   t1 = MPI_WTIME()
   call decomp_2d_adios_write_var(io, u1, 'io.dat')
   call io%end_close
   t2 = MPI_WTIME()

   if (nrank == 0) then
      write (*, *) 'I/O time: ', t2 - t1
      write (*, *) '   '
      write (*, *) 'IO_bench completed '
      write (*, *) '   '
   end if

   deallocate (u1)
   call decomp_2d_io_fin
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_bench

