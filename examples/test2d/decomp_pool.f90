!! SPDX-License-Identifier: BSD-3-Clause
program prog_decomp_pool

   use iso_fortran_env, only: output_unit, error_unit
   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_testing
   use m_decomp_pool, only: decomp_pool, decomp_pool_get, decomp_pool_free, decomp_pool_nblk

   implicit none

   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: iounit, ierror
   real(mytype), pointer, contiguous, dimension(:, :, :) :: a, b, c, d
#ifdef _GPU
   attributes(device) :: a, b, c, d
#endif

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
   ! The external code can change the number of blocks in the memory pool before calling decomp_2d_init
   !
   !    The default value is 2
   !    transpose operations often require a buffer to pack and another to unpack
   !
   decomp_pool_nblk = 1

   ! Initialize decomp_2d
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   ! Print a warning if there was an error in the inputs
   call decomp_2d_testing_log()

   ! Open the IO units
   iounit = d2d_log_get_unit()
   if (iounit == output_unit .and. nrank > 0) iounit = error_unit

   ! Print the content of the memory pool
   if (iounit /= error_unit) then
      write (iounit, *) ""
      write (iounit, *) "Content of the memory pool at initial stage"
      write (iounit, *) ""
   end if
   call decomp_pool%print(iounit)

   ! Get blocks from the memory pool
   !
   ! 3D arrays with index starting at 1
   call decomp_pool_get(a, xsize)
   call decomp_pool_get(b, ysize)
   call decomp_pool_get(c, zsize)
   ! 3D array with index starting at ystart
   call decomp_pool_get(d, ysize, ystart)

   ! Print the content of the memory pool
   if (iounit /= error_unit) then
      write (iounit, *) ""
      write (iounit, *) "Content of the memory pool after decomp_pool_get"
      write (iounit, *) ""
   end if
   call decomp_pool%print(iounit)

   ! Free blocks in any order
   call decomp_pool_free(d)
   call decomp_pool_free(b)
   call decomp_pool_free(c)
   call decomp_pool_free(a)

   ! Print the content of the memory pool
   if (iounit /= error_unit) then
      write (iounit, *) ""
      write (iounit, *) "Content of the memory pool after decomp_pool_free"
      write (iounit, *) ""
   end if
   call decomp_pool%print(iounit)

   ! Close the IO units
   call d2d_log_close_unit(iounit)

   ! Finalize decomp_2d and MPI
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program prog_decomp_pool
