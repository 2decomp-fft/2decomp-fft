!! SPDX-License-Identifier: BSD-3-Clause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example tests the allocation of arrays with halos, confirming
! that their extents and indexing are as expected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program alloc_halo_test

   use mpi

   use decomp_2d
   use decomp_2d_constants

   use decomp_2d_testing

   implicit none

   integer, parameter :: nx_base = 65, ny_base = 48, nz_base = 33
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer, parameter :: ih = 1, jh = 1, kh = 1 ! Halo depths

   real(mytype), allocatable, dimension(:, :, :) :: u1h, u2h, u3h

   integer :: irank, ierror

   logical :: passing, allpassing

   ! Initialisation
   call MPI_INIT(ierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror)
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

   passing = .true.

   ! Set default halo depths
   decomp_main%xlevel = [ih, jh, kh]
   decomp_main%ylevel = [ih, jh, kh]
   decomp_main%zlevel = [ih, jh, kh]

   ! Create and check arrays

   call alloc_x(u1h)
   if (any(lbound(u1h) /= [1, 1, 1] - [ih, jh, kh])) then
      print *, "Alloc_x: FAIL - expected halo-padded arrays to start from index 0, got ", lbound(u1h)
      passing = .false.
   end if
   if (any(ubound(u1h) /= xsize + [ih, jh, kh])) then
      print *, "Alloc_x: FAIL - expected halo-padded arrays to end at ", xsize + [ih, jh, kh], ", got ", ubound(u1h)
      passing = .false.
   end if
   deallocate (u1h)
   call alloc_x(u1h, opt_global=.true.)
   if (any(lbound(u1h) /= xstart - [ih, jh, kh])) then
      print *, "Alloc_x: FAIL - expected halo-padded arrays to start from index 0, got ", lbound(u1h)
      passing = .false.
   end if
   if (any(ubound(u1h) /= xend + [ih, jh, kh])) then
      print *, "Alloc_x: FAIL - expected halo-padded arrays to end at ", xsize + [ih, jh, kh], ", got ", ubound(u1h)
      passing = .false.
   end if
   deallocate (u1h)

   call alloc_y(u2h)
   if (any(lbound(u2h) /= [1, 1, 1] - [ih, jh, kh])) then
      print *, "Alloc_y: FAIL - expected halo-padded arrays to start from index 0, got ", lbound(u2h)
      passing = .false.
   end if
   if (any(ubound(u2h) /= ysize + [ih, jh, kh])) then
      print *, "Alloc_y: FAIL - expected halo-padded arrays to end at ", ysize + [ih, jh, kh], ", got ", ubound(u2h)
      passing = .false.
   end if
   deallocate (u2h)
   call alloc_y(u2h, opt_global=.true.)
   if (any(lbound(u2h) /= ystart - [ih, jh, kh])) then
      print *, "Alloc_y: FAIL - expected halo-padded arrays to start from index 0, got ", lbound(u2h)
      passing = .false.
   end if
   if (any(ubound(u2h) /= yend + [ih, jh, kh])) then
      print *, "Alloc_y: FAIL - expected halo-padded arrays to end at ", ysize + [ih, jh, kh], ", got ", ubound(u2h)
      passing = .false.
   end if
   deallocate (u2h)

   call alloc_z(u3h)
   if (any(lbound(u3h) /= [1, 1, 1] - [ih, jh, kh])) then
      print *, "Alloc_z: FAIL - expected halo-padded arrays to start from index 0, got ", lbound(u3h)
      passing = .false.
   end if
   if (any(ubound(u3h) /= zsize + [ih, jh, kh])) then
      print *, "Alloc_z: FAIL - expected halo-padded arrays to end at ", zsize + [ih, jh, kh], ", got ", ubound(u3h)
      passing = .false.
   end if
   deallocate (u3h)
   call alloc_z(u3h, opt_global=.true.)
   if (any(lbound(u3h) /= zstart - [ih, jh, kh])) then
      print *, "Alloc_z: FAIL - expected halo-padded arrays to start from index 0, got ", lbound(u3h)
      passing = .false.
   end if
   if (any(ubound(u3h) /= zend + [ih, jh, kh])) then
      print *, "Alloc_z: FAIL - expected halo-padded arrays to end at ", zsize + [ih, jh, kh], ", got ", ubound(u3h)
      passing = .false.
   end if
   deallocate (u3h)

   call MPI_Allreduce(passing, allpassing, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)

   call decomp_2d_finalize()

   ! Report test results
   if (.not. allpassing) then
      call decomp_2d_abort(1, "Error in allocate with halos test")
   end if

   call MPI_Finalize(ierror)

end program alloc_halo_test

