!! SPDX-License-Identifier: BSD-3-Clause
program test2d

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

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3

   integer :: i, j, k, ierror
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3
   integer :: yst1, yst2, yst3
   integer :: yen1, yen2, yen3
   integer :: zst1, zst2, zst3
   integer :: zen1, zen2, zen3
   logical :: error_flag
   real(mytype) :: m

   ! Init
   error_flag = .false.
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

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   ! Fill the local index
   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)
   yst1 = ystart(1); yen1 = yend(1)
   yst2 = ystart(2); yen2 = yend(2)
   yst3 = ystart(3); yen3 = yend(3)
   zst1 = zstart(1); zen1 = zend(1)
   zst2 = zstart(2); zen2 = zend(2)
   zst3 = zstart(3); zen3 = zend(3)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Testing the swap routines
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call alloc_x(u1, opt_global=.true.)
   call alloc_y(u2, opt_global=.true.)
   call alloc_z(u3, opt_global=.true.)

   !$acc data copy(u1,u2,u3)
   ! original x-pencil based data
   !$acc parallel loop default(present) collapse(3) private(m)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            m = real(i + (j - 1) * nx + (k - 1) * nx * ny, mytype)
            u1(i, j, k) = m
         end do
      end do
   end do
   !$acc end loop

10 format(15I5)

#ifdef DEBUG
   if (nrank == 0) then
      !$acc update self(u1)
      write (*, *) 'Numbers held on Rank 0'
      write (*, *) ' '
      write (*, *) 'X-pencil'
      write (*, 10) int(u1)
   end if
#endif

   ! call decomp_2d_write_one(1,u1,'u1.dat')

  !!!!!!!!!!!!!!!!!!!!!!!
   ! x-pencil ==> y-pencil
   call transpose_x_to_y(u1, u2)

#ifdef DEBUG
   if (nrank == 0) then
      !$acc update self(u2)
      write (*, *) ' '
      write (*, *) 'Y-pencil'
      write (*, 10) int(u2)
   end if
#endif

   ! call decomp_2d_write_one(2,u2,'u2.dat')
   ! 'u1.dat' and 'u2.dat' should be identical byte-by-byte

   ! also check the transposition this way
   !$acc parallel loop default(present) collapse(3) private(m) reduction(.or.:error_flag)
   do k = yst3, yen3
      do j = yst2, yen2
         do i = yst1, yen1
            m = real(i + (j - 1) * nx + (k - 1) * nx * ny, mytype)
            if (abs(u2(i, j, k) - m) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(1, "error swaping x->y")
   if (nrank == 0) write (*, *) 'SWAPPING X->Y OK'

   !!!!!!!!!!!!!!!!!!!!!!!
   ! y-pencil ==> z-pencil
   call transpose_y_to_z(u2, u3)

#ifdef DEBUG
   if (nrank == 0) then
      !$acc update self(u3)
      write (*, *) ' '
      write (*, *) 'Z-pencil'
      write (*, 10) int(u3)
   end if
#endif

   ! call decomp_2d_write_one(3,u3,'u3.dat')
   ! 'u1.dat','u2.dat' and 'u3.dat' should be identical

   !$acc parallel loop default(present) collapse(3) private(m) reduction(.or.:error_flag)
   do k = zst3, zen3
      do j = zst2, zen2
         do i = zst1, zen1
            m = real(i + (j - 1) * nx + (k - 1) * nx * ny, mytype)
            if (abs(u3(i, j, k) - m) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(2, "error swaping y->z")
   if (nrank == 0) write (*, *) 'SWAPPING Y->Z OK'

  !!!!!!!!!!!!!!!!!!!!!!!
   ! z-pencil ==> y-pencil
   call transpose_z_to_y(u3, u2)
   ! call decomp_2d_write_one(2,u2,'u2b.dat')

   !$acc parallel loop default(present) collapse(3) private(m) reduction(.or.:error_flag)
   do k = yst3, yen3
      do j = yst2, yen2
         do i = yst1, yen1
            m = real(i + (j - 1) * nx + (k - 1) * nx * ny, mytype)
            if (abs(u2(i, j, k) - m) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(3, "error swaping z->y")
   if (nrank == 0) write (*, *) 'SWAPPING Z->Y OK'

   !!!!!!!!!!!!!!!!!!!!!!!
   ! y-pencil ==> x-pencil
   call transpose_y_to_x(u2, u1)
   ! call decomp_2d_write_one(1,u1,'u1b.dat')

   !$acc parallel loop default(present) collapse(3) private(m) reduction(.or.:error_flag)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            m = real(i + (j - 1) * nx + (k - 1) * nx * ny, mytype)
            if (abs(u1(i, j, k) - m) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(4, "error swaping y->x")
   if (nrank == 0) write (*, *) 'SWAPPING Y->Z OK'

   if (nrank == 0) then
      write (*, *) " "
      write (*, *) "test2d completed"
      write (*, *) " "
   end if

   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)
   !$acc end data
   deallocate (u1, u2, u3)

end program test2d

