program test2d

   use mpi

   use decomp_2d
   ! use decomp_2d_io

   implicit none

   integer, parameter :: nx = 17, ny = 13, nz = 11
   integer :: p_row = 0, p_col = 0

   real(mytype), dimension(nx, ny, nz) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3

   integer :: i, j, k, m, ierror
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3
   integer :: yst1, yst2, yst3
   integer :: yen1, yen2, yen3
   integer :: zst1, zst2, zst3
   integer :: zen1, zen2, zen3
   logical :: error_flag

   ! Init
   error_flag = .false.
   call MPI_INIT(ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_INIT")
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   ! ***** global data *****
   m = 1
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            data1(i, j, k) = float(m)
            m = m + 1
         end do
      end do
   end do

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

   !$acc data copyin(data1) copy(u1,u2,u3)
   ! original x-pensil based data
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
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
   ! x-pensil ==> y-pensil
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
   !$acc parallel loop default(present)
   do k = yst3, yen3
      do j = yst2, yen2
         do i = yst1, yen1
            if (abs(u2(i, j, k) - data1(i, j, k)) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(1, "error swaping x->y")

  !!!!!!!!!!!!!!!!!!!!!!!
   ! y-pensil ==> z-pensil
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

   !$acc parallel loop default(present)
   do k = zst3, zen3
      do j = zst2, zen2
         do i = zst1, zen1
            if (abs(u3(i, j, k) - data1(i, j, k)) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(2, "error swaping y->z")

  !!!!!!!!!!!!!!!!!!!!!!!
   ! z-pensil ==> y-pensil
   call transpose_z_to_y(u3, u2)
   ! call decomp_2d_write_one(2,u2,'u2b.dat')

   !$acc parallel loop default(present)
   do k = yst3, yen3
      do j = yst2, yen2
         do i = yst1, yen1
            if (abs(u2(i, j, k) - data1(i, j, k)) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(3, "error swaping z->y")

   !!!!!!!!!!!!!!!!!!!!!!!
   ! y-pensil ==> x-pensil
   call transpose_y_to_x(u2, u1)
   ! call decomp_2d_write_one(1,u1,'u1b.dat')

   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            if (abs(u1(i, j, k) - data1(i, j, k)) > 0) error_flag = .true.
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
   if (error_flag) call decomp_2d_abort(4, "error swaping y->x")

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

