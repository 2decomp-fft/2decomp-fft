program io_test

   use mpi
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d
   use decomp_2d_io
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
   integer :: nargin, arg, FNLength, status, DecInd
   character(len=80) :: InputFN

#ifdef COMPLEX_TEST
   complex(mytype), allocatable, dimension(:, :, :) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   integer :: i, j, k, m, ierror
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot/4) + 1
   nx = nx_base*resize_domain
   ny = ny_base*resize_domain
   nz = nz_base*resize_domain
   ! Now we can check if user put some inputs
   ! Handle input file like a boss -- GD
   nargin=command_argument_count()
   if ((nargin==0).or.(nargin==2).or.(nargin==5)) then
      do arg = 1, nargin
         call get_command_argument(arg, InputFN, FNLength, status)
         read(InputFN, *, iostat=status) DecInd
         if (arg.eq.1) then
            p_row = DecInd
         elseif (arg.eq.2) then
            p_col = DecInd
         elseif (arg.eq.3) then
            nx = DecInd
         elseif (arg.eq.4) then
            ny = DecInd
         elseif (arg.eq.5) then
            nz = DecInd
         endif
      enddo
   else
      ! nrank not yet computed we need to avoid write
      ! for every rank
      call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
      if (nrank==0) then
         print *, "This Test takes no inputs or 2 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "or 5 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "  3) nx "
         print *, "  4) ny "
         print *, "  5) nz "
         print *, "Number of inputs is not correct and the defult settings"
         print *, "will be used"
      endif
   endif
   
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   ! ***** global data *****
   allocate(data1(nx,ny,nz))
   m = 1
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
#ifdef COMPLEX_TEST
            data1(i, j, k) = cmplx(real(m, mytype), real(nx*ny*nz - m, mytype))
#else
            data1(i, j, k) = real(m, mytype)
#endif
            m = m + 1
         end do
      end do
   end do

   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   call alloc_x(u1b, .true.)
   call alloc_y(u2b, .true.)
   call alloc_z(u3b, .true.)

   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)
   ! original x-pencil based data
   !$acc data copyin(data1) copy(u1,u2,u3)
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
         end do
      end do
   end do
   !$acc end loop

   ! transpose
   call transpose_x_to_y(u1, u2)
   call transpose_y_to_z(u2, u3)
   !$acc update self(u1)
   !$acc update self(u2)
   !$acc update self(u3)
   !$acc end data

   ! write to disk
   call decomp_2d_write_one(1, u1, '.', 'u1.dat', 0, 'test')
   call decomp_2d_write_one(2, u2, '.', 'u2.dat', 0, 'test')
   call decomp_2d_write_one(3, u3, '.', 'u3.dat', 0, 'test')

   ! read back to different arrays
   call decomp_2d_read_one(1, u1b, '.', 'u1.dat', 'test', reduce_prec=.false.)
   call decomp_2d_read_one(2, u2b, '.', 'u2.dat', 'test', reduce_prec=.false.)
   call decomp_2d_read_one(3, u3b, '.', 'u3.dat', 'test', reduce_prec=.false.)

   ! compare
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            if (abs((u1(i, j, k) - u1b(i, j, k))) > eps) stop 1
         end do
      end do
   end do

   do k = ystart(3), yend(3)
      do j = ystart(2), yend(2)
         do i = ystart(1), yend(1)
            if (abs((u2(i, j, k) - u2b(i, j, k))) > eps) stop 2
         end do
      end do
   end do

   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            if (abs((u3(i, j, k) - u3b(i, j, k))) > eps) stop 3
         end do
      end do
   end do

   ! Also check against the global data array
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            if (abs(data1(i, j, k) - u1b(i, j, k)) > eps) stop 4
         end do
      end do
   end do

   do k = ystart(3), yend(3)
      do j = ystart(2), yend(2)
         do i = ystart(1), yend(1)
            if (abs((data1(i, j, k) - u2b(i, j, k))) > eps) stop 5
         end do
      end do
   end do

   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            if (abs((data1(i, j, k) - u3b(i, j, k))) > eps) stop 6
         end do
      end do
   end do

   deallocate (u1, u2, u3)
   deallocate (u1b, u2b, u3b)
   deallocate (data1)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_test
