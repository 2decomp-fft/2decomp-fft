program timing2d_real

   use mpi
   use decomp_2d
   use decomp_2d_constants
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

   real(mytype), allocatable, dimension(:, :, :) :: data1

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

   real(mytype) :: t1, t2, t3, t4, t5, t6, t7, t8
   integer :: iter, niter =10

   ! Init
   error_flag = .false.
   call MPI_INIT(ierror)
   if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_INIT")
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

   !! ***** global data *****
   !allocate(data1(nx,ny,nz))
   !m = 1
   !do k = 1, nz
   !   do j = 1, ny
   !      do i = 1, nx
   !         data1(i, j, k) = float(m)
   !         m = m + 1
   !      end do
   !   end do
   !end do

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
   ! original x-pensil based data
   !$acc parallel loop default(present) private(m) 
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            m = real(i+(j-1)*nx+(k-1)*nx*ny,mytype)
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

   t2=0._mytype
   t4=0._mytype
   t6=0._mytype
   t8=0._mytype
   do iter=1,niter
      !!!!!!!!!!!!!!!!!!!!!!!
      ! x-pensil ==> y-pensil
      t1 = MPI_WTIME()
      call transpose_x_to_y(u1, u2)
      t2 = t2 + MPI_WTIME()-t1

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
      !$acc parallel loop default(present) private(m)
      do k = yst3, yen3
         do j = yst2, yen2
            do i = yst1, yen1
               m = real(i+(j-1)*nx+(k-1)*nx*ny,mytype)
               if (abs(u2(i, j, k) - m) > 0) error_flag = .true.
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
      if (error_flag) call decomp_2d_abort(1, "error swaping x->y")

      !!!!!!!!!!!!!!!!!!!!!!!
      ! y-pensil ==> z-pensil
      t3 = MPI_WTIME()
      call transpose_y_to_z(u2, u3)
      t4 = t4 + MPI_WTIME()-t3

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

      !$acc parallel loop default(present) private(m)
      do k = zst3, zen3
         do j = zst2, zen2
            do i = zst1, zen1
               m = real(i+(j-1)*nx+(k-1)*nx*ny,mytype)
               if (abs(u3(i, j, k) - m) > 0) error_flag = .true.
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
      if (error_flag) call decomp_2d_abort(2, "error swaping y->z")

      !!!!!!!!!!!!!!!!!!!!!!!
      ! z-pensil ==> y-pensil
      t5 = MPI_WTIME()
      call transpose_z_to_y(u3, u2)
      t6 = t6 + MPI_WTIME()-t5
      ! call decomp_2d_write_one(2,u2,'u2b.dat')

      !$acc parallel loop default(present) private(m)
      do k = yst3, yen3
         do j = yst2, yen2
            do i = yst1, yen1
               m = real(i+(j-1)*nx+(k-1)*nx*ny,mytype)
               if (abs(u2(i, j, k) - m) > 0) error_flag = .true.
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
      if (error_flag) call decomp_2d_abort(3, "error swaping z->y")

      !!!!!!!!!!!!!!!!!!!!!!!
      ! y-pensil ==> x-pensil
      t7 = MPI_WTIME()
      call transpose_y_to_x(u2, u1)
      t8 = t8 + MPI_WTIME()-t7
      ! call decomp_2d_write_one(1,u1,'u1b.dat')

      !$acc parallel loop default(present) private(m)
      do k = xst3, xen3
         do j = xst2, xen2
            do i = xst1, xen1
               m = real(i+(j-1)*nx+(k-1)*nx*ny,mytype)
               if (abs(u1(i, j, k) - m) > 0) error_flag = .true.
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(MPI_IN_PLACE, error_flag, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
      if (error_flag) call decomp_2d_abort(4, "error swaping y->x")
   enddo

   call MPI_ALLREDUCE(t2, t1, 1, real_type, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t1 = t1/real(nproc, mytype)
   call MPI_ALLREDUCE(t4, t3, 1, real_type, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t3 = t3/real(nproc, mytype)
   call MPI_ALLREDUCE(t6, t5, 1, real_type, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t5 = t5/real(nproc, mytype)
   call MPI_ALLREDUCE(t8, t7, 1, real_type, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t7 = t7/real(nproc, mytype)
   t8 = t1+t3+t5+t7
   if (nrank == 0) then
     write(*,*) 'Time X->Y ', t1
     write(*,*) 'Time Y->Z ', t3
     write(*,*) 'Time Z->Y ', t5
     write(*,*) 'Time Y->X ', t7
     write(*,*) 'Time TOT  ', t8
   endif

   if (nrank == 0) then
      write (*, *) " "
      write (*, *) "Real transpose completed"
      write (*, *) " "
   end if

   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)
   !$acc end data
   deallocate (u1, u2, u3)

end program timing2d_real

