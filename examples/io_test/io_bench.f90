program io_bench

   use decomp_2d
   use decomp_2d_io
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
   integer :: nargin, arg, FNLength, status, DecInd
   character(len=80) :: InputFN

   real(mytype), allocatable, dimension(:, :, :) :: u1

   double precision :: t1, t2
   integer :: ierror

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

   call alloc_x(u1, .true.)
   call random_number(u1)

   t1 = MPI_WTIME()
   call decomp_2d_write_one(1, u1, 'io.dat')
   t2 = MPI_WTIME()

   if (nrank == 0) write (*, *) 'I/O time: ', t2 - t1

   deallocate (u1)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_bench

