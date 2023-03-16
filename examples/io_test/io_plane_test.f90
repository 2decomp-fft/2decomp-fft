!!
!! FIXME The issue below is specific to GPU and should be discussed in a dedicated github issue
!!
!! NB in case of GPU only the writing in the aligned pencil (i.e. X for a 1 array) is performed.
!! IO subrotines needs update for non managed GPU case
!!
program io_plane_test

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

   real(mytype), allocatable, dimension(:, :, :) :: data1
   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3

   real(mytype), allocatable, dimension(:, :, :) :: work

   integer :: i, j, k, m, ierror, iol
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3
   logical :: found

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
            data1(i, j, k) = real(m, mytype)
            m = m + 1
         end do
      end do
   end do

   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   ! For GPU we port the global data create the different pencil arrays
   ! Move back to host the arrays for writing on disk
   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)

   !$acc data copyin(data1) copy(u1,u2,u3)
   ! original X-pensil based data
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
         end do
      end do
   end do
   !$acc end loop
   call transpose_x_to_y(u1, u2)
   call transpose_y_to_z(u2, u3)
   !$acc update self(u1)
   !$acc update self(u2)
   !$acc update self(u3)
   !$acc end data
   ! X-pencil data
   call decomp_2d_write_plane(1, u1, 1, nx/2, '.', 'x_pencil-x_plane.dat', 'test')
   call decomp_2d_write_plane(1, u1, 2, ny/2, '.', 'x_pencil-y_plane.dat', 'test')
   call decomp_2d_write_plane(1, u1, 3, nz/2, '.', 'x_pencil-z_plane.dat', 'test')
   ! Y-pencil data
   call decomp_2d_write_plane(2, u2, 2, ny/2, '.', 'y_pencil-y_plane.dat', 'test')
   call decomp_2d_write_plane(2, u2, 1, nx/2, '.', 'y_pencil-x_plane.dat', 'test')
   call decomp_2d_write_plane(2, u2, 3, nz/2, '.', 'y_pencil-z_plane.dat', 'test')
   ! Z-pencil data
   call decomp_2d_write_plane(3, u3, 1, nx/2, '.', 'z_pencil-x_plane.dat', 'test')
   call decomp_2d_write_plane(3, u3, 2, ny/2, '.', 'z_pencil-y_plane.dat', 'test')
   call decomp_2d_write_plane(3, u3, 3, nz/2, '.', 'z_pencil-z_plane.dat', 'test')
   ! Attemp to read the files
   if (nrank == 0) then
      inquire (iolength=iol) data1(1, 1, 1)

      ! X-plane
      inquire (file='x_pencil-x_plane.dat', exist=found)
      if (found) then
         allocate (work(1, ny, nz))
         open (10, FILE='x_pencil-x_plane.dat', FORM='unformatted', &
               ACCESS='DIRECT', RECL=iol)
         m = 1
         do k = 1, nz
            do j = 1, ny
               read (10, rec=m) work(1, j, k)
               m = m + 1
            end do
         end do
!        write(*,*) ' '
!        write(*,'(15I5)') int(work)
         close (10)
         deallocate (work)
         write (*, *) 'passed self test x-plane'
      else
         write (*, *) "Warning : x_pencil-x_plane.dat is missing"
      end if

      ! Y-plane
      inquire (file='x_pencil-y_plane.dat', exist=found)
      if (found) then
         allocate (work(nx, 1, nz))
         open (10, FILE='x_pencil-y_plane.dat', FORM='unformatted', &
               ACCESS='DIRECT', RECL=iol)
         m = 1
         do k = 1, nz
            do i = 1, nx
               read (10, rec=m) work(i, 1, k)
               m = m + 1
            end do
         end do
!        write(*,*) ' '
!        write(*,'(15I5)') int(work)
         close (10)
         deallocate (work)
         write (*, *) 'passed self test y-plane'
      else
         write (*, *) 'Warning : x_pencil-y_plane.dat is missing'
      end if

      ! Z-plane
      inquire (file='x_pencil-z_plane.dat', exist=found)
      if (found) then
         allocate (work(nx, ny, 1))
         open (10, FILE='x_pencil-z_plane.dat', FORM='unformatted', &
               ACCESS='DIRECT', RECL=iol)
         m = 1
         do j = 1, ny
            do i = 1, nx
               read (10, rec=m) work(i, j, 1)
               m = m + 1
            end do
         end do
!        write(*,*) ' '
!        write(*,'(15I5)') int(work)
         close (10)
         deallocate (work)
         write (*, *) 'passed self test z-plane'
      else
         write (*, *) 'Warning : x_pencil-z_plane.dat is missing'
      end if

   end if

   deallocate (u1, u2, u3)
   deallocate (data1)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_plane_test
