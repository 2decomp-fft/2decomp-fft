!! SPDX-License-Identifier: BSD-3-Clause
!!!=====================================================
!!!! init_test.f90
!!!  Tests initialising the 2decomp&fft library.
!!!=====================================================

program init_test

   use MPI
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_testing
   use decomp_2d
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
   integer :: nargin
   integer :: nexpect
   integer :: ierror

   call MPI_Init(ierror)

   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   nargin = command_argument_count()
   if ((nargin == 0) .or. (nargin == 2) .or. (nargin == 5)) then
      call decomp_2d_testing_init(p_row, p_col, nx, ny, nz)
   else
      ! nrank not yet computed we need to avoid write
      ! for every rank
      call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
      if (nrank == 0) then
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
      end if
   end if
   nexpect = nx * ny * nz
   call run(p_row, p_col)

   call MPI_Finalize(ierror)

contains

   subroutine run(p_row, p_col)

      integer, intent(inout) :: p_row, p_col

      call decomp_2d_init(nx, ny, nz, p_row, p_col)

      call check_axis("X")
      call check_axis("Y")
      call check_axis("Z")

      call decomp_2d_finalize()

   end subroutine run

   subroutine check_axis(axis)

      character(len=*), intent(in) :: axis

      integer :: suml
      integer :: sumg
      integer, dimension(3) :: sizes

      if (axis == "X") then
         sizes = xsize
      else if (axis == "Y") then
         sizes = ysize
      else if (axis == "Z") then
         sizes = zsize
      else
         sizes = 0
         if (nrank == 0) print *, "ERROR: unknown axis requested!"
         stop 1
      end if

      suml = product(sizes)
      call MPI_Allreduce(suml, sumg, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)

      if (sumg /= nexpect) then
         if (nrank == 0) print *, "ERROR: got ", sumg, " nodes, expected ", nexpect
         stop 1
      else
         if (nrank == 0) print *, "Init Test pass for axis ", axis
      end if

   end subroutine check_axis

end program init_test
