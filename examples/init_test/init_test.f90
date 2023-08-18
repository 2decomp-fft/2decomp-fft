!!!! init_test.f90
!!!  Tests initialising the 2decomp&fft library.
!!!=====================================================

program init_test

   use MPI
   use decomp_2d

   implicit none

   integer, parameter :: nx = 5
   integer, parameter :: ny = 6
   integer, parameter :: nz = 7
   integer, parameter :: nexpect = nx*ny*nz

   integer :: p_row, p_col

   integer :: ierr

   call MPI_Init(ierr)

   p_row = 0; p_col = 0
   call run(p_row, p_col)

   call MPI_Finalize(ierr)

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
         print *, "ERROR: unknown axis requested!"
         stop 1
      end if

      suml = product(sizes)
      call MPI_Allreduce(suml, sumg, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      if (sumg /= nexpect) then
         print *, "ERROR: got ", sumg, " nodes, expected ", nexpect
         stop 1
      end if

   end subroutine check_axis

end program init_test
