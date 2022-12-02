!!!! utilities.f90
!!!  Utilities for the initialisation tests.
!!!=====================================================

module init_utils

  use MPI
  use decomp_2d
  
  implicit none

contains

  subroutine check_axis(axis, nexpect)

    character(len=*), intent(in) :: axis
    integer, intent(in) :: nexpect

    integer :: suml
    integer :: sumg
    integer, dimension(3) :: sizes

    integer :: ierr

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

end module init_utils
