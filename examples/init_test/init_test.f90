!!!! init_test.f90
!!!  Tests initialising the 2decomp&fft library.
!!!=====================================================

program init_test

  use MPI
  use decomp_2d

  use init_utils
  
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

    call check_axis("X", nexpect)
    call check_axis("Y", nexpect)
    call check_axis("Z", nexpect)
    
    call decomp_2d_finalize()
    
  end subroutine run
  
end program init_test
