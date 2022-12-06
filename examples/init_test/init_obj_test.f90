!!!! init_obj_test.f90
!!!  Tests initialising the 2decomp&fft library with
!!!  objects.
!!!=====================================================

program init_obj_test

  use MPI
  use decomp_2d

  use init_utils
  
  implicit none

  type compute_grid
     integer, dimension(:), allocatable :: n
  end type compute_grid

  integer, parameter :: nx = 5
  integer, parameter :: ny = 6
  integer, parameter :: nz = 7
  integer, parameter :: nexpect = nx * ny * nz

  integer :: p_row, p_col

  integer :: ierr

  call MPI_Init(ierr)
  
  p_row = 0; p_col = 0
  call run(p_row, p_col)

  call MPI_Finalize(ierr)
  
contains

  subroutine run(p_row, p_col)

    integer, intent(inout) :: p_row, p_col
    type(compute_grid) :: grid

    call make_compute_grid(3, grid)
    call set_compute_grid(1, nx, grid)
    call set_compute_grid(2, ny, grid)
    call set_compute_grid(3, nz, grid)

    call decomp_2d_init_new(grid, p_row, p_col)

    call check_axis("X", nexpect)
    call check_axis("Y", nexpect)
    call check_axis("Z", nexpect)

    call free_compute_grid(grid) ! TODO: Make this type-bound finaliser
    
    call decomp_2d_finalize()
    
  end subroutine run

  subroutine decomp_2d_init_new(grid, p_row, p_col)

    type(compute_grid), intent(in) :: grid
    integer, intent(inout) :: p_row, p_col

    associate(nx => grid%n(1), &
         ny => grid%n(2), &
         nz => grid%n(3))
      call decomp_2d_init(nx, ny, nz, p_row, p_col)
    end associate

  end subroutine decomp_2d_init_new
  
  subroutine make_compute_grid(ndim, grid)

    integer, intent(in) :: ndim
    type(compute_grid), intent(out) :: grid

    allocate(grid%n(ndim))
    grid%n(:) = 0
    
  end subroutine make_compute_grid

  subroutine set_compute_grid(dim, n, grid)

    integer, intent(in) :: dim                ! Which dimension to set
    integer, intent(in) :: n                  ! Grid size
    type(compute_grid), intent(inout) :: grid ! The computational grid

    grid%n(dim) = n

  end subroutine set_compute_grid

  subroutine free_compute_grid(grid)

    type(compute_grid), intent(inout) :: grid

    deallocate(grid%n)
    
  end subroutine free_compute_grid
  
end program init_obj_test
