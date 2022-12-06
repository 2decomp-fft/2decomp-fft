!!!! init_obj_test.f90
!!!  Tests initialising the 2decomp&fft library with
!!!  objects.
!!!=====================================================

program init_obj_test

  use MPI
  use decomp_2d

  use init_utils
  
  implicit none

  type grid_extents
     integer, dimension(:), allocatable :: g
  end type grid_extents
  
  type, extends(grid_extents) :: compute_grid_extents
  end type compute_grid_extents

  type, extends(grid_extents) :: processor_grid_extents
  end type processor_grid_extents

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
    type(compute_grid_extents) :: cgrid_extents
    type(processor_grid_extents) :: pgrid_extents

    call make_grid_extents(3, cgrid_extents)
    call set_grid_extents(1, nx, cgrid_extents)
    call set_grid_extents(2, ny, cgrid_extents)
    call set_grid_extents(3, nz, cgrid_extents)

    call make_grid_extents(2, pgrid_extents)
    call set_grid_extents(1, p_row, pgrid_extents)
    call set_grid_extents(2, p_col, pgrid_extents)

    call decomp_2d_init_new(cgrid_extents, pgrid_extents)

    call check_axis("X", nexpect)
    call check_axis("Y", nexpect)
    call check_axis("Z", nexpect)

    call free_grid_extents(cgrid_extents) ! TODO: Make this type-bound finaliser
    call free_grid_extents(pgrid_extents) ! TODO: Make this type-bound finaliser
    
    call decomp_2d_finalize()
    
  end subroutine run

  subroutine decomp_2d_init_new(cgrid, pgrid)

    type(compute_grid_extents), intent(in) :: cgrid
    type(processor_grid_extents), intent(inout) :: pgrid

    associate(nx => get_grid_extents(1, cgrid), &
         ny => get_grid_extents(2, cgrid), &
         nz => get_grid_extents(3, cgrid), &
         p_row => pgrid%g(1), &
         p_col => pgrid%g(2))
      call decomp_2d_init(nx, ny, nz, p_row, p_col)
    end associate

  end subroutine decomp_2d_init_new
  
  subroutine make_grid_extents(ndim, grid)

    integer, intent(in) :: ndim
    class(grid_extents), intent(out) :: grid

    allocate(grid%g(ndim))
    grid%g(:) = 0
    
  end subroutine make_grid_extents

  subroutine set_grid_extents(dim, n, grid)

    integer, intent(in) :: dim                 ! Which dimension to set
    integer, intent(in) :: n                   ! Grid_Extents size
    class(grid_extents), intent(inout) :: grid ! The grid_extents

    grid%g(dim) = n

  end subroutine set_grid_extents

  integer function get_grid_extents(dim, grid)

    integer, intent(in) :: dim
    class(grid_extents), intent(in) :: grid

    get_grid_extents = grid%g(dim)

  end function get_grid_extents

  subroutine free_grid_extents(grid)

    class(grid_extents), intent(inout) :: grid

    deallocate(grid%g)
    
  end subroutine free_grid_extents
  
end program init_obj_test
