!!!! init_obj_test.f90
!!!  Tests initialising the 2decomp&fft library with
!!!  objects.
!!!=====================================================

program init_obj_test

  use MPI
  use decomp_2d

  use init_utils
  
  implicit none

  type compute_grid_extents
     integer, dimension(:), allocatable :: n
  end type compute_grid_extents

  type processor_grid_extents
     integer, dimension(:), allocatable :: p
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

    call make_compute_grid_extents(3, cgrid_extents)
    call set_compute_grid_extents(1, nx, cgrid_extents)
    call set_compute_grid_extents(2, ny, cgrid_extents)
    call set_compute_grid_extents(3, nz, cgrid_extents)

    call make_processor_grid_extents(2, pgrid_extents)
    call set_processor_grid_extents(1, p_row, pgrid_extents)
    call set_processor_grid_extents(2, p_col, pgrid_extents)

    call decomp_2d_init_new(cgrid_extents, pgrid_extents)

    call check_axis("X", nexpect)
    call check_axis("Y", nexpect)
    call check_axis("Z", nexpect)

    call free_compute_grid_extents(cgrid_extents)   ! TODO: Make this type-bound finaliser
    call free_processor_grid_extents(pgrid_extents) ! TODO: Make this type-bound finaliser
    
    call decomp_2d_finalize()
    
  end subroutine run

  subroutine decomp_2d_init_new(cgrid_extents, pgrid_extents)

    type(compute_grid_extents), intent(in) :: cgrid_extents
    type(processor_grid_extents), intent(inout) :: pgrid_extents

    associate(nx => cgrid_extents%n(1), &
         ny => cgrid_extents%n(2), &
         nz => cgrid_extents%n(3), &
         p_row => pgrid_extents%p(1), &
         p_col => pgrid_extents%p(2))
      call decomp_2d_init(nx, ny, nz, p_row, p_col)
    end associate

  end subroutine decomp_2d_init_new
  
  subroutine make_compute_grid_extents(ndim, grid_extents)

    integer, intent(in) :: ndim
    type(compute_grid_extents), intent(out) :: grid_extents

    allocate(grid_extents%n(ndim))
    grid_extents%n(:) = 0
    
  end subroutine make_compute_grid_extents

  subroutine make_processor_grid_extents(ndim, grid_extents)

    integer, intent(in) :: ndim
    type(processor_grid_extents), intent(out) :: grid_extents

    allocate(grid_extents%p(ndim))
    grid_extents%p(:) = 0

  end subroutine make_processor_grid_extents

  subroutine set_compute_grid_extents(dim, n, grid_extents)

    integer, intent(in) :: dim                                ! Which dimension to set
    integer, intent(in) :: n                                  ! Grid_Extents size
    type(compute_grid_extents), intent(inout) :: grid_extents ! The computational grid_extents

    grid_extents%n(dim) = n

  end subroutine set_compute_grid_extents

  subroutine free_compute_grid_extents(grid_extents)

    type(compute_grid_extents), intent(inout) :: grid_extents

    deallocate(grid_extents%n)
    
  end subroutine free_compute_grid_extents

  subroutine set_processor_grid_extents(dim, n, grid_extents)

    integer, intent(in) :: dim                                  ! Which dimension to set
    integer, intent(in) :: n                                    ! Grid_Extents size
    type(processor_grid_extents), intent(inout) :: grid_extents ! The processor grid_extents

    grid_extents%p(dim) = n
    
  end subroutine set_processor_grid_extents

  subroutine free_processor_grid_extents(grid_extents)

    type(processor_grid_extents), intent(inout) :: grid_extents

    deallocate(grid_extents%p)
    
  end subroutine free_processor_grid_extents
  
end program init_obj_test
