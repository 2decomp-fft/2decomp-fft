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

  type processor_grid
     integer, dimension(:), allocatable :: p
  end type processor_grid

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
    type(compute_grid) :: cgrid
    type(processor_grid) :: pgrid

    call make_compute_grid(3, cgrid)
    call set_compute_grid(1, nx, cgrid)
    call set_compute_grid(2, ny, cgrid)
    call set_compute_grid(3, nz, cgrid)

    call make_processor_grid(2, pgrid)
    call set_processor_grid(1, p_row, pgrid)
    call set_processor_grid(2, p_col, pgrid)

    call decomp_2d_init_new(cgrid, pgrid)

    call check_axis("X", nexpect)
    call check_axis("Y", nexpect)
    call check_axis("Z", nexpect)

    call free_compute_grid(cgrid)   ! TODO: Make this type-bound finaliser
    call free_processor_grid(pgrid) ! TODO: Make this type-bound finaliser
    
    call decomp_2d_finalize()
    
  end subroutine run

  subroutine decomp_2d_init_new(cgrid, pgrid)

    type(compute_grid), intent(in) :: cgrid
    type(processor_grid), intent(inout) :: pgrid

    associate(nx => cgrid%n(1), &
         ny => cgrid%n(2), &
         nz => cgrid%n(3), &
         p_row => pgrid%p(1), &
         p_col => pgrid%p(2))
      call decomp_2d_init(nx, ny, nz, p_row, p_col)
    end associate

  end subroutine decomp_2d_init_new
  
  subroutine make_compute_grid(ndim, grid)

    integer, intent(in) :: ndim
    type(compute_grid), intent(out) :: grid

    allocate(grid%n(ndim))
    grid%n(:) = 0
    
  end subroutine make_compute_grid

  subroutine make_processor_grid(ndim, grid)

    integer, intent(in) :: ndim
    type(processor_grid), intent(out) :: grid

    allocate(grid%p(ndim))
    grid%p(:) = 0

  end subroutine make_processor_grid

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

  subroutine set_processor_grid(dim, n, grid)

    integer, intent(in) :: dim                  ! Which dimension to set
    integer, intent(in) :: n                    ! Grid size
    type(processor_grid), intent(inout) :: grid ! The computational grid

    grid%p(dim) = n
    
  end subroutine set_processor_grid

  subroutine free_processor_grid(grid)

    type(processor_grid), intent(inout) :: grid

    deallocate(grid%p)
    
  end subroutine free_processor_grid
  
end program init_obj_test
