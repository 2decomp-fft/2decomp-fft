!! SPDX-License-Identifier: BSD-3-Clause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This tests that the counts, lengths and strides of the halo exchange buffers are computed
! correctly.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_exchange_dims

   use m_info, only: decomp_info
   use m_halo, only: halo_extents_t, init_halo_extents

   implicit none

   type(decomp_info) :: decomp ! We need an info object, but it's not actually used unless doing global indices
   logical, parameter :: global = .false. ! Use local indexing

   logical :: passing

   passing = .true.

   print *, "Testing the halo/exchange buffer sizes"

   call test_x()
   call test_y()
   call test_z()

   if (.not. passing) then
      error stop "Test failed"
   else
      print *, "Pass"
   end if

contains

   subroutine test_x()

      integer, parameter :: ipencil = 1
      type(halo_extents_t) :: halo_extents

      integer, dimension(3) :: halos ! The X/Y/Z halo depths
      integer, dimension(3) :: dims  ! The grid dimensions

      integer :: count, expect_count   ! Number of buffers (and expected value)
      integer :: length, expect_length ! Length of buffers (and expected value)
      integer :: stride, expect_stride ! Stride between buffers (and expected value)

      logical :: x_passing

      print *, "+ Testing X-pencil buffers"
      x_passing = .true.

    !! Test North/South buffer
      print *, "-- North/South buffers"

      ! zero halo in pencil axis
      halos = [0, 3, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [0, 3, 2] halos
      expect_count = 12
      expect_length = 33
      expect_stride = 528

      count = halo_extents%buffer_count(2)
      length = halo_extents%buffer_length(2)
      stride = halo_extents%buffer_stride(2)
      x_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims)

      ! halo in pencil axis
      halos = [5, 3, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [5, 3, 2] halos
      expect_count = 12
      expect_length = 63
      expect_stride = 1008

      count = halo_extents%buffer_count(2)
      length = halo_extents%buffer_length(2)
      stride = halo_extents%buffer_stride(2)
      x_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. x_passing

    !! Test Top/Bottom buffer
      print *, "-- Top/Bottom buffers"

      ! zero halo in pencil axis
      halos = [0, 3, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [0, 3, 2] halos
      ! buffer is contiguous in memory so count is 1 and stride is undefined
      expect_count = 1
      expect_length = 1056
      expect_stride = -1

      count = halo_extents%buffer_count(3)
      length = halo_extents%buffer_length(3)
      stride = halo_extents%buffer_stride(3)
      x_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. x_passing

      ! halo in pencil axis
      halos = [5, 3, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [5, 3, 2] halos
      ! buffer is contiguous in memory so count is 1 and stride is undefined
      expect_count = 1
      expect_length = 2016
      expect_stride = -1

      count = halo_extents%buffer_count(3)
      length = halo_extents%buffer_length(3)
      stride = halo_extents%buffer_stride(3)
      x_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. x_passing

      if (x_passing) then
         print *, "-- PASS"
      else
         passing = .false.
      end if

   end subroutine test_x

   subroutine test_y()

      integer, parameter :: ipencil = 2
      type(halo_extents_t) :: halo_extents

      integer, dimension(3) :: halos ! The X/Y/Z halo depths
      integer, dimension(3) :: dims  ! The grid dimensions

      integer :: count, expect_count   ! Number of buffers (and expected value)
      integer :: length, expect_length ! Length of buffers (and expected value)
      integer :: stride, expect_stride ! Stride between buffers (and expected value)

      logical :: y_passing

      print *, "+ Testing Y-pencil buffers"
      y_passing = .true.

    !! Test East/West buffer
      print *, "-- East/West buffers"

      ! zero halo in pencil axis
      halos = [3, 0, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [3, 0, 2] halos
      expect_count = 504
      expect_length = 3
      expect_stride = 17

      count = halo_extents%buffer_count(1)
      length = halo_extents%buffer_length(1)
      stride = halo_extents%buffer_stride(1)
      y_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims)

      ! halo in pencil axis
      halos = [3, 5, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [3, 5, 2] halos
      expect_count = 624
      expect_length = 3
      expect_stride = 17

      count = halo_extents%buffer_count(1)
      length = halo_extents%buffer_length(1)
      stride = halo_extents%buffer_stride(1)
      y_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. y_passing

    !! Test Top/Bottom buffer
      print *, "-- Top/Bottom buffers"

      ! zero halo in pencil axis
      halos = [3, 0, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [3, 0, 2] halos
      ! buffer is contiguous in memory so count is 1 and stride is undefined
      expect_count = 1
      expect_length = 1428
      expect_stride = -1

      count = halo_extents%buffer_count(3)
      length = halo_extents%buffer_length(3)
      stride = halo_extents%buffer_stride(3)
      y_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. y_passing

      ! halo in pencil axis
      halos = [3, 5, 2]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [3, 5, 2] halos
      ! buffer is contiguous in memory so count is 1 and stride is undefined
      expect_count = 1
      expect_length = 1768
      expect_stride = -1

      count = halo_extents%buffer_count(3)
      length = halo_extents%buffer_length(3)
      stride = halo_extents%buffer_stride(3)
      y_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. y_passing

      if (y_passing) then
         print *, "-- PASS"
      else
         passing = .false.
      end if

   end subroutine test_y

   subroutine test_z()

      integer, parameter :: ipencil = 3
      type(halo_extents_t) :: halo_extents

      integer, dimension(3) :: halos ! The X/Y/Z halo depths
      integer, dimension(3) :: dims  ! The grid dimensions

      integer :: count, expect_count   ! Number of buffers (and expected value)
      integer :: length, expect_length ! Length of buffers (and expected value)
      integer :: stride, expect_stride ! Stride between buffers (and expected value)

      logical :: z_passing

      print *, "+ Testing Z-pencil buffers"
      z_passing = .true.

    !! Test East/West buffer
      print *, "-- East/West buffers"

      ! zero halo in pencil axis
      halos = [2, 3, 0]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [2, 3, 0] halos
      expect_count = 384
      expect_length = 2
      expect_stride = 15

      count = halo_extents%buffer_count(1)
      length = halo_extents%buffer_length(1)
      stride = halo_extents%buffer_stride(1)
      z_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims)

      ! halo in pencil axis
      halos = [2, 3, 5]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [2, 3, 5] halos
      expect_count = 864
      expect_length = 2
      expect_stride = 15

      count = halo_extents%buffer_count(1)
      length = halo_extents%buffer_length(1)
      stride = halo_extents%buffer_stride(1)
      z_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. z_passing

    !! Test North/South buffer
      print *, "-- North/South buffers"

      ! zero halo in pencil axis
      halos = [2, 3, 0]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [2, 3, 0] halos
      expect_count = 8
      expect_length = 45
      expect_stride = 720

      count = halo_extents%buffer_count(2)
      length = halo_extents%buffer_length(2)
      stride = halo_extents%buffer_stride(2)
      z_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. z_passing

      ! halo in pencil axis
      halos = [2, 3, 5]
      dims = [11, 42, 8]
      halo_extents = init_halo_extents(ipencil, dims, decomp, halos, global)

      ! Expected values based on [11, 42, 8] local grid with [2, 3, 5] halos
      expect_count = 18
      expect_length = 45
      expect_stride = 720

      count = halo_extents%buffer_count(2)
      length = halo_extents%buffer_length(2)
      stride = halo_extents%buffer_stride(2)
      z_passing = check_buffer([count, length, stride], &
                               [expect_count, expect_length, expect_stride], &
                               halos, dims) &
                  .and. z_passing

      if (z_passing) then
         print *, "-- PASS"
      else
         passing = .false.
      end if

   end subroutine test_z

   logical function check_buffer(buffer, expected_buffer, halos, dims) result(check)
      integer, dimension(3), intent(in) :: buffer ! [count, length, stride]
      integer, dimension(3), intent(in) :: expected_buffer
      integer, dimension(3), intent(in) :: halos  ! The halo definition
      integer, dimension(3), intent(in) :: dims   ! The grid definition

      if (any(buffer /= expected_buffer)) then
         print *, "--- FAIL"
         print *, "--- halo definition: ", halos
         print *, "--- grid definition: ", dims
         print *, "--- expected count: ", expected_buffer(1), " got ", buffer(1)
         print *, "--- expected length: ", expected_buffer(2), " got ", buffer(2)
         print *, "--- expected stride: ", expected_buffer(3), " got ", buffer(3)
         check = .false.
      else
         print *, "--- PASS"
         check = .true.
      end if
   end function check_buffer

end program test_exchange_dims
