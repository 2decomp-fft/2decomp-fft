!! SPDX-License-Identifier: BSD-3-Clause

! Testing module for 2decomp&fft library

module decomp_2d_testing

   use decomp_2d_mpi

   implicit none

   logical, save :: use_default = .false.

   public :: decomp_2d_testing_init, &
             decomp_2d_testing_log

   interface decomp_2d_testing_init
      module procedure decomp_2d_testing_init_2int
      module procedure decomp_2d_testing_init_5int
      module procedure decomp_2d_testing_init_6int
   end interface decomp_2d_testing_init

contains

   !
   ! Process command-line arguments, 2 integers
   !
   subroutine decomp_2d_testing_init_2int(p_row, p_col)

      implicit none

      integer, intent(inout) :: p_row, p_col

      integer :: arg, nargin, DecInd, status, FNLength
      character(len=80) :: InputFN

      nargin = command_argument_count()

      if (nargin == 0) return

      if (nargin == 2) then
         do arg = 1, nargin
            call get_command_argument(arg, InputFN, FNLength, status)
            read (InputFN, *, iostat=status) DecInd
            if (arg == 1) then
               p_row = DecInd
            elseif (arg == 2) then
               p_col = DecInd
            end if
         end do
      else
         use_default = .true.
      end if

   end subroutine decomp_2d_testing_init_2int

   !
   ! Process command-line arguments, 5 integers
   !
   subroutine decomp_2d_testing_init_5int(p_row, p_col, nx, ny, nz)

      implicit none

      integer, intent(inout) :: p_row, p_col, nx, ny, nz

      integer :: arg, nargin, DecInd, status, FNLength
      character(len=80) :: InputFN

      nargin = command_argument_count()

      if (nargin == 0) return

      if (nargin == 2 .or. nargin == 5) then
         do arg = 1, nargin
            call get_command_argument(arg, InputFN, FNLength, status)
            read (InputFN, *, iostat=status) DecInd
            if (arg == 1) then
               p_row = DecInd
            elseif (arg == 2) then
               p_col = DecInd
            elseif (arg == 3) then
               nx = DecInd
            elseif (arg == 4) then
               ny = DecInd
            elseif (arg == 5) then
               nz = DecInd
            end if
         end do
      else
         use_default = .true.
      end if

   end subroutine decomp_2d_testing_init_5int

   !
   ! Process command-line arguments, 6 integers
   !
   subroutine decomp_2d_testing_init_6int(p_row, p_col, nx, ny, nz, ntest)

      implicit none

      integer, intent(inout) :: p_row, p_col, nx, ny, nz, ntest

      integer :: arg, nargin, DecInd, status, FNLength
      character(len=80) :: InputFN

      nargin = command_argument_count()

      if (nargin == 0) return

      if (nargin == 2 .or. nargin == 5 .or. nargin == 6) then
         do arg = 1, nargin
            call get_command_argument(arg, InputFN, FNLength, status)
            read (InputFN, *, iostat=status) DecInd
            if (arg == 1) then
               p_row = DecInd
            elseif (arg == 2) then
               p_col = DecInd
            elseif (arg == 3) then
               nx = DecInd
            elseif (arg == 4) then
               ny = DecInd
            elseif (arg == 5) then
               nz = DecInd
            elseif (arg == 6) then
               ntest = DecInd
            end if
         end do
      else
         use_default = .true.
      end if

   end subroutine decomp_2d_testing_init_6int

   !
   ! Print some warning if the input was incorrect
   !
   subroutine decomp_2d_testing_log()

      implicit none

      ! If the command-line was processed correctly, return
      if (.not. use_default) return

      ! Otherwise, print a warning
      call decomp_2d_warning(0, "Invalid command-line input. Using default values")

   end subroutine decomp_2d_testing_log

end module decomp_2d_testing
