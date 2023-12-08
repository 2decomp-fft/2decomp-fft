!! SPDX-License-Identifier: BSD-3-Clause

!
! Dummy submodule when there is no profiler
!
submodule(decomp_2d) d2d_profiler_none

   implicit none

contains

   !
   ! Dummy initialize
   !
   module subroutine decomp_profiler_init_noarg

      implicit none

      decomp_profiler = decomp_profiler_none

   end subroutine decomp_profiler_init_noarg

   !
   ! Dummy finalize
   !
   module subroutine decomp_profiler_fin_noarg

      implicit none

      decomp_profiler = decomp_profiler_none

   end subroutine decomp_profiler_fin_noarg

   !
   ! Dummy log setup
   !
   module subroutine decomp_profiler_log_int(io_unit)

      implicit none

      integer, intent(in) :: io_unit

      write (io_unit, *) "No profiling"

   end subroutine decomp_profiler_log_int

   !
   ! Dummy setup
   !
   module subroutine decomp_profiler_prep_bool(profiler_setup)

      implicit none

      logical, dimension(4), intent(in), optional :: profiler_setup

      associate (tmp => profiler_setup); end associate

   end subroutine decomp_profiler_prep_bool

   !
   ! Dummy start a timer
   !
   module subroutine decomp_profiler_start_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      associate (tmp => timer_name); end associate

   end subroutine decomp_profiler_start_char

   !
   ! Dummy stop a timer
   !
   module subroutine decomp_profiler_end_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      associate (tmp => timer_name); end associate

   end subroutine decomp_profiler_end_char

end submodule d2d_profiler_none
