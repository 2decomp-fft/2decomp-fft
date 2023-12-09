!! SPDX-License-Identifier: BSD-3-Clause

!
! Dummy submodule when there is no profiler
!
module decomp_2d_profiler

   use decomp_2d_constants, only: decomp_profiler_none

   implicit none

   private

   ! public user routines
   public :: decomp_profiler_init, &
             decomp_profiler_fin, &
             decomp_profiler_prep, &
             decomp_profiler_log, &
             decomp_profiler_start, &
             decomp_profiler_end

   ! Generic interface to initialize the profiler
   interface decomp_profiler_init
      module procedure decomp_profiler_init_noarg
   end interface decomp_profiler_init

   ! Generic interface to finalize the profiler
   interface decomp_profiler_fin
      module procedure decomp_profiler_fin_int
   end interface decomp_profiler_fin

   ! Generic interface for the profiler to log setup
   interface decomp_profiler_log
      module procedure decomp_profiler_log_int
   end interface decomp_profiler_log

   ! Generic interface to prepare the profiler before init.
   interface decomp_profiler_prep
      module procedure decomp_profiler_prep_bool
   end interface decomp_profiler_prep

   ! Generic interface for the profiler to start a given timer
   interface decomp_profiler_start
      module procedure decomp_profiler_start_char
   end interface decomp_profiler_start

   ! Generic interface for the profiler to end a given timer
   interface decomp_profiler_end
      module procedure decomp_profiler_end_char
   end interface decomp_profiler_end

contains

   !
   ! Dummy initialize
   !
   subroutine decomp_profiler_init_noarg()

      implicit none

   end subroutine decomp_profiler_init_noarg

   !
   ! Dummy finalize
   !
   subroutine decomp_profiler_fin_int(decomp_profiler)

      implicit none

      integer, intent(out) :: decomp_profiler

      decomp_profiler = decomp_profiler_none

   end subroutine decomp_profiler_fin_int

   !
   ! Dummy log setup
   !
   subroutine decomp_profiler_log_int(io_unit)

      implicit none

      integer, intent(in) :: io_unit

      write (io_unit, *) "No profiling"

   end subroutine decomp_profiler_log_int

   !
   ! Dummy setup
   !
   subroutine decomp_profiler_prep_bool(decomp_profiler, b_transpose, b_io, b_fft, b_2d, profiler_setup)

      implicit none

      integer, intent(out) :: decomp_profiler
      logical, intent(out) :: b_transpose, b_io, b_fft, b_2d
      logical, dimension(4), intent(in), optional :: profiler_setup

      decomp_profiler = decomp_profiler_none

      associate (tmp => b_transpose); end associate
      associate (tmp => b_io); end associate
      associate (tmp => b_fft); end associate
      associate (tmp => b_2d); end associate
      associate (tmp => profiler_setup); end associate

   end subroutine decomp_profiler_prep_bool

   !
   ! Dummy start a timer
   !
   subroutine decomp_profiler_start_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      associate (tmp => timer_name); end associate

   end subroutine decomp_profiler_start_char

   !
   ! Dummy stop a timer
   !
   subroutine decomp_profiler_end_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      associate (tmp => timer_name); end associate

   end subroutine decomp_profiler_end_char

end module decomp_2d_profiler
