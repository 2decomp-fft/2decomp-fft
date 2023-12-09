!! SPDX-License-Identifier: BSD-3-Clause

!
! Dummy submodule when there is no profiler
!
module decomp_2d_profiler

   use decomp_2d_constants, only: decomp_profiler_none

   implicit none

   !
   ! Integer to select the profiling tool
   !    0 => no profiling, default
   !    1 => Caliper (https://github.com/LLNL/Caliper)
   !
   integer, save, public :: decomp_profiler = decomp_profiler_none
   ! Default : profile everything
   logical, save, public :: decomp_profiler_transpose = .true.
   logical, save, public :: decomp_profiler_io = .true.
   logical, save, public :: decomp_profiler_fft = .true.
   logical, save, public :: decomp_profiler_d2d = .true.

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
      module procedure decomp_profiler_fin_noarg
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
   subroutine decomp_profiler_fin_noarg()

      implicit none

      decomp_profiler = decomp_profiler_none

   end subroutine decomp_profiler_fin_noarg

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
   subroutine decomp_profiler_prep_bool(profiler_setup)

      implicit none

      logical, dimension(4), intent(in), optional :: profiler_setup

      decomp_profiler = decomp_profiler_none

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
