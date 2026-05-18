!! SPDX-License-Identifier: BSD-3-Clause

! Preprocessor macro to deal with unused variables
#define unused(x) associate(tmp => x); end associate

!
! Dummy module when there is no profiler
!
module decomp_2d_profiler

   use, intrinsic :: iso_fortran_env, only: real64
   use decomp_2d_constants, only: decomp_profiler_none
   use decomp_2d_mpi, only: nrank, nproc, decomp_2d_abort, decomp_2d_mpi_allreduce
   use MPI, only: MPI_WTIME, MPI_MAX, MPI_MIN, MPI_SUM

   implicit none

   !
   ! Integer to select the profiling tool
   !    0 => no profiling, default
   !    1 => Caliper (https://github.com/LLNL/Caliper)
   !
   integer, save, public :: decomp_profiler = decomp_profiler_none
   ! Default : profile everything
   logical, parameter :: default_profiler = .false.
   logical, save, public :: decomp_profiler_transpose = default_profiler
   logical, save, public :: decomp_profiler_io = default_profiler
   logical, save, public :: decomp_profiler_fft = default_profiler
   logical, save, public :: decomp_profiler_d2d = default_profiler

   ! Number of timers
   integer, save :: ncur_timers = 0
   integer, save :: nmax_timers = 0
   ! Duration and starting time
   double precision, save, allocatable, dimension(:) :: timer, timer_start
   ! Number of calls
   integer, save, allocatable, dimension(:) :: timer_n
   ! Name of the timers
   character(len=128), save, allocatable, dimension(:) :: timer_name

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
   ! Initialize the basic timer module
   !
   subroutine decomp_profiler_init_noarg()

      implicit none

      ! Estimate the number of timers
      if (decomp_profiler_transpose) nmax_timers = nmax_timers + 8
      if (decomp_profiler_io) nmax_timers = nmax_timers + 26
      if (decomp_profiler_fft) nmax_timers = nmax_timers + 7
      if (decomp_profiler_d2d) nmax_timers = nmax_timers + 2

      ! Allocate memory if needed
      if (nmax_timers > 0) then
         allocate (timer(nmax_timers))
         timer = 0.d0
         allocate (timer_start(nmax_timers))
         timer_start = 0.d0
         allocate (timer_n(nmax_timers))
         timer_n = 0
         allocate (timer_name(nmax_timers))
      end if

   end subroutine decomp_profiler_init_noarg

   !
   ! Finalize the basic timer module
   !
   subroutine decomp_profiler_fin_noarg()

      implicit none

      call timer_print()

      decomp_profiler = decomp_profiler_none

      ! Free memory if needed
      if (nmax_timers > 0) then
         ncur_timers = 0
         nmax_timers = 0
         deallocate (timer)
         deallocate (timer_start)
         deallocate (timer_n)
         deallocate (timer_name)
      end if

   end subroutine decomp_profiler_fin_noarg

   !
   ! Dummy log setup
   !
   subroutine decomp_profiler_log_int(io_unit)

      implicit none

      ! Argument
      integer, intent(in) :: io_unit

      if (nmax_timers == 0) then
         write (io_unit, *) "No profiling"
      else
         write (io_unit, *) "Generic profiling active"
      end if

   end subroutine decomp_profiler_log_int

   !
   ! Dummy setup
   !
   subroutine decomp_profiler_prep_bool(profiler_setup)

      implicit none

      logical, dimension(4), intent(in), optional :: profiler_setup

      decomp_profiler = decomp_profiler_none

      unused(profiler_setup)

   end subroutine decomp_profiler_prep_bool

   !
   ! Dummy start a timer
   !
   subroutine decomp_profiler_start_char(timer_name)

      implicit none

      ! Argument
      character(len=*), intent(in) :: timer_name

      timer_start(timer_find_or_create(timer_name)) = MPI_WTIME()

   end subroutine decomp_profiler_start_char

   !
   ! Dummy stop a timer
   !
   subroutine decomp_profiler_end_char(timer_name)

      implicit none

      ! Argument
      character(len=*), intent(in) :: timer_name

      ! Local variables
      integer :: id
      double precision :: deltaT

      ! Get the ID of the provided timer
      id = timer_find(timer_name)

      ! Update the timer
      deltaT = MPI_WTIME() - timer_start(id)
      timer(id) = timer(id) + deltaT
      timer_n(id) = timer_n(id) + 1

   end subroutine decomp_profiler_end_char

   !
   ! Try to find a timer with the provided name
   !
   ! Output :
   !    >0 : Success. Id of the timer.
   !    -1 : Failure
   !
   function timer_search(name) result(output)

      implicit none

      ! Arguments
      character(len=*), intent(in) :: name
      integer :: output

      ! Local variables
      integer :: id

      ! Safety check
      if (nmax_timers <= 0) call decomp_2d_abort(__FILE__, __LINE__, nmax_timers, "Invalid number of timers")

      ! Default value
      output = -1

      ! Try to find the given name
      do id = 1, ncur_timers
         if (trim(timer_name(id)) == trim(name)) then
            output = id
            return
         end if
      end do

   end function timer_search

   !
   ! Find or create a timer using the provided name
   !
   function timer_find_or_create(name) result(output)

      implicit none

      ! Arguments
      character(len=*), intent(in) :: name
      integer :: output

      ! Safety check
      if (nmax_timers <= 0) call decomp_2d_abort(__FILE__, __LINE__, nmax_timers, "Invalid number of timers")

      ! Try to find the given name
      output = timer_search(name)
      if (output > 0) return

      ! Create a new timer
      ncur_timers = ncur_timers + 1
      if (ncur_timers > nmax_timers) call decomp_2d_abort(__FILE__, __LINE__, ncur_timers, "Invalid number of timers")
      output = ncur_timers
      timer_name(output) (:) = ''
      timer_name(output) = trim(name)

   end function timer_find_or_create

   !
   ! Find a timer using the provided name
   !
   function timer_find(name) result(output)

      implicit none

      ! Arguments
      character(len=*), intent(in) :: name
      integer :: output

      ! Safety check
      if (nmax_timers <= 0) call decomp_2d_abort(__FILE__, __LINE__, nmax_timers, "Invalid number of timers")

      ! Try to find the given name
      output = timer_search(name)
      if (output > 0) return

      ! Timer not found, error
      call decomp_2d_abort(__FILE__, __LINE__, output, "Timer "//trim(name)//"not available")

   end function timer_find

   !
   ! Print the average / min / max for each timer
   !
   subroutine timer_print()

      implicit none

      ! Local variables
      integer :: io_unit, id
      double precision :: time, timer_min, timer_max, timer_avg

      ! Safety check
      if (ncur_timers <= 0) return

      ! Get the IO unit
      if (nrank == 0) open (newunit=io_unit, file='decomp_2d_perf.log', form='formatted')

      do id = 1, ncur_timers
         ! Compute min, max and average
         time = timer(id)
         call decomp_2d_mpi_allreduce(time, timer_min, MPI_MIN)
         call decomp_2d_mpi_allreduce(time, timer_max, MPI_MAX)
         call decomp_2d_mpi_allreduce(time, timer_avg, MPI_SUM)
         timer_avg = timer_avg / real(timer_n(id), real64) / real(nproc, real64)
         ! Print
         if (nrank == 0) then
            write (io_unit, *) "Timer "//trim(timer_name(id))//" avg, min, max"
            write (io_unit, *) "   ", real(timer_avg, 4), real(timer_min, 4), real(timer_max, 4)
         end if
      end do

      ! Close the IO unit
      if (nrank == 0) close (io_unit)

   end subroutine timer_print

end module decomp_2d_profiler
