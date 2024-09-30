!! SPDX-License-Identifier: BSD-3-Clause

! Preprocessor macro to deal with unused variables
#define unused(x) associate(tmp => x); end associate

!
! Dummy module when there is no ADIOS2 IO available
!
module decomp_2d_io_adios

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io_utilities
   use decomp_2d_mpi
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

   ! The external code can use this variable to check if adios2 is available
   logical, parameter, public :: decomp_2d_with_adios2 = .false.

   private

   public :: decomp_2d_io_adios_init, &
             decomp_2d_io_adios_fin

contains

   !
   ! Initialize the dummy ADIOS2 IO module
   !
   subroutine decomp_2d_io_adios_init(adios_xml)

      implicit none

      character(len=*), intent(in), optional :: adios_xml

      unused(adios_xml)

   end subroutine decomp_2d_io_adios_init

   !
   ! Finalize the dummy ADIOS2 IO module
   !
   subroutine decomp_2d_io_adios_fin()

      implicit none

   end subroutine decomp_2d_io_adios_fin

end module decomp_2d_io_adios
