!! SPDX-License-Identifier: BSD-3-Clause

!
! This is the module for simple IO utilities
!

module decomp_2d_io_utilities

   use decomp_2d_constants
   use decomp_2d_mpi
   use m_info
   use MPI

   implicit none

   private

   public :: io_get_size

contains

   !
   ! Extract the size / start / subsize
   !
   subroutine io_get_size(ipencil, decomp, sizes, starts, subsizes, nplanes)

      implicit none

      integer, intent(in) :: ipencil
      class(info), intent(in) :: decomp
      integer, dimension(3), intent(out) :: sizes, subsizes, starts
      integer, intent(in), optional :: nplanes

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if
      if (present(nplanes)) then
         if (nplanes < 1) then
            call decomp_2d_abort(nplanes, "Invalid value for nplanes")
         end if
      end if

      sizes(1) = decomp%xsz(1)
      sizes(2) = decomp%ysz(2)
      sizes(3) = decomp%zsz(3)
      if (ipencil == 1) then
         starts(:) = decomp%xst(:) - 1
         subsizes(:) = decomp%xsz(:)
      else if (ipencil == 2) then
         starts(:) = decomp%yst(:) - 1
         subsizes(:) = decomp%ysz(:)
      else
         starts(:) = decomp%zst(:) - 1
         subsizes(:) = decomp%zsz(:)
      end if

      if (present(nplanes)) then
         sizes(ipencil) = nplanes
         subsizes(ipencil) = nplanes
      end if

   end subroutine io_get_size

end module decomp_2d_io_utilities
