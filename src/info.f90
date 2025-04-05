!! SPDX-License-Identifier: BSD-3-Clause

!
! Provide a simplified decomp_info object
!
module m_info

   implicit none

   ! Abstract type, decomp_info will extend it
   type, abstract :: info
      ! staring/ending index and size of data held by current processor
      integer, dimension(3) :: xst, xen, xsz ! x-pencil
      integer, dimension(3) :: yst, yen, ysz ! y-pencil
      integer, dimension(3) :: zst, zen, zsz ! z-pencil
   end type info

   ! Default : public
   public

end module m_info
