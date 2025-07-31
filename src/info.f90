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

   ! derived type to store decomposition info for a given global data size
   TYPE, extends(info), public :: DECOMP_INFO
      ! starting/ending index and size are defined in the parent type, see info.f90

      ! in addition to local information, processors also need to know
      ! some global information for global communications to work

      ! how each dimension is distributed along pencils
      integer, allocatable, dimension(:) :: &
         x1dist, y1dist, y2dist, z2dist

      ! send/receive buffer counts and displacements for MPI_ALLTOALLV
      integer, allocatable, dimension(:) :: &
         x1cnts, y1cnts, y2cnts, z2cnts
      integer, allocatable, dimension(:) :: &
         x1disp, y1disp, y2disp, z2disp

#ifdef EVEN
      ! buffer counts for MPI_ALLTOALL for padded-alltoall
      integer :: x1count, y1count, y2count, z2count
      ! evenly distributed data
      logical :: even
#endif

      ! Default halo extents per pencil
      integer, dimension(3) :: xlevel = [0, 0, 0]
      integer, dimension(3) :: ylevel = [0, 0, 0]
      integer, dimension(3) :: zlevel = [0, 0, 0]

   END TYPE DECOMP_INFO

   ! main (default) decomposition information for global size nx*ny*nz
   TYPE(DECOMP_INFO), target, save, public :: decomp_main

   ! Default : public
   public

end module m_info
