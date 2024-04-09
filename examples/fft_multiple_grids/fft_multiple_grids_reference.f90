!! SPDX-License-Identifier: BSD-3-Clause

!
! Small module to provide the reference solution
!
! The input coordinate (i, j, k) is the global one
!
module reference_data

   public

contains

   ! This subroutines provides the reference real data
   pure function reference_rdata(i, j, k, nx, ny, nz)

      !$acc routine seq
      use decomp_2d_constants, only: mytype

      implicit none

      integer, intent(in) :: i, j, k, nx, ny, nz
      real(mytype) :: reference_rdata

      reference_rdata = real(i, mytype) / real(nx, mytype) &
                        * real(j, mytype) / real(ny, mytype) &
                        * real(k, mytype) / real(nz, mytype)

   end function reference_rdata

   ! This subroutine provides the reference complex data
   pure function reference_cdata(i, j, k, nx, ny, nz)

      !$acc routine seq
      use decomp_2d_constants, only: mytype

      implicit none

      integer, intent(in) :: i, j, k, nx, ny, nz
      complex(mytype) :: reference_cdata

      reference_cdata = cmplx(reference_rdata(i, j, k, nx, ny, nz), &
                              -3 * reference_rdata(i, j, k, nx, ny, nz), &
                              mytype)

   end function reference_cdata

end module reference_data
