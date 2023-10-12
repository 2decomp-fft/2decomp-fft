!! SPDX-License-Identifier: BSD-3-Clause

! Module for the cuda aware MPI

module decomp_2d_cumpi

   use decomp_2d_constants
   use decomp_2d_mpi

   implicit none

   private        ! Make everything private unless declared public

   ! Real/complex pointers to GPU buffers
   real(mytype), dimension(:), device, public, pointer :: work1_r_d, work2_r_d
   complex(mytype), dimension(:), device, public, pointer :: work1_c_d, work2_c_d

   public :: decomp_2d_cumpi_init, &
             decomp_2d_cumpi_fin

contains
   !
   ! init of the arrays
   !
   subroutine decomp_2d_cumpi_init(buf_size, wk1, wk2)

      implicit none

      integer, intent(in) :: buf_size
      real(mytype), target, dimension(:), device :: wk1, wk2

      if (associated(work1_r_d)) nullify (work1_r_d)
      if (associated(work2_r_d)) nullify (work2_r_d)
      if (associated(work1_c_d)) nullify (work1_c_d)
      if (associated(work2_c_d)) nullify (work2_c_d)
      call c_f_pointer(c_loc(wk1), work1_r_d, [buf_size])
      call c_f_pointer(c_loc(wk2), work2_r_d, [buf_size])
      call c_f_pointer(c_loc(wk1), work1_c_d, [buf_size])                                            
      call c_f_pointer(c_loc(wk2), work2_c_d, [buf_size])

   end subroutine decomp_2d_cumpi_init
   !
   ! init of the arrays
   !
   subroutine decomp_2d_cumpi_fin

      implicit none

      if (associated(work1_r_d)) nullify (work1_r_d)
      if (associated(work2_r_d)) nullify (work2_r_d)                                                 
      if (associated(work1_c_d)) nullify (work1_c_d)
      if (associated(work2_c_d)) nullify (work2_c_d)

   end subroutine decomp_2d_cumpi_fin

end module decomp_2d_cumpi

