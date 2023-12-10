!! SPDX-License-Identifier: BSD-3-Clause

! This file contains 3D c2c/r2c/c2r transform subroutines which are
! identical for several FFT engines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D FFT - complex to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_c2c(in, out, isign)

   implicit none

   ! Arguments
   complex(mytype), dimension(:, :, :), intent(INOUT) :: in
   complex(mytype), dimension(:, :, :), intent(OUT) :: out
   integer, intent(IN) :: isign

   ! Local variables
   integer :: i, j, k
   complex(mytype), allocatable, dimension(:, :, :) :: wk1

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_start("fft_c2c")
#endif

   if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
       format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

      ! ===== 1D FFTs in X =====
      if (inplace) then
         call c2c_1m_x(in, isign, ph)
      else
         call alloc_x(wk1, ph)
         do concurrent(k=1:ph%xsz(3), j=1:ph%xsz(2), i=1:ph%xsz(1))
            wk1(i, j, k) = in(i, j, k)
         end do
         call c2c_1m_x(wk1, isign, ph)
      end if

      ! ===== Swap X --> Y; 1D FFTs in Y =====
      if (dims(1) > 1) then
         if (inplace) then
            call transpose_x_to_y(in, wk2_c2c, ph)
         else
            call transpose_x_to_y(wk1, wk2_c2c, ph)
         end if
         call c2c_1m_y(wk2_c2c, isign, ph)
      else
         if (inplace) then
            call c2c_1m_y(in, isign, ph)
         else
            call c2c_1m_y(wk1, isign, ph)
         end if
      end if

      ! ===== Swap Y --> Z; 1D FFTs in Z =====
      if (dims(1) > 1) then
         call transpose_y_to_z(wk2_c2c, out, ph)
      else
         if (inplace) then
            call transpose_y_to_z(in, out, ph)
         else
            call transpose_y_to_z(wk1, out, ph)
         end if
      end if
      call c2c_1m_z(out, isign, ph)

   else if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_BACKWARD &
            .OR. &
            format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_FORWARD) then

      ! ===== 1D FFTs in Z =====
      if (inplace) then
         call c2c_1m_z(in, isign, ph)
      else
         call alloc_z(wk1, ph)
         do concurrent(k=1:ph%zsz(3), j=1:ph%zsz(2), i=1:ph%zsz(1))
            wk1(i, j, k) = in(i, j, k)
         end do
         call c2c_1m_z(wk1, isign, ph)
      end if

      ! ===== Swap Z --> Y; 1D FFTs in Y =====
      if (dims(1) > 1) then
         if (inplace) then
            call transpose_z_to_y(in, wk2_c2c, ph)
         else
            call transpose_z_to_y(wk1, wk2_c2c, ph)
         end if
         call c2c_1m_y(wk2_c2c, isign, ph)
      else  ! out==wk2_c2c if 1D decomposition
         if (inplace) then
            call transpose_z_to_y(in, out, ph)
         else
            call transpose_z_to_y(wk1, out, ph)
         end if
         call c2c_1m_y(out, isign, ph)
      end if

      ! ===== Swap Y --> X; 1D FFTs in X =====
      if (dims(1) > 1) then
         call transpose_y_to_x(wk2_c2c, out, ph)
      end if
      call c2c_1m_x(out, isign, ph)

   end if

   ! Free memory
   if (allocated(wk1)) deallocate (wk1)

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_end("fft_c2c")
#endif

   return
end subroutine fft_3d_c2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_r2c(in_r, out_c)

   implicit none

   real(mytype), dimension(:, :, :), intent(IN) :: in_r
   complex(mytype), dimension(:, :, :), intent(OUT) :: out_c

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")
#endif

   if (format == PHYSICAL_IN_X) then

      ! ===== 1D FFTs in X =====
      call r2c_1m_x(in_r, wk13)

      ! ===== Swap X --> Y; 1D FFTs in Y =====
      if (dims(1) > 1) then
         call transpose_x_to_y(wk13, wk2_r2c, sp)
         call c2c_1m_y(wk2_r2c, -1, sp)
      else
         call c2c_1m_y(wk13, -1, sp)
      end if

      ! ===== Swap Y --> Z; 1D FFTs in Z =====
      if (dims(1) > 1) then
         call transpose_y_to_z(wk2_r2c, out_c, sp)
      else
         call transpose_y_to_z(wk13, out_c, sp)
      end if
      call c2c_1m_z(out_c, -1, sp)

   else if (format == PHYSICAL_IN_Z) then

      ! ===== 1D FFTs in Z =====
      call r2c_1m_z(in_r, wk13)

      ! ===== Swap Z --> Y; 1D FFTs in Y =====
      if (dims(1) > 1) then
         call transpose_z_to_y(wk13, wk2_r2c, sp)
         call c2c_1m_y(wk2_r2c, -1, sp)
      else  ! out_c==wk2_r2c if 1D decomposition
         call transpose_z_to_y(wk13, out_c, sp)
         call c2c_1m_y(out_c, -1, sp)
      end if

      ! ===== Swap Y --> X; 1D FFTs in X =====
      if (dims(1) > 1) then
         call transpose_y_to_x(wk2_r2c, out_c, sp)
      end if
      call c2c_1m_x(out_c, -1, sp)

   end if

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_end("fft_r2c")
#endif

   return
end subroutine fft_3d_r2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_c2r(in_c, out_r)

   implicit none

   ! Arguments
   complex(mytype), dimension(:, :, :), intent(INOUT) :: in_c
   real(mytype), dimension(:, :, :), intent(OUT) :: out_r

   ! Local variables
   integer :: i, j, k
   complex(mytype), allocatable, dimension(:, :, :) :: wk1

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_start("fft_c2r")
#endif

   if (format == PHYSICAL_IN_X) then

      ! ===== 1D FFTs in Z =====
      if (inplace) then
         call c2c_1m_z(in_c, 1, sp)
      else
         call alloc_z(wk1, sp)
         do concurrent(k=1:sp%zsz(3), j=1:sp%zsz(2), i=1:sp%zsz(1))
            wk1(i, j, k) = in_c(i, j, k)
         end do
         call c2c_1m_z(wk1, 1, sp)
      end if

      ! ===== Swap Z --> Y; 1D FFTs in Y =====
      if (inplace) then
         call transpose_z_to_y(in_c, wk2_r2c, sp)
      else
         call transpose_z_to_y(wk1, wk2_r2c, sp)
      end if
      call c2c_1m_y(wk2_r2c, 1, sp)

      ! ===== Swap Y --> X; 1D FFTs in X =====
      if (dims(1) > 1) then
         call transpose_y_to_x(wk2_r2c, wk13, sp)
         call c2r_1m_x(wk13, out_r)
      else
         call c2r_1m_x(wk2_r2c, out_r)
      end if

   else if (format == PHYSICAL_IN_Z) then

      ! ===== 1D FFTs in X =====
      if (inplace) then
         call c2c_1m_x(in_c, 1, sp)
      else
         call alloc_x(wk1, sp)
         do concurrent(k=1:sp%xsz(3), j=1:sp%xsz(2), i=1:sp%xsz(1))
            wk1(i, j, k) = in_c(i, j, k)
         end do
         call c2c_1m_x(wk1, 1, sp)
      end if

      ! ===== Swap X --> Y; 1D FFTs in Y =====
      if (dims(1) > 1) then
         if (inplace) then
            call transpose_x_to_y(in_c, wk2_r2c, sp)
         else
            call transpose_x_to_y(wk1, wk2_r2c, sp)
         end if
         call c2c_1m_y(wk2_r2c, 1, sp)
      else  ! in_c==wk2_r2c if 1D decomposition
         if (inplace) then
            call c2c_1m_y(in_c, 1, sp)
         else
            call c2c_1m_y(wk1, 1, sp)
         end if
      end if

      ! ===== Swap Y --> Z; 1D FFTs in Z =====
      if (dims(1) > 1) then
         call transpose_y_to_z(wk2_r2c, wk13, sp)
      else
         if (inplace) then
            call transpose_y_to_z(in_c, wk13, sp)
         else
            call transpose_y_to_z(wk1, wk13, sp)
         end if
      end if
      call c2r_1m_z(wk13, out_r)

   end if

   ! Free memory
   if (allocated(wk1)) deallocate (wk1)

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_end("fft_c2r")
#endif

   return
end subroutine fft_3d_c2r
