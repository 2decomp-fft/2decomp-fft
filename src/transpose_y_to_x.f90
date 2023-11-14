!! SPDX-License-Identifier: BSD-3-Clause

! This file contains the routines that transpose data from Y to X pencil
submodule(decomp_2d) d2d_transpose_y_to_x

   use decomp_2d_constants, only: mytype

   use transpose_common
   
   implicit none

contains

   ! Short transpose interface for real numbers with implicit decomposition.
   module subroutine transpose_y_to_x_real_short(src, dst)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_y_to_x(src, dst, decomp_main)

   end subroutine transpose_y_to_x_real_short

   ! Full transpose interface for real numbers with explicit decomposition
   module subroutine transpose_y_to_x_real_long(src, dst, decomp)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      real(mytype), dimension(:), pointer, contiguous :: wk1
      real(mytype), dimension(:), pointer, contiguous :: wk2
#ifdef _GPU
      attributes(device) :: wk1, wk2
#endif
      

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_start("transp_y_x_r")
#endif

#ifdef _GPU
      wk1 => work1_r_d
      wk2 => work2_r_d
#else
      wk1 => work1_r
      wk2 => work2_r
#endif

      if (dims(1) == 1) then
         call transpose_real_fast(src, dst, product(decomp%ysz))
      else
         call transpose_y_to_x_real(src, dst, decomp, wk1, wk2)
      end if

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_end("transp_y_x_r")
#endif

    end subroutine transpose_y_to_x_real_long

    ! Default implementation of transpose for real numbers
    subroutine transpose_y_to_x_real(src, dst, decomp, wk1, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(out) :: wk1, wk2
#if defined(_GPU) && defined(_NCCL)
      attributes(device) :: wk1, wk2
#endif

      integer :: s1, s2, s3, d1, d2, d3

      s1 = SIZE(src, 1)
      s2 = SIZE(src, 2)
      s3 = SIZE(src, 3)
      d1 = SIZE(dst, 1)
      d2 = SIZE(dst, 2)
      d3 = SIZE(dst, 3)

      call rearrange_sendbuf_real(src, wk1, decomp)
      call perform_transpose_real(decomp, wk1, wk2)
      call rearrange_recvbuf_real(dst, decomp, wk2)

   end subroutine transpose_y_to_x_real

   module subroutine transpose_y_to_x_complex_short(src, dst)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_y_to_x(src, dst, decomp_main)

   end subroutine transpose_y_to_x_complex_short

   module subroutine transpose_y_to_x_complex_long(src, dst, decomp)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp

#if defined(_GPU)
      integer :: istat, nsize
#endif

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_start("transp_y_x_c")
#endif

      if (dims(1) == 1) then
#if defined(_GPU)
         nsize = product(decomp%ysz)
         !$acc host_data use_device(src,dst)
         istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
         dst = src
#endif
      else
#if defined(_GPU)
         call transpose_y_to_x_complex(src, dst, decomp, work1_c_d, work2_c_d)
#else
         call transpose_y_to_x_complex(src, dst, decomp, work1_c, work2_c)
#endif
      end if

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_end("transp_y_x_c")
#endif

   end subroutine transpose_y_to_x_complex_long

   subroutine transpose_y_to_x_complex(src, dst, decomp, wk1, wk2)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      complex(mytype), dimension(:), intent(OUT) :: wk1, wk2
#if defined(_GPU)
      attributes(device) :: wk1, wk2
#endif

      integer :: s1, s2, s3, d1, d2, d3
      integer :: ierror

      s1 = SIZE(src, 1)
      s2 = SIZE(src, 2)
      s3 = SIZE(src, 3)
      d1 = SIZE(dst, 1)
      d2 = SIZE(dst, 2)
      d3 = SIZE(dst, 3)

      ! rearrange source array as send buffer
      call mem_split_yx_complex(src, s1, s2, s3, wk1, dims(1), &
                                decomp%y1dist, decomp)

      ! define receive buffer
      ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
      call MPI_ALLTOALL(wk1, decomp%y1count, complex_type, &
                        wk2, decomp%x1count, complex_type, &
                        DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU) && defined(_NCCL)
      call decomp_2d_nccl_send_recv_col(wk2, &
                                        wk1, &
                                        decomp%y1disp, &
                                        decomp%y1cnts, &
                                        decomp%x1disp, &
                                        decomp%x1cnts, &
                                        dims(1), &
                                        decomp_buf_size)
#else
      call MPI_ALLTOALLV(wk1, decomp%y1cnts, decomp%y1disp, complex_type, &
                         wk2, decomp%x1cnts, decomp%x1disp, complex_type, &
                         DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

      ! rearrange receive buffer
      call mem_merge_yx_complex(wk2, d1, d2, d3, dst, dims(1), &
                                decomp%x1dist, decomp)

   end subroutine transpose_y_to_x_complex

   subroutine mem_split_yx_real(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      real(mytype), dimension(n1, n2, n3), intent(IN) :: in
      real(mytype), dimension(*), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: out
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%y1count + 1
#else
         pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(in)
         istat = cudaMemcpy2D(out(pos), n1 * (i2 - i1 + 1), in(1, i1, 1), n1 * n2, n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  out(pos) = in(i, j, k)
                  pos = pos + 1
               end do
            end do
         end do
#endif
      end do

   end subroutine mem_split_yx_real

   subroutine mem_split_yx_complex(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(n1, n2, n3), intent(IN) :: in
      complex(mytype), dimension(*), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: out
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%y1count + 1
#else
         pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(in)
         istat = cudaMemcpy2D(out(pos), n1 * (i2 - i1 + 1), in(1, i1, 1), n1 * n2, n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  out(pos) = in(i, j, k)
                  pos = pos + 1
               end do
            end do
         end do
#endif
      end do

   end subroutine mem_split_yx_complex

   subroutine mem_merge_yx_real(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      real(mytype), dimension(*), intent(IN) :: in
      real(mytype), dimension(n1, n2, n3), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: in
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%x1count + 1
#else
         pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(out)
         istat = cudaMemcpy2D(out(i1, 1, 1), n1, in(pos), i2 - i1 + 1, i2 - i1 + 1, n2 * n3, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         do k = 1, n3
            do j = 1, n2
               do i = i1, i2
                  out(i, j, k) = in(pos)
                  pos = pos + 1
               end do
            end do
         end do
#endif
      end do

   end subroutine mem_merge_yx_real

   subroutine mem_merge_yx_complex(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(*), intent(IN) :: in
      complex(mytype), dimension(n1, n2, n3), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: in
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         pos = m * decomp%x1count + 1
#else
         pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(out)
         istat = cudaMemcpy2D(out(i1, 1, 1), n1, in(pos), i2 - i1 + 1, i2 - i1 + 1, n2 * n3, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         do k = 1, n3
            do j = 1, n2
               do i = i1, i2
                  out(i, j, k) = in(pos)
                  pos = pos + 1
               end do
            end do
         end do
#endif
      end do

   end subroutine mem_merge_yx_complex

!!! Utility subroutines

   ! Pack the array to be transposed into the send buffer as necessary.
   subroutine rearrange_sendbuf_real(src, wk1, decomp)

     real(mytype), dimension(:, :, :), intent(IN) :: src
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     real(mytype), dimension(:), intent(out) :: wk1
#if defined(_GPU)
     attributes(device) :: wk1
#endif

     integer :: s1, s2, s3

     s1 = SIZE(src, 1)
     s2 = SIZE(src, 2)
     s3 = SIZE(src, 3)

     call mem_split_yx_real(src, s1, s2, s3, wk1, dims(1), &
                            decomp%y1dist, decomp)
      
   end subroutine rearrange_sendbuf_real

   ! Subroutine to do the transpose communications.
   subroutine perform_transpose_real(decomp, wk1, wk2)

      implicit none

      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2
#if defined(_GPU) && defined(_NCCL)
      attributes(device) :: wk1, wk2
#endif
      
#ifdef EVEN
      call perform_transpose_real_even(decomp, wk1, wk2)
#else

#if defined(_GPU) && defined(_NCCL)
      call perform_transpose_real_gpu(decomp, wk1, wk2)
#else
      call perform_transpose_real_std(decomp, wk1, wk2)
#endif

#endif
      
   end subroutine perform_transpose_real

   subroutine rearrange_recvbuf_real(dst, decomp, wk2)

      implicit none
     
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk2
#if defined(_GPU)
      attributes(device) :: wk2
#endif

      integer :: d1, d2, d3

      d1 = SIZE(dst, 1)
      d2 = SIZE(dst, 2)
      d3 = SIZE(dst, 3)

      call mem_merge_yx_real(wk2, d1, d2, d3, dst, dims(1), &
                             decomp%x1dist, decomp)
      
   end subroutine rearrange_recvbuf_real

!!! Option-specific utility subroutines

#ifndef EVEN

#if defined(_GPU) && defined(_NCCL)

   ! Subroutine to do the transpose communications on GPUs.
   subroutine perform_transpose_real_gpu(decomp, wk1, wk2)

      implicit none

      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2
      attributes(device) :: wk1, wk2

      call decomp_2d_nccl_send_recv_col(wk2, &
                                        wk1, &
                                        decomp%y1disp, &
                                        decomp%y1cnts, &
                                        decomp%x1disp, &
                                        decomp%x1cnts, &
                                        dims(1))
      
   end subroutine perform_transpose_real_gpu
   
#else
   ! At least one of _GPU and _NCCL is not defined

   ! Subroutine to do the transpose communications in the default case.
   subroutine perform_transpose_real_std(decomp, wk1, wk2)

      implicit none

      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(IN) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2

      integer :: ierror

      call MPI_ALLTOALLV(wk1, decomp%y1cnts, decomp%y1disp, real_type, &
                         wk2, decomp%x1cnts, decomp%x1disp, real_type, &
                         DECOMP_2D_COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
      
   end subroutine perform_transpose_real_std
   
#endif
   ! ifdef _GPU
   
#else
   ! EVEN is defined

   ! Subroutine to do the transpose communications for even communications
   subroutine perform_transpose_real_even(decomp, wk1, wk2)

      implicit none

      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2

      integer :: ierror

     call MPI_ALLTOALL(wk1, decomp%y1count, real_type, &
                       wk2, decomp%x1count, real_type, &
                       DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")

    end subroutine perform_transpose_real_even
   
#endif
   ! ifdef EVEN
   
end submodule d2d_transpose_y_to_x
