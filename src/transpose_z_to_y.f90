!! SPDX-License-Identifier: BSD-3-Clause

! This file contains the routines that transpose data from Z to Y pencil
submodule(decomp_2d) d2d_transpose_z_to_y

   use decomp_2d_constants, only: mytype

   implicit none

contains

   ! Short transpose interface for real numbers with implicit decomposition.
   module subroutine transpose_z_to_y_real_short(src, dst)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_z_to_y(src, dst, decomp_main)

   end subroutine transpose_z_to_y_real_short
    
   ! Full transpose interface for real numbers with explicit decomposition
   module subroutine transpose_z_to_y_real_long(src, dst, decomp)

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
      if (decomp_profiler_transpose) call decomp_profiler_start("transp_z_y_r")
#endif

#if defined(_GPU)
      wk1 => work1_r_d
      wk2 => work2_r_d
#else
      wk1 => work1_r
      wk2 => work2_r
#endif

      if (dims(2) == 1) then
         call transpose_z_to_y_real_fast(src, dst, decomp)
      else
         call transpose_z_to_y_real(src, dst, decomp, wk1, wk2)
      end if

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_end("transp_z_y_r")
#endif

   end subroutine transpose_z_to_y_real_long

   ! Fast implementation of transpose for real numbers avoiding communication
   subroutine transpose_z_to_y_real_fast(src, dst, decomp)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: nsize
#if defined(_GPU)
      integer :: istat
#endif

      nsize = product(decomp%zsz)
#if defined(_GPU)
      !$acc host_data use_device(src,dst)
      istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
      !$acc end host_data
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
      dst = src
#endif

   end subroutine transpose_z_to_y_real_fast

   ! Default implementation of transpose for real numbers
   subroutine transpose_z_to_y_real(src, dst, decomp, wk1, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(out) :: wk1, wk2
#if defined(_GPU)
      attributes(device) :: wk1, wk2
#endif

      call rearrange_sendbuf_real(src, wk1, decomp)
      call perform_transpose_real(src, decomp, wk1, wk2)
      call rearrange_recvbuf_real(dst, decomp, wk2)

   end subroutine transpose_z_to_y_real

   ! Short transpose interface for complex numbers with implicit decomposition
   module subroutine transpose_z_to_y_complex_short(src, dst)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_z_to_y(src, dst, decomp_main)

   end subroutine transpose_z_to_y_complex_short

   ! Full transpose interface for complex numbers with explicit decomposition
   module subroutine transpose_z_to_y_complex_long(src, dst, decomp)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      integer :: istat, nsize
#endif

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_start("transp_z_y_c")
#endif

      if (dims(2) == 1) then
#if defined(_GPU)
         nsize = product(decomp%zsz)
         !$acc host_data use_device(src,dst)
         istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
         dst = src
#endif
      else
#if defined(_GPU)
         call transpose_z_to_y_complex(src, dst, decomp, work1_c_d, work2_c_d)
#else
         call transpose_z_to_y_complex(src, dst, decomp, work1_c, work2_c)
#endif
      end if

#ifdef PROFILER
      if (decomp_profiler_transpose) call decomp_profiler_end("transp_z_y_c")
#endif

   end subroutine transpose_z_to_y_complex_long

   ! Transpose implementation for complex numbers
   subroutine transpose_z_to_y_complex(src, dst, decomp, wk1, wk2)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      complex(mytype), dimension(:), intent(out) :: wk1, wk2
#if defined(_GPU)
      attributes(device) :: wk1, wk2

      integer :: istat
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
#ifdef EVEN
      if (.not. decomp%even) then
         call mem_split_zy_complex(src, s1, s2, s3, wk1, dims(2), &
                                   decomp%z2dist, decomp)
      end if
#else
      ! note the src array is suitable to be a send buffer
      ! so no split operation needed

#if defined(_GPU)
      !$acc host_data use_device(src)
      istat = cudaMemcpy(wk1, src, s1 * s2 * s3, cudaMemcpyDeviceToDevice)
      !$acc end host_data
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#endif

#endif

      ! define receive buffer
#ifdef EVEN
      if (decomp%even) then
         call MPI_ALLTOALL(src, decomp%z2count, complex_type, &
                           wk2, decomp%y2count, complex_type, &
                           DECOMP_2D_COMM_ROW, ierror)
      else
         call MPI_ALLTOALL(wk1, decomp%z2count, complex_type, &
                           wk2, decomp%y2count, complex_type, &
                           DECOMP_2D_COMM_ROW, ierror)
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
      call decomp_2d_nccl_send_recv_row(wk2, &
                                        wk1, &
                                        decomp%z2disp, &
                                        decomp%z2cnts, &
                                        decomp%y2disp, &
                                        decomp%y2cnts, &
                                        dims(2), &
                                        decomp_buf_size)
#else
      call MPI_ALLTOALLV(wk1, decomp%z2cnts, decomp%z2disp, complex_type, &
                         wk2, decomp%y2cnts, decomp%y2disp, complex_type, &
                         DECOMP_2D_COMM_ROW, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
      associate (wk => wk1)
      end associate
      call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, complex_type, &
                         wk2, decomp%y2cnts, decomp%y2disp, complex_type, &
                         DECOMP_2D_COMM_ROW, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

      ! rearrange receive buffer
      call mem_merge_zy_complex(wk2, d1, d2, d3, dst, dims(2), &
                                decomp%y2dist, decomp)

   end subroutine transpose_z_to_y_complex

   ! pack/unpack ALLTOALL(V) buffers

   subroutine mem_split_zy_real(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      real(mytype), dimension(n1, n2, n3), intent(IN) :: in
      real(mytype), dimension(*), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

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
         pos = m * decomp%z2count + 1
#else
         pos = decomp%z2disp(m) + 1
#endif

         do k = i1, i2
            do j = 1, n2
               do i = 1, n1
                  out(pos) = in(i, j, k)
                  pos = pos + 1
               end do
            end do
         end do
      end do

   end subroutine mem_split_zy_real

   subroutine mem_split_zy_complex(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(n1, n2, n3), intent(IN) :: in
      complex(mytype), dimension(*), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

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
         pos = m * decomp%z2count + 1
#else
         pos = decomp%z2disp(m) + 1
#endif

         do k = i1, i2
            do j = 1, n2
               do i = 1, n1
                  out(pos) = in(i, j, k)
                  pos = pos + 1
               end do
            end do
         end do
      end do

   end subroutine mem_split_zy_complex

   subroutine mem_merge_zy_real(in, n1, n2, n3, out, iproc, dist, decomp)

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
         pos = m * decomp%y2count + 1
#else
         pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(out)
         istat = cudaMemcpy2D(out(1, i1, 1), n1 * n2, in(pos), n1 * (i2 - i1 + 1), n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  out(i, j, k) = in(pos)
                  pos = pos + 1
               end do
            end do
         end do
#endif
      end do

   end subroutine mem_merge_zy_real

   subroutine mem_merge_zy_complex(in, n1, n2, n3, out, iproc, dist, decomp)

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
         pos = m * decomp%y2count + 1
#else
         pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(out)
         istat = cudaMemcpy2D(out(1, i1, 1), n1 * n2, in(pos), n1 * (i2 - i1 + 1), n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  out(i, j, k) = in(pos)
                  pos = pos + 1
               end do
            end do
         end do
#endif
      end do

   end subroutine mem_merge_zy_complex

!!! Utility subroutines

   ! Pack the array to be transposed into the send buffer as necessary.
   subroutine rearrange_sendbuf_real(src, wk1, decomp)

      real(mytype), dimension(:, :, :), intent(IN) :: src
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(out) :: wk1
#if defined(_GPU)
      attributes(device) :: wk1

      integer :: istat
#endif

      integer :: s1, s2, s3

      s1 = SIZE(src, 1)
      s2 = SIZE(src, 2)
      s3 = SIZE(src, 3)

#ifdef EVEN
      if (.not. decomp%even) then
         call mem_split_zy_real(src, s1, s2, s3, wk1, dims(2), &
                                decomp%z2dist, decomp)
      end if
#else
      ! note the src array is suitable to be a send buffer
      ! so no split operation needed

#if defined(_GPU)
      !$acc host_data use_device(src)
      istat = cudaMemcpy(wk1, src, s1 * s2 * s3, cudaMemcpyDeviceToDevice)
      !$acc end host_data
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
      associate(foo => wk1, bar => decomp)
      end associate
#endif

#endif
   end subroutine rearrange_sendbuf_real

   ! Subroutine to do the transpose communications.
   subroutine perform_transpose_real(src, decomp, wk1, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2
#if defined(_GPU)
      attributes(device) :: wk1, wk2
#endif
     
#ifdef EVEN
      call perform_transpose_real_even(src, decomp, wk1, wk2)
#else

#if defined(_GPU)
      call perform_transpose_real_gpu(decomp, wk1, wk2)
#else
      associate (wk => wk1)
      end associate
      call perform_transpose_real_std(src, decomp, wk2)
#endif
#endif

   end subroutine perform_transpose_real

   ! Unpack the transposed data from the receive buffer into the destination array.
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

      call mem_merge_zy_real(wk2, d1, d2, d3, dst, dims(2), &
           decomp%y2dist, decomp)
      
   end subroutine rearrange_recvbuf_real

!!!! Option-specific utility subroutines
   
#ifndef EVEN
   
#ifndef _GPU

   ! Subroutine to do the transpose communications in the default case.
   subroutine perform_transpose_real_std(src, decomp, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(out) :: wk2

      integer :: ierror

      call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, real_type, &
                         wk2, decomp%y2cnts, decomp%y2disp, real_type, &
                         DECOMP_2D_COMM_ROW, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
      
   end subroutine perform_transpose_real_std

#else
   ! _GPU is defined

   ! Subroutine to do the transpose communications on GPUs.
   subroutine perform_transpose_real_gpu(decomp, wk1, wk2)

      implicit none

      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2
      attributes(device) :: wk1, wk2

#ifndef _NCCL
      integer :: istat
      integer :: ierror
      
      call MPI_ALLTOALLV(wk1, decomp%z2cnts, decomp%z2disp, real_type, &
                         wk2, decomp%y2cnts, decomp%y2disp, real_type, &
                         DECOMP_2D_COMM_ROW, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#else
      call decomp_2d_nccl_send_recv_row(wk2, &
                                        wk1, &
                                        decomp%z2disp, &
                                        decomp%z2cnts, &
                                        decomp%y2disp, &
                                        decomp%y2cnts, &
                                        dims(2))
#endif

   end subroutine perform_transpose_real_gpu

   ! END IFDEF _GPU
#endif

#else
   ! EVEN is defined

   ! Subroutine to do the transpose communications for even communications.
   subroutine perform_transpose_real_even(src, decomp, wk1, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(in) :: wk1
      real(mytype), dimension(:), intent(out) :: wk2

      integer :: ierror

      if (decomp%even) then
         call MPI_ALLTOALL(src, decomp%z2count, real_type, &
                           wk2, decomp%y2count, real_type, &
                           DECOMP_2D_COMM_ROW, ierror)
      else
         call MPI_ALLTOALL(wk1, decomp%z2count, real_type, &
                           wk2, decomp%y2count, real_type, &
                           DECOMP_2D_COMM_ROW, ierror)
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
      
   end subroutine perform_transpose_real_even

   ! END IFDEF EVEN
#endif

end submodule d2d_transpose_z_to_y
