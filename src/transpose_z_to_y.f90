!! SPDX-License-Identifier: BSD-3-Clause

! This file contains the routines that transpose data from Z to Y pencil
submodule(decomp_2d) d2d_transpose_z_to_y

   use decomp_2d_constants, only: mytype

   implicit none

contains

   module subroutine transpose_z_to_y_real_short(src, dst)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_z_to_y(src, dst, decomp_main)

   end subroutine transpose_z_to_y_real_short

   module subroutine transpose_z_to_y_real_long(src, dst, decomp)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp

#if defined(_GPU)
      integer :: istat, nsize
#endif

      if (decomp_profiler_transpose) call decomp_profiler_start("transp_z_y_r")

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
         call transpose_z_to_y_real(src, dst, decomp, work1_r_d, work2_r_d)
#else
         if (use_pool) then
            call decomp_pool_get(work1_r)
            call decomp_pool_get(work2_r)
         end if
         call transpose_z_to_y_real(src, dst, decomp, work1_r, work2_r)
         if (use_pool) then
            call decomp_pool_free(work1_r)
            call decomp_pool_free(work2_r)
         end if
#endif
      end if

      if (decomp_profiler_transpose) call decomp_profiler_end("transp_z_y_r")

   end subroutine transpose_z_to_y_real_long

   subroutine transpose_z_to_y_real(src, dst, decomp, wk1, wk2)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN) :: src
      real(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(mytype), dimension(:), intent(out) :: wk1, wk2
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
      ! note the src array is suitable to be a send buffer
      ! so no split operation needed
#ifdef EVEN
      if (.not. decomp%even) then
         call mem_split_zy_real(src, s1, s2, s3, wk1, dims(2), &
                                decomp%z2dist, decomp)
      end if
# else
    ! Not EVEN branch
#   if defined(_NCCL)
      ! For NCCL comms we need to use a device array
      !$acc host_data use_device(src)
      istat = cudaMemcpy(wk1, src, s1 * s2 * s3, cudaMemcpyDeviceToDevice)
      !$acc end host_data
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#   endif
#endif


      ! define receive buffer
#ifdef EVEN
      if (decomp%even) then
         !$acc host_data use_device(src)
         call MPI_ALLTOALL(src, decomp%z2count, real_type, &
                           wk2, decomp%y2count, real_type, &
                           DECOMP_2D_COMM_ROW, ierror)
         !$acc end host_data
      else
         call MPI_ALLTOALL(wk1, decomp%z2count, real_type, &
                           wk2, decomp%y2count, real_type, &
                           DECOMP_2D_COMM_ROW, ierror)
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")

#else
#  if defined(_NCCL)
      call decomp_2d_nccl_send_recv_row(wk2, &
                                        wk1, &
                                        decomp%z2disp, &
                                        decomp%z2cnts, &
                                        decomp%y2disp, &
                                        decomp%y2cnts, &
                                        dims(2))
#  else
      associate (wk => wk1); end associate
      !$acc host_data use_device(src)
      call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, real_type, &
                         wk2, decomp%y2cnts, decomp%y2disp, real_type, &
                         DECOMP_2D_COMM_ROW, ierror)
      !$acc end host_data
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#  endif
#endif

      ! rearrange receive buffer
      call mem_merge_zy_real(wk2, d1, d2, d3, dst, dims(2), &
                             decomp%y2dist, decomp)

   end subroutine transpose_z_to_y_real

   module subroutine transpose_z_to_y_complex_short(src, dst)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst

      call transpose_z_to_y(src, dst, decomp_main)

   end subroutine transpose_z_to_y_complex_short

   module subroutine transpose_z_to_y_complex_long(src, dst, decomp)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: src
      complex(mytype), dimension(:, :, :), intent(OUT) :: dst
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      integer :: istat, nsize
#endif

      if (decomp_profiler_transpose) call decomp_profiler_start("transp_z_y_c")

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
         if (use_pool) then
            call decomp_pool_get(work1_c)
            call decomp_pool_get(work2_c)
         end if
         call transpose_z_to_y_complex(src, dst, decomp, work1_c, work2_c)
         if (use_pool) then
            call decomp_pool_free(work1_c)
            call decomp_pool_free(work2_c)
         end if
#endif
      end if

      if (decomp_profiler_transpose) call decomp_profiler_end("transp_z_y_c")

   end subroutine transpose_z_to_y_complex_long

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
#   if defined(_NCCL)
      ! For NCCL comms we need to use a device array
      !$acc host_data use_device(src)
      istat = cudaMemcpy(wk1, src, s1 * s2 * s3, cudaMemcpyDeviceToDevice)
      !$acc end host_data
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#   endif
#endif


      ! define receive buffer
#ifdef EVEN
      if (decomp%even) then
         !$acc host_data use_device(src)
         call MPI_ALLTOALL(src, decomp%z2count, complex_type, &
                           wk2, decomp%y2count, complex_type, &
                           DECOMP_2D_COMM_ROW, ierror)
         !$acc end host_data
      else
         call MPI_ALLTOALL(wk1, decomp%z2count, complex_type, &
                           wk2, decomp%y2count, complex_type, &
                           DECOMP_2D_COMM_ROW, ierror)
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else
    ! NO EVEN BRANCH
#   if defined(_NCCL)
      call decomp_2d_nccl_send_recv_row(wk2, &
                                        wk1, &
                                        decomp%z2disp, &
                                        decomp%z2cnts, &
                                        decomp%y2disp, &
                                        decomp%y2cnts, &
                                        dims(2), &
                                        decomp_buf_size)
#   else
      ! MPI
      associate (wk => wk1); end associate
      !$acc host_data use_device(src)
      call MPI_ALLTOALLV(src, decomp%z2cnts, decomp%z2disp, complex_type, &
                         wk2, decomp%y2cnts, decomp%y2disp, complex_type, &
                         DECOMP_2D_COMM_ROW, ierror)
      !$acc end host_data
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#   endif
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
      real(mytype), dimension(:), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: out
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos, init_pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         init_pos = m * decomp%z2count + 1
#else
         init_pos = decomp%z2disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(in)
         istat = cudaMemcpy2D(out(init_pos),      & !dst_addr
                              n1,                 & !dst_pitch
                              in(1,1,i1),         & !src_addr
                              n1,                 & !src_pitch
                              n1,                 & !width
                              n2 * (i2 - i1 + 1), & !height
                              cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         !$omp parallel do private(pos) collapse(3)
         do k = i1, i2
            do j = 1, n2
               do i = 1, n1
                  pos = init_pos + (i-1) + (j-1)*n1 + (k-i1)*n2*n1
                  out(pos) = in(i, j, k)
               end do
            end do
         end do
         !$omp end parallel do
#endif
      end do

   end subroutine mem_split_zy_real

   subroutine mem_split_zy_complex(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(n1, n2, n3), intent(IN) :: in
      complex(mytype), dimension(:), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: out
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos, init_pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         init_pos = m * decomp%z2count + 1
#else
         init_pos = decomp%z2disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(in)
         istat = cudaMemcpy2D(out(init_pos),      & !dst_addr
                              n1,                 & !dst_pitch
                              in(1,1,i1),         & !src_addr
                              n1,                 & !src_pitch
                              n1,                 & !width
                              n2 * (i2 - i1 + 1), & !height
                              cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         !$omp parallel do private(pos) collapse(3)
         do k = i1, i2
            do j = 1, n2
               do i = 1, n1
                  pos = init_pos + (i-1) + (j-1)*n1 + (k-i1)*n2*n1
                  out(pos) = in(i, j, k)
               end do
            end do
         end do
         !$omp end parallel do
#endif
      end do

   end subroutine mem_split_zy_complex

   subroutine mem_merge_zy_real(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      real(mytype), dimension(:), intent(IN) :: in
      real(mytype), dimension(n1, n2, n3), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
      attributes(device) :: in
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos, init_pos
      
      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         init_pos = m * decomp%y2count + 1
#else
         init_pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(out)
         istat = cudaMemcpy2D(out(1, i1, 1),      &
                              n1 * n2,            &
                              in(init_pos),       &
                              n1 * (i2 - i1 + 1), &
                              n1 * (i2 - i1 + 1), &
                              n3,                 &
                              cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         !$omp parallel do private(pos) collapse(3)
         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  pos = init_pos + (i-1) + (j-i1)*n1 + (k-1)*(i2-i1+1)*n1
                  out(i, j, k) = in(pos)
               end do
            end do
         end do
         !$omp end parallel do
#endif
      end do

   end subroutine mem_merge_zy_real

   subroutine mem_merge_zy_complex(in, n1, n2, n3, out, iproc, dist, decomp)

      implicit none

      integer, intent(IN) :: n1, n2, n3
      complex(mytype), dimension(:), intent(IN) :: in
      complex(mytype), dimension(n1, n2, n3), intent(OUT) :: out
      integer, intent(IN) :: iproc
      integer, dimension(0:iproc - 1), intent(IN) :: dist
      TYPE(DECOMP_INFO), intent(IN) :: decomp

#if defined(_GPU)
      attributes(device) :: in
      integer :: istat
#endif

      integer :: i, j, k, m, i1, i2, pos, init_pos

      do m = 0, iproc - 1
         if (m == 0) then
            i1 = 1
            i2 = dist(0)
         else
            i1 = i2 + 1
            i2 = i1 + dist(m) - 1
         end if

#ifdef EVEN
         init_pos = m * decomp%y2count + 1
#else
         init_pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
         !$acc host_data use_device(out)
         istat = cudaMemcpy2D(out(1, i1, 1),      &
                              n1 * n2,            &
                              in(init_pos),       &
                              n1 * (i2 - i1 + 1), &
                              n1 * (i2 - i1 + 1), &
                              n3,                 &
                              cudaMemcpyDeviceToDevice)
         !$acc end host_data
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
         !$omp parallel do private(pos) collapse(3)
         do k = 1, n3
            do j = i1, i2
               do i = 1, n1
                  pos = init_pos + (i-1) + (j-i1)*n1 + (k-1)*(i2-i1+1)*n1
                  out(i, j, k) = in(pos)
               end do
            end do
         end do
         !$omp end parallel do
#endif
      end do

   end subroutine mem_merge_zy_complex

end submodule d2d_transpose_z_to_y
