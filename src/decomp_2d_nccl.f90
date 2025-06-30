!! SPDX-License-Identifier: BSD-3-Clause

! Module for the cuda aware MPI

module decomp_2d_nccl

   use mpi
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_cumpi
   use cudafor
   use nccl

   implicit none

   private        ! Make everything private unless declared public

#ifdef DOUBLE_PREC
   type(ncclDataType), parameter, public :: ncclType = ncclDouble
#else
   type(ncclDataType), parameter, public :: ncclType = ncclFloat
#endif

   integer, save, public :: row_rank, col_rank

   integer, save, public :: col_comm_size, row_comm_size
   integer, allocatable, dimension(:), save, public :: local_to_global_col, local_to_global_row
   type(ncclUniqueId), save, public :: nccl_uid_2decomp
   type(ncclComm), save, public :: nccl_comm_2decomp
   integer(kind=cuda_stream_kind), save, public :: cuda_stream_2decomp

   ! Extra pointers for nccl complex transpose
   real(mytype), target, device, allocatable, dimension(:) :: work3, work4
   real(mytype), dimension(:), device, pointer, contiguous :: work3_r_d, work4_r_d

   public :: decomp_2d_nccl_init, &
             decomp_2d_nccl_fin, &
             decomp_2d_nccl_mem_init, &
             decomp_2d_nccl_mem_fin, &
             decomp_2d_nccl_send_recv_col, &
             decomp_2d_nccl_send_recv_row

   interface decomp_2d_nccl_send_recv_col
      module procedure decomp_2d_nccl_send_recv_real_col
      module procedure decomp_2d_nccl_send_recv_cmplx_col
   end interface decomp_2d_nccl_send_recv_col

   interface decomp_2d_nccl_send_recv_row
      module procedure decomp_2d_nccl_send_recv_real_row
      module procedure decomp_2d_nccl_send_recv_cmplx_row
   end interface decomp_2d_nccl_send_recv_row

contains
   !
   ! init of the arrays
   !
   subroutine decomp_2d_nccl_init(COMM_COL, COMM_ROW)

      implicit none

      integer, intent(in) :: COMM_COL, COMM_ROW
      integer             :: ierror
      integer :: cuda_stat
      type(ncclResult) :: nccl_stat

      call MPI_COMM_RANK(COMM_COL, col_rank, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
      call MPI_COMM_RANK(COMM_ROW, row_rank, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
      call MPI_COMM_SIZE(COMM_COL, col_comm_size, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
      call MPI_COMM_SIZE(COMM_ROW, row_comm_size, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")

      allocate (local_to_global_col(col_comm_size), local_to_global_row(row_comm_size))

      local_to_global_col(:) = 0
      local_to_global_row(:) = 0
      local_to_global_col(col_rank + 1) = nrank
      local_to_global_row(row_rank + 1) = nrank

      call mpi_allreduce(MPI_IN_PLACE, local_to_global_col, col_comm_size, MPI_INTEGER, MPI_SUM, COMM_COL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")
      call mpi_allreduce(MPI_IN_PLACE, local_to_global_row, row_comm_size, MPI_INTEGER, MPI_SUM, COMM_ROW, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")

      if (nrank == 0) then
         nccl_stat = ncclGetUniqueId(nccl_uid_2decomp)
         if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclGetUniqueId")
      end if
      call MPI_Bcast(nccl_uid_2decomp, int(sizeof(ncclUniqueId)), MPI_BYTE, 0, decomp_2d_comm, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

      nccl_stat = ncclCommInitRank(nccl_comm_2decomp, nproc, nccl_uid_2decomp, nrank)
      if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclCommInitRank")
      cuda_stat = cudaStreamCreate(cuda_stream_2decomp)
      if (cuda_stat /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "cudaStreamCreate")

   end subroutine decomp_2d_nccl_init

   !
   ! Finalize the module (release nccl communicator)
   !
   subroutine decomp_2d_nccl_fin()

      implicit none

      integer :: cuda_stat
      type(ncclResult) :: nccl_stat

      nccl_stat = ncclCommDestroy(nccl_comm_2decomp)
      if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclCommDestroy")
      cuda_stat = cudaStreamDestroy(cuda_stream_2decomp)
      if (cuda_stat /= 0) call decomp_2d_abort(__FILE__, __LINE__, cuda_stat, "cudaStreamDestroy")

   end subroutine decomp_2d_nccl_fin
   !
   ! Allocate the arrays
   !
   subroutine decomp_2d_nccl_mem_init(buf_size)

      use, intrinsic:: iso_c_binding, only: c_f_pointer, c_loc

      implicit none

      integer, intent(in) :: buf_size
      integer :: status, errorcode

      allocate (work3(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work4(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      if (associated(work3_r_d)) nullify (work3_r_d)
      if (associated(work4_r_d)) nullify (work4_r_d)
      call c_f_pointer(c_loc(work3), work3_r_d, [buf_size])
      call c_f_pointer(c_loc(work4), work4_r_d, [buf_size])

   end subroutine decomp_2d_nccl_mem_init
   !
   ! Free the arrays
   !
   subroutine decomp_2d_nccl_mem_fin

      implicit none

      if (associated(work3_r_d)) nullify (work3_r_d)
      if (associated(work4_r_d)) nullify (work4_r_d)
      if (allocated(work3)) deallocate (work3)
      if (allocated(work4)) deallocate (work4)

   end subroutine decomp_2d_nccl_mem_fin
   !
   ! Send-Recv Real Col
   !
   subroutine decomp_2d_nccl_send_recv_real_col(dst_d, &
                                                src_d, &
                                                disp_s, &
                                                cnts_s, &
                                                disp_r, &
                                                cnts_r, &
                                                dime)

      implicit none

      integer, intent(in) :: dime
      real(mytype), dimension(:), intent(in), device :: src_d
      real(mytype), dimension(:), intent(out), device :: dst_d
      integer, dimension(0:dime - 1), intent(in) :: disp_s, cnts_s
      integer, dimension(0:dime - 1), intent(in) :: disp_r, cnts_r

      integer :: col_rank_id, cuda_stat
      type(ncclResult) :: nccl_stat

      nccl_stat = ncclGroupStart()
      if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclGroupStart")
      do col_rank_id = 0, (col_comm_size - 1)
         nccl_stat = ncclSend(src_d(disp_s(col_rank_id) + 1), cnts_s(col_rank_id), &
                              ncclType, local_to_global_col(col_rank_id + 1), nccl_comm_2decomp, cuda_stream_2decomp)
         if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclSend")
         nccl_stat = ncclRecv(dst_d(disp_r(col_rank_id) + 1), cnts_r(col_rank_id), &
                              ncclType, local_to_global_col(col_rank_id + 1), nccl_comm_2decomp, cuda_stream_2decomp)
         if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclRecv")
      end do
      nccl_stat = ncclGroupEnd()
      if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclGroupEnd")
      cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
      if (cuda_stat /= 0) call decomp_2d_abort(__FILE__, __LINE__, cuda_stat, "cudaStreamSynchronize")

   end subroutine decomp_2d_nccl_send_recv_real_col
   !
   ! Send-Recv complex
   !
   subroutine decomp_2d_nccl_send_recv_cmplx_col(dst_d, &
                                                 src_d, &
                                                 disp_s, &
                                                 cnts_s, &
                                                 disp_r, &
                                                 cnts_r, &
                                                 dime, &
                                                 buf_size)

      implicit none

      integer, intent(in) :: dime, buf_size
      complex(mytype), dimension(buf_size), intent(in), device :: src_d
      complex(mytype), dimension(buf_size), intent(out), device :: dst_d
      integer, dimension(0:dime - 1), intent(in) :: disp_s, cnts_s
      integer, dimension(0:dime - 1), intent(in) :: disp_r, cnts_r

      integer :: ii

      ! Send-Recv Real part
      !$acc kernels default(present)
      do ii = 1, buf_size
         work3_r_d(ii) = real(src_d(ii), mytype)
      end do
      !$acc end kernels
      call decomp_2d_nccl_send_recv_col(work4_r_d, &
                                        work3_r_d, &
                                        disp_s, &
                                        cnts_s, &
                                        disp_r, &
                                        cnts_r, &
                                        dime)
      !$acc kernels default(present)
      do ii = 1, buf_size
         dst_d(ii) = cmplx(work4_r_d(ii), 0._mytype, mytype)
      end do
      !$acc end kernels
      ! Send-Recv Immaginary Part
      !$acc kernels default(present)
      do ii = 1, buf_size
         work3_r_d(ii) = aimag(src_d(ii))
      end do
      !$acc end kernels
      call decomp_2d_nccl_send_recv_col(work4_r_d, &
                                        work3_r_d, &
                                        disp_s, &
                                        cnts_s, &
                                        disp_r, &
                                        cnts_r, &
                                        dime)
      !$acc kernels default(present)
      do ii = 1, buf_size
         dst_d(ii) = cmplx(real(dst_d(ii), mytype), work4_r_d(ii), mytype)
      end do
      !$acc end kernels
   end subroutine decomp_2d_nccl_send_recv_cmplx_col
   !
   ! Send-Recv Real Row
   !
   subroutine decomp_2d_nccl_send_recv_real_row(dst_d, &
                                                src_d, &
                                                disp_s, &
                                                cnts_s, &
                                                disp_r, &
                                                cnts_r, &
                                                dime)

      implicit none

      integer, intent(in) :: dime
      real(mytype), dimension(:), intent(in), device :: src_d
      real(mytype), dimension(:), intent(out), device :: dst_d
      integer, dimension(0:dime - 1), intent(in) :: disp_s, cnts_s
      integer, dimension(0:dime - 1), intent(in) :: disp_r, cnts_r

      integer :: row_rank_id, cuda_stat
      type(ncclResult) :: nccl_stat

      nccl_stat = ncclGroupStart()
      if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclGroupStart")
      do row_rank_id = 0, (row_comm_size - 1)
         nccl_stat = ncclSend(src_d(disp_s(row_rank_id) + 1), cnts_s(row_rank_id), &
                              ncclType, local_to_global_row(row_rank_id + 1), nccl_comm_2decomp, cuda_stream_2decomp)
         if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclSend")
         nccl_stat = ncclRecv(dst_d(disp_r(row_rank_id) + 1), cnts_r(row_rank_id), &
                              ncclType, local_to_global_row(row_rank_id + 1), nccl_comm_2decomp, cuda_stream_2decomp)
         if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclRecv")
      end do
      nccl_stat = ncclGroupEnd()
      if (nccl_stat /= ncclSuccess) call decomp_2d_abort(__FILE__, __LINE__, nccl_stat, "ncclGroupEnd")
      cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
      if (cuda_stat /= 0) call decomp_2d_abort(__FILE__, __LINE__, cuda_stat, "cudaStreamSynchronize")

   end subroutine decomp_2d_nccl_send_recv_real_row
   !
   ! Send-Recv complex
   !
   subroutine decomp_2d_nccl_send_recv_cmplx_row(dst_d, &
                                                 src_d, &
                                                 disp_s, &
                                                 cnts_s, &
                                                 disp_r, &
                                                 cnts_r, &
                                                 dime, &
                                                 buf_size)

      implicit none

      integer, intent(in) :: dime, buf_size
      complex(mytype), dimension(buf_size), intent(in), device :: src_d
      complex(mytype), dimension(buf_size), intent(out), device :: dst_d
      integer, dimension(0:dime - 1), intent(in) :: disp_s, cnts_s
      integer, dimension(0:dime - 1), intent(in) :: disp_r, cnts_r

      integer :: ii

      ! Send-Recv Real part
      !$acc kernels default(present)
      do ii = 1, buf_size
         work3_r_d(ii) = real(src_d(ii), mytype)
      end do
      !$acc end kernels
      call decomp_2d_nccl_send_recv_row(work4_r_d, &
                                        work3_r_d, &
                                        disp_s, &
                                        cnts_s, &
                                        disp_r, &
                                        cnts_r, &
                                        dime)
      !$acc kernels default(present)
      do ii = 1, buf_size
         dst_d(ii) = cmplx(work4_r_d(ii), 0._mytype, mytype)
      end do
      !$acc end kernels
      ! Send-Recv Immaginary Part
      !$acc kernels default(present)
      do ii = 1, buf_size
         work3_r_d(ii) = aimag(src_d(ii))
      end do
      !$acc end kernels
      call decomp_2d_nccl_send_recv_row(work4_r_d, &
                                        work3_r_d, &
                                        disp_s, &
                                        cnts_s, &
                                        disp_r, &
                                        cnts_r, &
                                        dime)
      !$acc kernels default(present)
      do ii = 1, buf_size
         dst_d(ii) = cmplx(real(dst_d(ii), mytype), work4_r_d(ii), mytype)
      end do
      !$acc end kernels
   end subroutine decomp_2d_nccl_send_recv_cmplx_row

end module decomp_2d_nccl

