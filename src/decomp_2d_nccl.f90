!=======================================================================
! This is part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
! Copyright (C) 2022               Science and Technology Facilities Council(STFC)
!
!=======================================================================

! Module for the cuda aware MPI

module decomp_2d_nccl

   use mpi
   use decomp_2d_constants  
   use decomp_2d_cumpi
   use nccl

   implicit none
   
   private        ! Make everything private unless declared public

#ifdef DOUBLE_PREC
   type(ncclDataType),public :: ncclType = ncclDouble
#else
   type(ncclDataType),public :: ncclType = ncclFloat
#endif
   
   integer, save, public :: row_rank, col_rank
   
   integer, save, public :: col_comm_size, row_comm_size
   integer, allocatable, dimension(:), save, public :: local_to_global_col, local_to_global_row
   type(ncclUniqueId), save, public :: nccl_uid_2decomp
   type(ncclComm), save, public :: nccl_comm_2decomp
   integer(kind=cuda_stream_kind), save, public :: cuda_stream_2decomp
             
   public :: decomp_2d_nccl_init, &
             decomp_2d_nccl_fin, &
             decomp_2d_nccl_send_recv_real_col, &
             decomp_2d_nccl_send_recv_real_row

contains
   !
   ! init of the arrays
   ! 
   subroutine decomp_2d_nccl_init(COMM_COL,COMM_ROW)
      
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
   ! init of the arrays
   ! 
   subroutine decomp_2d_nccl_fin()
      
      implicit none
      
      type(ncclResult) :: nccl_stat

      nccl_stat = ncclCommDestroy(nccl_comm_2decomp)

   end subroutine decomp_2d_nccl_fin
  !
  ! Send-Recv Real Col 
  !
  subroutine decomp_2d_nccl_send_recv_real_col(dst_d,  &
                                               src_d,  &
                                               disp_s, &
                                               cnts_s, &
                                               disp_r, &
                                               cnts_r, &
                                               dime    )
     
     implicit none

     integer, intent(in) :: dime      
     real(mytype), dimension(:), intent(in), device :: src_d
     real(mytype), dimension(:), intent(out),device :: dst_d
     integer, dimension(0:dime-1), intent(in) :: disp_s, cnts_s
     integer, dimension(0:dime-1), intent(in) :: disp_r, cnts_r


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
  
  end subroutine decomp_2d_nccl_send_recv_real_col
  !
  ! Send-Recv Real Row
  !
  subroutine decomp_2d_nccl_send_recv_real_row(dst_d,  &
                                               src_d,  &
                                               disp_s, &
                                               cnts_s, &
                                               disp_r, &
                                               cnts_r, &
                                               dime    )
     
     implicit none
     
     integer, intent(in) :: dime      
     real(mytype), dimension(:), intent(in), device :: src_d
     real(mytype), dimension(:), intent(out),device :: dst_d
     integer, dimension(0:dime-1), intent(in) :: disp_s, cnts_s
     integer, dimension(0:dime-1), intent(in) :: disp_r, cnts_r


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
  
  end subroutine decomp_2d_nccl_send_recv_real_row
  
  

end module decomp_2d_nccl

