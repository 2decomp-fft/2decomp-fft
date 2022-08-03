!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to X pencil

  subroutine transpose_y_to_x_real(src, dst, opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#if defined(_GPU)
#if defined(_NCCL)
    integer :: col_rank_id
#endif
#endif

#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%COL_INFO%SND_P
    call mem_split_yx_real(src, s1, s2, s3, work1, dims(1), &
         decomp%y1dist, decomp)
#else

#if defined(_GPU)
    call mem_split_yx_real(src, s1, s2, s3, work1_r_d, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%y1dist, decomp)
#endif

#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%COL_INFO%RCV_P
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BARRIER")
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y1cnts_s, decomp%y1disp_s, &
            real_type, work2, decomp%x1cnts_s, decomp%x1disp_s, &
            real_type, decomp%COL_INFO%SMP_COMM, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
    end if
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%y1count, &
         real_type, work2_r, decomp%x1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
    nccl_stat = ncclGroupStart()
    do col_rank_id = 0, (col_comm_size - 1)
     nccl_stat = ncclSend(work1_r_d( decomp%y1disp(col_rank_id)+1 ), decomp%y1cnts(col_rank_id), ncclDouble, local_to_global_col(col_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
     nccl_stat = ncclRecv(work2_r_d( decomp%x1disp(col_rank_id)+1 ), decomp%x1cnts(col_rank_id), ncclDouble, local_to_global_col(col_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
    end do
    nccl_stat = ncclGroupEnd()
    cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
#else
    call MPI_ALLTOALLV(work1_r_d, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r_d, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BARRIER")
    call mem_merge_yx_real(work2, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else

#if defined(_GPU)
    call mem_merge_yx_real(work2_r_d, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif

#endif
    
    return
  end subroutine transpose_y_to_x_real


  subroutine transpose_y_to_x_complex(src, dst, opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#if defined(_GPU)
#if defined(_NCCL)
    integer :: col_rank_id
#endif
#endif

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%COL_INFO%SND_P_c
    call mem_split_yx_complex(src, s1, s2, s3, work1, dims(1), &
         decomp%y1dist, decomp)
#else

#if defined(_GPU)
    call mem_split_yx_complex(src, s1, s2, s3, work1_c_d, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_complex(src, s1, s2, s3, work1_c, dims(1), &
         decomp%y1dist, decomp)
#endif

#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p = decomp%COL_INFO%RCV_P_c
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BARRIER")
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y1cnts_s, decomp%y1disp_s, &
            complex_type, work2, decomp%x1cnts_s, decomp%x1disp_s, &
            complex_type, decomp%COL_INFO%SMP_COMM, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
    end if
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1_c, decomp%y1count, &
         complex_type, work2_c, decomp%x1count, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
    call MPI_ALLTOALLV(work1_c_d, decomp%y1cnts, decomp%y1disp, &
         complex_type, work2_c_d, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#else
    call MPI_ALLTOALLV(work1_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, work2_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BARRIER")
    call mem_merge_yx_complex(work2, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else

#if defined(_GPU)
    call mem_merge_yx_complex(work2_c_d, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_complex(work2_c, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif

#endif

    return
  end subroutine transpose_y_to_x_complex


  subroutine mem_split_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: out
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(pos), n1*(i2-i1+1), in(1,i1,1), n1*n2, n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_split_yx_real


  subroutine mem_split_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: out
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(pos), n1*(i2-i1+1), in(1,i1,1), n1*n2, n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_split_yx_complex


  subroutine mem_merge_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: in
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%x1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(i1,1,1), n1, in(pos), i2-i1+1, i2-i1+1, n2*n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_merge_yx_real


  subroutine mem_merge_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: in
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%x1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(i1,1,1), n1, in(pos), i2-i1+1, i2-i1+1, n2*n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_merge_yx_complex
