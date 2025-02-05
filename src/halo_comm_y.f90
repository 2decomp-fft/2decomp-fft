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

! This file contain common code to be included by subroutines
! 'exchange_halo_y_...' in halo.f90

    if (present(opt_ylevel)) then
       level_x = opt_ylevel(1)
       level_y = opt_ylevel(2)
       level_z = opt_ylevel(3)
    else
       level_x = decomp%ylevel(1)
       level_y = decomp%ylevel(2)
       level_z = decomp%ylevel(3)
    end if

    s1 = decomp%ysz(1)
    s2 = decomp%ysz(2)
    s3 = decomp%ysz(3)

    xs = 1
    xe = s1 + 2 * level_x
    ys = 1 + level_y
    zs = 1
    ze = s3 + 2 * level_z

    ! if (decomp%halos_for_pencil) then
    !   ! don't communicate lower halo (south boundary)
    !   ys = 1 + level_y
    !   ! note north boundary not communicated either due to value of s2
    ! else
    !   ys = 1
    ! end if

#ifdef HALO_DEBUG
    if (nrank == 4) then
       write (*, *) 'Y-pencil input'
       write (*, *) '=============='
       write (*, *) 'Data on a x-z plane is shown'
       write (*, *) 'Before halo exchange'
       do i = xe, xs, -1
          write (*, '(10F4.0)') (inout(i, 1, k), k=zs, ze)
       end do
    end if
#endif

    ! *** east/west ***
    tag_w = coord(1)
    if (coord(1) == dims(1) - 1 .AND. periodic_x) then
       tag_e = 0
    else
       tag_e = coord(1) + 1
    end if
    icount = s2 * (s3 + 2 * level_z)
    ilength = level_x
    ijump = s1 + 2 * level_x
    call MPI_TYPE_VECTOR(icount, ilength, ijump, &
                         data_type, halo21, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_VECTOR")
    call MPI_TYPE_COMMIT(halo21, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
    ! receive from west
    call MPI_IRECV(inout(xs, ys, zs), 1, halo21, &
                   neighbour(2, 2), tag_w, DECOMP_2D_COMM_CART_Y, &
                   requests(1), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
    ! receive from east
    call MPI_IRECV(inout(xe - level_x + 1, ys, zs), 1, halo21, &
                   neighbour(2, 1), tag_e, DECOMP_2D_COMM_CART_Y, &
                   requests(2), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
    ! send to west
    call MPI_ISSEND(inout(xs + level_x, ys, zs), 1, halo21, &
                    neighbour(2, 2), tag_w, DECOMP_2D_COMM_CART_Y, &
                    requests(3), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
    ! send to east
    call MPI_ISSEND(inout(xe - level_x - level_x + 1, ys, zs), 1, halo21, &
                    neighbour(2, 1), tag_e, DECOMP_2D_COMM_CART_Y, &
                    requests(4), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
    call MPI_WAITALL(4, requests, status, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
    call MPI_TYPE_FREE(halo21, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")
#ifdef HALO_DEBUG
    if (nrank == 4) then
       write (*, *) 'After exchange in X'
       do i = xe, xs, -1
          write (*, '(10F4.0)') (inout(i, 1, k), k=zs, ze)
       end do
    end if
#endif

    ! *** north/south ***
    ! all data in local memory already, no halo exchange

    ! *** top/bottom ***
    ! no need to define derived data type as data on xy-planes
    ! all contiguous in memory, which can be sent/received using
    ! MPI directly
    tag_b = coord(2)
    if (coord(2) == dims(2) - 1 .AND. periodic_z) then
       tag_t = 0
    else
       tag_t = coord(2) + 1
    end if
    icount = (s2 * (s1 + 2 * level_x)) * level_z
    ! receive from bottom
    call MPI_IRECV(inout(xs, ys, zs), icount, data_type, &
                   neighbour(2, 6), tag_b, DECOMP_2D_COMM_CART_Y, &
                   requests(1), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
    ! receive from top
    call MPI_IRECV(inout(xs, ys, ze - level_z + 1), icount, data_type, &
                   neighbour(2, 5), tag_t, DECOMP_2D_COMM_CART_Y, &
                   requests(2), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
    ! send to bottom
    call MPI_ISSEND(inout(xs, ys, zs + level_z), icount, data_type, &
                    neighbour(2, 6), tag_b, DECOMP_2D_COMM_CART_Y, &
                    requests(3), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
    ! send to top
    call MPI_ISSEND(inout(xs, ys, ze - level_z - level_z + 1), icount, data_type, &
                    neighbour(2, 5), tag_t, DECOMP_2D_COMM_CART_Y, &
                    requests(4), ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
    call MPI_WAITALL(4, requests, status, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
#ifdef HALO_DEBUG
    if (nrank == 4) then
       write (*, *) 'After exchange in Z'
       do i = xe, xs, -1
          write (*, '(10F4.0)') (inout(i, 1, k), k=zs, ze)
       end do
    end if
#endif
