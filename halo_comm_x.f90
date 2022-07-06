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
! 'exchange_halo_x_...' in halo.f90

    if (present(opt_decomp)) then
      decomp = opt_decomp
    else
      decomp = decomp_main
    end if

    if (present(opt_level)) then
      level_x = 0
      level_y = opt_level
      level_z = opt_level
    else
      level_x = decomp%xlevel(1)
      level_y = decomp%xlevel(2)
      level_z = decomp%xlevel(3)
    end if

    s1 = decomp%xsz(1)
    s2 = decomp%xsz(2)
    s3 = decomp%xsz(3)

    xs = 1 + level_x
    ys = 1
    ye = s2 + 2*level_y
    zs = 1
    ze = s3 + 2*level_z

    ! if (decomp%halos_for_pencil) then
    !   ! don't communicate lower halo (west boundary)
    !   xs = 1 + level_x
    !   ! note east boundary not communicated either due to value of s1
    ! else
    !   xs = 1
    ! end if

#ifdef HALO_DEBUG
    if (nrank==4) then
      write(*,*) 'X-pencil input'
      write(*,*) '=============='
      write(*,*) 'Data on a y-z plane is shown'
      write(*,*) 'Before halo exchange'
      do j=ye,ys,-1
        write(*,'(10F4.0)') (inout(1,j,k),k=zs,ze)
      end do
    end if
#endif

    ! *** east/west ***
    ! all data in local memory already, no halo exchange

    ! *** north/south ***
    tag_s = coord(1)
    if (coord(1)==dims(1)-1 .AND. periodic_y) then
      tag_n = 0
    else
      tag_n = coord(1) + 1
    end if
    icount = s3 + 2*level_z
    ilength = level_y * s1
    ijump = s1 * (s2 + 2*level_y)
    call MPI_TYPE_VECTOR(icount, ilength, ijump, &
      data_type, halo12, ierror)
    call MPI_TYPE_COMMIT(halo12, ierror)
    ! receive from south
    call MPI_IRECV(inout(xs,ys,zs), 1, halo12, &
      neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
      requests(1), ierror)
    ! receive from north
    call MPI_IRECV(inout(xs,ye-level_y+1,zs), 1, halo12, &
      neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
      requests(2), ierror)
    ! send to south
    call MPI_ISSEND(inout(xs,ys+level_y,zs), 1, halo12, &
      neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
      requests(3), ierror)
    ! send to north
    call MPI_ISSEND(inout(xs,ye-level_y-level_y+1,zs), 1, halo12, &
      neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
      requests(4), ierror)
    call MPI_WAITALL(4, requests, status, ierror)
    call MPI_TYPE_FREE(halo12, ierror)
#ifdef HALO_DEBUG
    if (nrank==4) then
      write(*,*) 'After exchange in Y'
      do j=ye,ys,-1
        write(*,'(10F4.0)') (inout(1,j,k),k=zs,ze)
      end do
    end if
#endif

    ! *** top/bottom ***
    ! no need to define derived data type as data on xy-planes
    ! all contiguous in memory, which can be sent/received using
    ! MPI directly
    tag_b = coord(2)
    if (coord(2)==dims(2)-1 .AND. periodic_z) then
      tag_t = 0
    else
      tag_t = coord(2) + 1
    end if
    icount = (s1 * (s2 + 2*level_y)) * level_z
    ! receive from bottom
    call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
      neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
      requests(1), ierror)
    ! receive from top
    call MPI_IRECV(inout(xs,ys,ze-level_z+1), icount, data_type, &
      neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
      requests(2), ierror)
    ! send to bottom
    call MPI_ISSEND(inout(xs,ys,zs+level_z), icount, data_type, &
      neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
      requests(3), ierror)
    ! send to top
    call MPI_ISSEND(inout(xs,ys,ze-level_z-level_z+1), icount, data_type, &
      neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
      requests(4), ierror)
    call MPI_WAITALL(4, requests, status, ierror)

#ifdef HALO_DEBUG
    if (nrank==4) then
      write(*,*) 'After exchange in Z'
      do j=ye,ys,-1
        write(*,'(10F4.0)') (inout(1,j,k),k=zs,ze)
      end do
    end if
#endif
