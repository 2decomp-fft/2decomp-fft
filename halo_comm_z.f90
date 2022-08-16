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
! 'exchange_halo_z_...' in halo.f90

    if (present(opt_decomp)) then
      decomp = opt_decomp
    else
      decomp = decomp_main
    end if

    if (present(opt_zlevel)) then ! assume same level for all directions
      level_x = opt_zlevel(1)
      level_y = opt_zlevel(2)
      level_z = opt_zlevel(3)
    else
      level_x = decomp%zlevel(1)
      level_y = decomp%zlevel(2)
      level_z = decomp%zlevel(3)
      ! add checks to make sure level_x and level_y are sensible values (positive integer)?
    end if

    s1 = decomp%zsz(1)
    s2 = decomp%zsz(2)
    s3 = decomp%zsz(3)

    xs = 1
    xe = s1 + 2*level_x
    ys = 1
    ye = s2 + 2*level_y
    zs = 1 + level_z

    ! if (decomp%halos_for_pencil) then
    !   ! don't communicate lower halo (bottom boundary)
    !   zs = 1 + level_z
    !   ! note top boundary not communicated either due to value of s3
    ! else
    !   zs = 1
    ! end if

#ifdef HALO_DEBUG
    if (nrank==4) then
      write(*,*) 'Z-pencil input'
      write(*,*) '=============='
      write(*,*) 'Data on a x-y plane is shown'
      write(*,*) 'Before halo exchange'
      do i=xe,xs,-1
        write(*,'(10F4.0)') (inout(i,j,1),j=ys,ye)
      end do
    end if
#endif
    ! *** east/west ***
    tag_w = coord(1)
    if (coord(1)==dims(1)-1 .AND. periodic_x) then
      tag_e = 0
    else
      tag_e = coord(1) + 1
    end if
    icount = (s2 + 2*level_y) * s3
    ilength = level_x
    ijump = s1 + 2*level_x
    call MPI_TYPE_VECTOR(icount, ilength, ijump, &
        data_type, halo31, ierror)
    call MPI_TYPE_COMMIT(halo31, ierror)
    ! receive from west
    call MPI_IRECV(inout(xs,ys,zs), 1, halo31, &
        neighbour(3,2), tag_w, DECOMP_2D_COMM_CART_Z, &
        requests(1), ierror)
    ! receive from east
    call MPI_IRECV(inout(xe-level_x+1,ys,zs), 1, halo31, &
        neighbour(3,1), tag_e, DECOMP_2D_COMM_CART_Z, &
        requests(2), ierror)
    ! send to west
    call MPI_ISSEND(inout(xs+level_x,ys,zs), 1, halo31, &
        neighbour(3,2), tag_w, DECOMP_2D_COMM_CART_Z, &
        requests(3), ierror)
    ! send to east
    call MPI_ISSEND(inout(xe-level_x-level_x+1,ys,zs), 1, halo31, &
        neighbour(3,1), tag_e, DECOMP_2D_COMM_CART_Z, &
        requests(4), ierror)
    call MPI_WAITALL(4, requests, status, ierror)
    call MPI_TYPE_FREE(halo31, ierror)

#ifdef HALO_DEBUG
    if (nrank==4) then
      write(*,*) 'After exchange in X'
      do i=xe,xs,-1
        write(*,'(10F4.0)') (inout(i,j,1),j=ys,ye)
      end do
    end if
#endif

    ! *** north/south ***
    tag_s = coord(2)
    if (coord(2)==dims(2)-1 .AND. periodic_y) then
      tag_n = 0
    else
      tag_n = coord(2) + 1
    end if
    icount = s3
    ilength = level_y * (s1 + 2*level_x)
    ijump = (s1 + 2*level_x) * (s2 + 2*level_y)
    call MPI_TYPE_VECTOR(icount, ilength, ijump, &
        data_type, halo32, ierror)
    call MPI_TYPE_COMMIT(halo32, ierror)
    ! receive from south
    call MPI_IRECV(inout(xs,ys,zs), 1, halo32, &
        neighbour(3,4), tag_s, DECOMP_2D_COMM_CART_Z, &
        requests(1), ierror)
    ! receive from north
    call MPI_IRECV(inout(xs,ye-level_y+1,zs), 1, halo32, &
        neighbour(3,3), tag_n, DECOMP_2D_COMM_CART_Z, &
        requests(2), ierror)
    ! send to south
    call MPI_ISSEND(inout(xs,ys+level_y,zs), 1, halo32, &
        neighbour(3,4), tag_s, DECOMP_2D_COMM_CART_Z, &
        requests(3), ierror)
    ! send to north
    call MPI_ISSEND(inout(xs,ye-level_y-level_y+1,zs), 1, halo32, &
        neighbour(3,3), tag_n, DECOMP_2D_COMM_CART_Z, &
        requests(4), ierror)
    call MPI_WAITALL(4, requests, status, ierror)
    call MPI_TYPE_FREE(halo32, ierror)
#ifdef HALO_DEBUG
    if (nrank==4) then
      write(*,*) 'After exchange in Y'
      do i=xe,xs,-1
        write(*,'(10F4.0)') (inout(i,j,1),j=ys,ye)
      end do
    end if
#endif

    ! *** top/bottom ***
    ! all data in local memory already, no halo exchange
