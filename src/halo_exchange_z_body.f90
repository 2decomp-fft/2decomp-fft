  !! halo_exchange_z_body.f90
  !!
  
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'Z-pencil input'
     write (*, *) '=============='
     write (*, *) 'Data on a x-y plane is shown'
     write (*, *) 'Before halo exchange'
     do i = halo_extents%xe, halo_extents%xs, -1
        write (*, '(10F4.0)') (arr(i, j, 1), j=halo_extents%ys, halo_extents%ye)
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
  icount = (s2 + 2 * level) * s3
  ilength = level
  ijump = s1 + 2 * level
  call MPI_TYPE_VECTOR(icount, ilength, ijump, &
       data_type, halo31, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_VECTOR")
  call MPI_TYPE_COMMIT(halo31, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
  ! receive from west
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs), 1, halo31, &
       neighbour(3, 2), tag_w, DECOMP_2D_COMM_CART_Z, &
       requests(1), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! receive from east
  call MPI_IRECV(arr(halo_extents%xe - level + 1, halo_extents%ys, halo_extents%zs), 1, halo31, &
       neighbour(3, 1), tag_e, DECOMP_2D_COMM_CART_Z, &
       requests(2), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! send to west
  call MPI_ISSEND(arr(halo_extents%xs + level, halo_extents%ys, halo_extents%zs), 1, halo31, &
       neighbour(3, 2), tag_w, DECOMP_2D_COMM_CART_Z, &
       requests(3), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  ! send to east
  call MPI_ISSEND(arr(halo_extents%xe - level - level + 1, halo_extents%ys, halo_extents%zs), 1, halo31, &
       neighbour(3, 1), tag_e, DECOMP_2D_COMM_CART_Z, &
       requests(4), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  call MPI_WAITALL(4, requests, status, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
  call MPI_TYPE_FREE(halo31, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'After exchange in X'
     do i = halo_extents%xe, halo_extents%xs, -1
        write (*, '(10F4.0)') (arr(i, j, 1), j=halo_extents%ys, halo_extents%ye)
     end do
  end if
#endif

  ! *** north/south ***
  tag_s = coord(2)
  if (coord(2) == dims(2) - 1 .AND. periodic_y) then
     tag_n = 0
  else
     tag_n = coord(2) + 1
  end if
  icount = s3
  ilength = level * (s1 + 2 * level)
  ijump = (s1 + 2 * level) * (s2 + 2 * level)
  call MPI_TYPE_VECTOR(icount, ilength, ijump, &
       data_type, halo32, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_VECTOR")
  call MPI_TYPE_COMMIT(halo32, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
  ! receive from south
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs), 1, halo32, &
       neighbour(3, 4), tag_s, DECOMP_2D_COMM_CART_Z, &
       requests(1), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! receive from north
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ye - level + 1, halo_extents%zs), 1, halo32, &
       neighbour(3, 3), tag_n, DECOMP_2D_COMM_CART_Z, &
       requests(2), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! send to south
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ys + level, halo_extents%zs), 1, halo32, &
       neighbour(3, 4), tag_s, DECOMP_2D_COMM_CART_Z, &
       requests(3), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  ! send to north
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ye - level - level + 1, halo_extents%zs), 1, halo32, &
       neighbour(3, 3), tag_n, DECOMP_2D_COMM_CART_Z, &
       requests(4), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  call MPI_WAITALL(4, requests, status, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
  call MPI_TYPE_FREE(halo32, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'After exchange in Y'
     do i = halo_extents%xe, halo_extents%xs, -1
        write (*, '(10F4.0)') (arr(i, j, 1), j=halo_extents%ys, halo_extents%ye)
     end do
  end if
#endif

  ! *** top/bottom ***
  ! all data in local memory already, no halo exchange

