  !! halo_exchange_x_body.f90
  !!
  
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'X-pencil input'
     write (*, *) '=============='
     write (*, *) 'Data on a y-z plane is shown'
     write (*, *) 'Before halo exchange'
     do j = halo_extents%ye, halo_extents%ys, -1
        write (*, '(10F4.0)') (arr(1, j, k), k=halo_extents%zs, halo_extents%ze)
     end do
  end if
#endif

  ! *** east/west ***
  ! all data in local memory already, no halo exchange

  ! *** north/south ***
  tag_s = coord(1)
  if (coord(1) == dims(1) - 1 .AND. periodic_y) then
     tag_n = 0
  else
     tag_n = coord(1) + 1
  end if
  icount = s3 + 2 * level
  ilength = level * s1
  ijump = s1 * (s2 + 2 * level)
  call MPI_TYPE_VECTOR(icount, ilength, ijump, &
       data_type, halo12, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_VECTOR")
  call MPI_TYPE_COMMIT(halo12, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
  ! receive from south
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs), 1, halo12, &
       neighbour(1, 4), tag_s, DECOMP_2D_COMM_CART_X, &
       requests(1), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! receive from north
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ye - level + 1, halo_extents%zs), 1, halo12, &
       neighbour(1, 3), tag_n, DECOMP_2D_COMM_CART_X, &
       requests(2), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! send to south
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ys + level, halo_extents%zs), 1, halo12, &
       neighbour(1, 4), tag_s, DECOMP_2D_COMM_CART_X, &
       requests(3), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  ! send to north
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ye - level - level + 1, halo_extents%zs), 1, halo12, &
       neighbour(1, 3), tag_n, DECOMP_2D_COMM_CART_X, &
       requests(4), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  call MPI_WAITALL(4, requests, status, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
  call MPI_TYPE_FREE(halo12, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'After exchange in Y'
     do j = halo_extents%ye, halo_extents%ys, -1
        write (*, '(10F4.0)') (arr(1, j, k), k=halo_extents%zs, halo_extents%ze)
     end do
  end if
#endif

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
  icount = (s1 * (s2 + 2 * level)) * level
  ! receive from bottom
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs), icount, data_type, &
       neighbour(1, 6), tag_b, DECOMP_2D_COMM_CART_X, &
       requests(1), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! receive from top
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%ze - level + 1), icount, data_type, &
       neighbour(1, 5), tag_t, DECOMP_2D_COMM_CART_X, &
       requests(2), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! send to bottom
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs + level), icount, data_type, &
       neighbour(1, 6), tag_b, DECOMP_2D_COMM_CART_X, &
       requests(3), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  ! send to top
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ys, halo_extents%ze - level - level + 1), icount, data_type, &
       neighbour(1, 5), tag_t, DECOMP_2D_COMM_CART_X, &
       requests(4), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  call MPI_WAITALL(4, requests, status, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'After exchange in Z'
     do j = halo_extents%ye, halo_extents%ys, -1
        write (*, '(10F4.0)') (arr(1, j, k), k=halo_extents%zs, halo_extents%ze)
     end do
  end if
#endif


