  !! halo_exchange_y_body.f90
  !!

#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'Y-pencil input'
     write (*, *) '=============='
     write (*, *) 'Data on a x-z plane is shown'
     write (*, *) 'Before halo exchange'
     do i = halo_extents%xe, halo_extents%xs, -1
        write (*, '(10F4.0)') (arr(i, 1, k), k=halo_extents%zs, halo_extents%ze)
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
  icount = s2 * (s3 + 2 * levels(3))
  ilength = levels(1)
  ijump = s1 + 2 * levels(1)
  call MPI_TYPE_VECTOR(icount, ilength, ijump, &
                       data_type, halo21, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_VECTOR")
  call MPI_TYPE_COMMIT(halo21, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
  ! receive from west
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs), 1, halo21, &
                 neighbour(2, 2), tag_w, DECOMP_2D_COMM_CART_Y, &
                 requests(1), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! receive from east
  call MPI_IRECV(arr(halo_extents%xe - levels(1) + 1, halo_extents%ys, halo_extents%zs), 1, halo21, &
                 neighbour(2, 1), tag_e, DECOMP_2D_COMM_CART_Y, &
                 requests(2), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! send to west
  call MPI_ISSEND(arr(halo_extents%xs + levels(1), halo_extents%ys, halo_extents%zs), 1, halo21, &
                  neighbour(2, 2), tag_w, DECOMP_2D_COMM_CART_Y, &
                  requests(3), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  ! send to east
  call MPI_ISSEND(arr(halo_extents%xe - levels(1) - levels(1) + 1, halo_extents%ys, halo_extents%zs), 1, halo21, &
                  neighbour(2, 1), tag_e, DECOMP_2D_COMM_CART_Y, &
                  requests(4), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  call MPI_WAITALL(4, requests, status, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
  call MPI_TYPE_FREE(halo21, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'After exchange in X'
     do i = halo_extents%xe, halo_extents%xs, -1
        write (*, '(10F4.0)') (arr(i, 1, k), k=halo_extents%zs, halo_extents%ze)
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
  icount = (s2 * (s1 + 2 * levels(1))) * levels(3)
  ! receive from bottom
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs), icount, data_type, &
                 neighbour(2, 6), tag_b, DECOMP_2D_COMM_CART_Y, &
                 requests(1), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! receive from top
  call MPI_IRECV(arr(halo_extents%xs, halo_extents%ys, halo_extents%ze - levels(3) + 1), icount, data_type, &
                 neighbour(2, 5), tag_t, DECOMP_2D_COMM_CART_Y, &
                 requests(2), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_IRECV")
  ! send to bottom
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ys, halo_extents%zs + levels(3)), icount, data_type, &
                  neighbour(2, 6), tag_b, DECOMP_2D_COMM_CART_Y, &
                  requests(3), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  ! send to top
  call MPI_ISSEND(arr(halo_extents%xs, halo_extents%ys, halo_extents%ze - levels(3) - levels(3) + 1), icount, data_type, &
                  neighbour(2, 5), tag_t, DECOMP_2D_COMM_CART_Y, &
                  requests(4), ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ISSEND")
  call MPI_WAITALL(4, requests, status, ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WAITALL")
#ifdef HALO_DEBUG
  if (nrank == 0) then
     write (*, *) 'After exchange in Z'
     do i = halo_extents%xe, halo_extents%xs, -1
        write (*, '(10F4.0)') (arr(i, 1, k), k=halo_extents%zs, halo_extents%ze)
     end do
  end if
#endif

