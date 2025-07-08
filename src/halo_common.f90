!! SPDX-License-Identifier: BSD-3-Clause

! This file contain common code to be included by subroutines
! 'update_halo_...' in halo.f90

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    s1 = size(in, 1)
    s2 = size(in, 2)
    s3 = size(in, 3)

    ipencil = get_pencil([s1, s2, s3], decomp, opt_pencil)
    halo_extents = halo_extents_t(ipencil, [s1, s2, s3], decomp, level, global)
    allocate (out(halo_extents%xs:halo_extents%xe, &
                  halo_extents%ys:halo_extents%ye, &
                  halo_extents%zs:halo_extents%ze))
    !    out = -1.0_mytype ! fill the halo for debugging

    !$acc enter data create(requests,neighbour)
    ! copy input data to output
    if (global) then
       ! using global coordinate
       ! note the input array passed in always has index starting from 1
       ! need to work out the corresponding global index
       if (ipencil == 1) then
          kst = decomp%xst(3); ken = decomp%xen(3)
          jst = decomp%xst(2); jen = decomp%xen(2)
          !$acc kernels default(present)
          do k = kst, ken
             do j = jst, jen
                do i = 1, s1  ! x all local
                   out(i, j, k) = in(i, j - decomp%xst(2) + 1, k - decomp%xst(3) + 1)
                end do
             end do
          end do
          !$acc end kernels
       else if (ipencil == 2) then
          kst = decomp%yst(3); ken = decomp%yen(3)
          ist = decomp%yst(1); ien = decomp%yen(1)
          !$acc kernels default(present)
          do k = kst, ken
             do j = 1, s2  ! y all local
                do i = ist, ien
                   out(i, j, k) = in(i - decomp%yst(1) + 1, j, k - decomp%yst(3) + 1)
                end do
             end do
          end do
          !$acc end kernels
       else if (ipencil == 3) then
          jst = decomp%zst(2); jen = decomp%zen(2)
          ist = decomp%zst(1); ien = decomp%zen(1)
          !$acc kernels default(present)
          do k = 1, s3  ! z all local
             do j = jst, jen
                do i = ist, ien
                   out(i, j, k) = in(i - decomp%zst(1) + 1, j - decomp%zst(2) + 1, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
    else
       ! not using global coordinate
       !$acc kernels default(present)
       do k = 1, s3
          do j = 1, s2
             do i = 1, s1
                out(i, j, k) = in(i, j, k)
             end do
          end do
       end do
       !$acc end kernels
       !!! istat = cudaMemcpy(out,in,s1*s2*s3,cudaMemcpyDeviceToDevice)
    end if

    ! If needed, define MPI derived data type to pack halo data,
    ! then call MPI send/receive to exchange halo data

    if (ipencil == 1) then
#include "halo_exchange_x_body.f90"
    else if (ipencil == 2) then
#include "halo_exchange_y_body.f90"
    else if (ipencil == 3) then
#include "halo_exchange_z_body.f90"
    end if  ! pencil
