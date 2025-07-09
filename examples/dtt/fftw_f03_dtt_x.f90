!! SPDX-License-Identifier: BSD-3-Clause

program dtt_x

   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI

   implicit none

   integer, parameter :: nx_base = 8, ny_base = 9, nz_base = 10
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: dttxi, dttxf, dttyi, dttyf, dttzi, dttzf

   ! Decomp_info objects in the phuysical and spectral space
   type(decomp_info), pointer :: ph => null(), sp => null()
   ! Ouput when there is no periodicity
   real(mytype), target, allocatable, dimension(:, :, :) :: out_r
   ! Input
   real(mytype), target, allocatable, dimension(:, :, :) :: in_r
   ! Objects
   type(decomp_2d_fft_engine), target, save :: dtt_engine

   ! 3 directions
   ! 8 transforms (4 DCT, 4 DST)
   ! Default values for ifirst, ofirst and ndismiss
   integer, dimension(3, 9, 9, 9) :: DTT

   integer :: ierror, j, k, l, ii, jj, kk
   integer :: st1, st2, st3
   integer :: en1, en2, en3
   real(mytype) :: error, errl2, errli

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   dttxi = -1
   dttyi = -1
   dttzi = -1
   ! Now we can check if user put some inputs
   call decomp_2d_testing_init(p_row, p_col, nx, ny, nz, dttxi, dttyi, dttzi)
   if (dttxi == -1) then
      dttxi = 2
      dttxf = 9
   else
      dttxf = dttxi
   end if
   if (dttyi == -1) then
      dttyi = 2
      dttyf = 9
   else
      dttyf = dttyi
   end if
   if (dttzi == -1) then
      dttzi = 2
      dttzf = 9
   else
      dttzf = dttzi
   end if

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Define all combinations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do j = 1, 9
      do k = 1, 9
         do l = 1, 9
            DTT(:, j, k, l) = (/j - 1, k - 1, l - 1/)
         end do
      end do
   end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Test the dtt/idtt interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do j = dttxi, dttxf
      do k = dttyi, dttyf
         do l = dttzi, dttzf

            ! Init the FFT engine
            call dtt_engine%init(pencil=PHYSICAL_IN_X, &
                                 nx=nx, ny=ny, nz=nz, &
                                 opt_DTT=DTT(:, j, k, l))
            call dtt_engine%use_it()

            ! Get the decomp_info objects
            ph => decomp_2d_fft_get_ph()
            sp => ph

            ! Allocate memory and set to zero
            call alloc_x(in_r, ph, opt_global=.true.)
            in_r = 0._mytype
            call alloc_z(out_r, sp, opt_global=.true.)
            out_r = 0._mytype

            ! Local size
            st1 = ph%xst(1); en1 = ph%xen(1)
            st2 = ph%xst(2); en2 = ph%xen(2)
            st3 = ph%xst(3); en3 = ph%xen(3)

            ! Set input
            do kk = st3, en3
               do jj = st2, en2
                  do ii = st1, en1
                     in_r(ii, jj, kk) = real(ii, mytype) / real(nx, mytype) &
                                        * real(jj, mytype) / real(ny, mytype) &
                                        * real(kk, mytype) / real(nz, mytype)
                  end do
               end do
            end do

            ! Perform forward and backward transform once
            ! Forward, real input, real output
            call decomp_2d_dtt_3d_r2r(in_r, out_r, DECOMP_2D_FFT_FORWARD)
            ! Backward, real input, real output
            call decomp_2d_dtt_3d_r2r(out_r, in_r, DECOMP_2D_FFT_BACKWARD)

            ! Rescale
            if (DTT(1, j, k, l) == 0) then
               in_r = in_r / real(nx, mytype)
            else if (DTT(1, j, k, l) == 1 .or. DTT(1, j, k, l) == 4) then
               in_r = in_r / real(2 * nx - 2, mytype)
            else if (DTT(1, j, k, l) == 5) then
               in_r = in_r / real(2 * nx + 2 - 2 * dtt_engine%dtt(7), mytype)
            else
               in_r = in_r / real(2 * nx - 2 * dtt_engine%dtt(7), mytype)
            end if
            if (DTT(2, j, k, l) == 0) then
               in_r = in_r / real(ny, mytype)
            else if (DTT(2, j, k, l) == 1 .or. DTT(2, j, k, l) == 4) then
               in_r = in_r / real(2 * ny - 2, mytype)
            else if (DTT(2, j, k, l) == 5) then
               in_r = in_r / real(2 * ny + 2 - 2 * dtt_engine%dtt(8), mytype)
            else
               in_r = in_r / real(2 * ny - 2 * dtt_engine%dtt(8), mytype)
            end if
            if (DTT(3, j, k, l) == 0) then
               in_r = in_r / real(nz, mytype)
            else if (DTT(3, j, k, l) == 1 .or. DTT(3, j, k, l) == 4) then
               in_r = in_r / real(2 * nz - 2, mytype)
            else if (DTT(3, j, k, l) == 5) then
               in_r = in_r / real(2 * nz + 2 - 2 * dtt_engine%dtt(9), mytype)
            else
               in_r = in_r / real(2 * nz - 2 * dtt_engine%dtt(9), mytype)
            end if

            ! Check the error
            errl2 = 0._mytype
            errli = 0._mytype
            ! Local size
            st1 = dtt_engine%dtt(4)
            en1 = ph%xen(1) - dtt_engine%dtt(7) + dtt_engine%dtt(4) - 1
            if (ph%xst(2) == 1) then
               st2 = dtt_engine%dtt(5)
            else
               st2 = ph%xst(2)
            end if
            if (ph%xen(2) == ny) then
               en2 = ph%xen(2) - dtt_engine%dtt(8) + dtt_engine%dtt(5) - 1
            else
               en2 = ph%xen(2)
            end if
            if (ph%xst(3) == 1) then
               st3 = dtt_engine%dtt(6)
            else
               st3 = ph%xst(3)
            end if
            if (ph%xen(3) == nz) then
               en3 = ph%xen(3) - dtt_engine%dtt(9) + dtt_engine%dtt(6) - 1
            else
               en3 = ph%xen(3)
            end if
            do kk = st3, en3
               do jj = st2, en2
                  do ii = st1, en1
                     error = abs(in_r(ii, jj, kk) - real(ii, mytype) / real(nx, mytype) &
                                 * real(jj, mytype) / real(ny, mytype) &
                                 * real(kk, mytype) / real(nz, mytype))
                     errl2 = errl2 + error**2
                     errli = max(errli, error)
                  end do
               end do
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE, errl2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
            errl2 = errl2 / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
            call MPI_ALLREDUCE(MPI_IN_PLACE, errli, 1, real_type, MPI_MAX, MPI_COMM_WORLD, ierror)
            if (nrank == 0) write (*, '(A,ES14.8,1x,A,1x,ES14.8)') 'L2 and Linf error: ', &
                                                                    sqrt(errl2), errli

            if (errli > 100 * epsilon(errli)) then
               if (nrank == 0) then
                  write (*, *) "check case ", j, k, l, DTT(:, j, k, l)
                  write (*, *) dtt_engine%dtt
                  write (*, *) real(in_r(:, st2, st3) * nx * ny * nz, kind=real32)
               end if
               call decomp_2d_abort(__FILE__, __LINE__, 1, "DTT: Error above limit")
            end if

            ! Free memory and objects
            if (allocated(out_r)) deallocate (out_r)
            deallocate (in_r)
            nullify (sp)
            nullify (ph)

            ! Clear the FFT engine
            call dtt_engine%fin()

         end do
      end do
   end do

   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program dtt_x
