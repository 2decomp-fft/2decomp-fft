!! SPDX-License-Identifier: BSD-3-Clause

!
! Small module to provide the reference solution
!
! The input coordinate (i, j, k) is the global one
!
module reference_data

   public

contains

   ! This subroutines provides the reference real data
   pure function reference_rdata(i, j, k, nx, ny, nz)

      !$acc routine seq
      use decomp_2d_constants, only: mytype

      implicit none

      integer, intent(in) :: i, j, k, nx, ny, nz
      real(mytype) :: reference_rdata

      reference_rdata = real(i, mytype) / real(nx, mytype) &
                        * real(j, mytype) / real(ny, mytype) &
                        * real(k, mytype) / real(nz, mytype)

   end function reference_rdata

   ! This subroutine provides the reference complex data
   pure function reference_cdata(i, j, k, nx, ny, nz)

      !$acc routine seq
      use decomp_2d_constants, only: mytype

      implicit none

      integer, intent(in) :: i, j, k, nx, ny, nz
      complex(mytype) :: reference_cdata

      reference_cdata = cmplx(reference_rdata(i, j, k, nx, ny, nz), &
                              -3 * reference_rdata(i, j, k, nx, ny, nz), &
                              mytype)

   end function reference_cdata

end module reference_data

!
! Small module to run c2c + inv. c2c and r2C + c2r
!
! Measure the error against reference solution and time the FFT
!
module mod_ntest

   public

contains

! Run a test for c2c transform several times
   subroutine ntest_c2c(nt)

      use decomp_2d
      use decomp_2d_constants
      use decomp_2d_fft
      use decomp_2d_mpi
      use MPI
      use reference_data
#if defined(_GPU)
      use cudafor
      use cufft
      use openacc
#endif

      implicit none

      integer, intent(in) :: nt

      ! Local variables
      integer :: nx, ny, nz, st1, en1, st2, en2, st3, en3
      integer :: i, j, k, itest, ierror, format
      type(decomp_info), pointer :: ph => null()
      complex(mytype), dimension(:, :, :), allocatable :: in, out
      double precision :: t1, t2, tt
      real(mytype) :: error

      ! Get the decomp_info object describing the arrays
      ph => decomp_2d_fft_get_ph()
      nx = ph%xsz(1); ny = ph%ysz(2); nz = ph%zsz(3)

      ! Get the format
      format = decomp_2d_fft_get_format()

      ! Allocate memory
      if (format == PHYSICAL_IN_X) then
         call alloc_x(in, ph, opt_global=.true.)
         call alloc_z(out, ph, opt_global=.true.)
      else
         call alloc_z(in, ph, opt_global=.true.)
         call alloc_x(out, ph, opt_global=.true.)
      end if

      ! Initialise input
      if (format == PHYSICAL_IN_X) then
         st1 = ph%xst(1); en1 = ph%xen(1)
         st2 = ph%xst(2); en2 = ph%xen(2)
         st3 = ph%xst(3); en3 = ph%xen(3)
      else
         st1 = ph%zst(1); en1 = ph%zen(1)
         st2 = ph%zst(2); en2 = ph%zen(2)
         st3 = ph%zst(3); en3 = ph%zen(3)
      end if
      do k = st3, en3
         do j = st2, en2
            do i = st1, en1
               in(i, j, k) = reference_cdata(i, j, k, nx, ny, nz)
            end do
         end do
      end do

      !$acc data copyin(in), copyout(out)

      ! First iteration is out of the loop with a dedicated timer
      t1 = MPI_WTIME()
      call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
      t1 = MPI_WTIME() - t1
      ! Reset the initial array
      !$acc kernels
      in = 0._mytype
      !$acc end kernels
      t2 = MPI_WTIME()
      call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
      t2 = MPI_WTIME() - t2
      ! Rescale input
      !$acc kernels
      in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
      !$acc end kernels
      ! MPI reduction of timers
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      if (nrank == 0) then
         write (*, *) '===== c2c interface ====='
         write (*, *) 'Grid ', nx, ny, nz
         if (format == PHYSICAL_IN_X) then
            write (*, *) 'Format PHYSICAL_IN_X'
         else
            write (*, *) 'Format PHYSICAL_IN_Z'
         end if
         if (decomp_2d_fft_get_inplace()) then
            write (*, *) 'In-place c2c'
         end if
         write (*, *) ''
         write (*, *) 'First iteration with dedicated timer'
         write (*, *) '     time (sec): ', t1 / dble(nproc), t2 / dble(nproc)
         write (*, *) ''
      end if

      ! Init the timer and run the tests
      t1 = 0.d0
      t2 = 0.d0
      do itest = 1, nt

         ! Forward FFT
         tt = MPI_WTIME()
         call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
         tt = MPI_WTIME() - tt
         t1 = t1 + tt

         ! Reset the initial array
         !$acc kernels
         in = 0._mytype
         !$acc end kernels

         ! Inverse FFT
         tt = MPI_WTIME()
         call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
         tt = MPI_WTIME() - tt
         t2 = t2 + tt

         ! Rescale input
         !$acc kernels
         in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
         !$acc end kernels

      end do

#if defined(_GPU)
      ierror = cudaDeviceSynchronize()
#endif

      ! MPI reduction of timers
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      t1 = t1 / dble(nproc) / dble(nt)
      t2 = t2 / dble(nproc) / dble(nt)

      ! Compute the sqrt of the L2 error
      error = 0._mytype
      if (format == PHYSICAL_IN_X) then
         st1 = ph%xst(1); en1 = ph%xen(1)
         st2 = ph%xst(2); en2 = ph%xen(2)
         st3 = ph%xst(3); en3 = ph%xen(3)
      else
         st1 = ph%zst(1); en1 = ph%zen(1)
         st2 = ph%zst(2); en2 = ph%zen(2)
         st3 = ph%zst(3); en3 = ph%zen(3)
      end if
      !$acc parallel loop default(present) reduction(+:error)
      do k = st3, en3
         do j = st2, en2
            do i = st1, en1
               error = error + (real(in(i, j, k), mytype) - reference_rdata(i, j, k, nx, ny, nz))**2 &
                       + (aimag(in(i, j, k)) + 3 * reference_rdata(i, j, k, nx, ny, nz))**2
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(MPI_IN_PLACE, error, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      error = sqrt(error) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      ! Print some stuff
      if (nrank == 0) then
         write (*, *) 'error / mesh point: ', error
         write (*, *) 'Avg time (sec): ', t1, t2
         write (*, *) ''
      end if

      !$acc end data

      ! Free memory, nullify pointers
      deallocate (in, out)
      nullify (ph)

      ! Abort if the error is too high
      if (error > epsilon(error) * 5 * nt) then
         call decomp_2d_abort(__FILE__, __LINE__, int(log10(error)), "c2c test")
      end if

   end subroutine ntest_c2c

! Run a test for r2c transform several times
   subroutine ntest_r2c(nt)

      use decomp_2d
      use decomp_2d_fft
      use decomp_2d_constants
      use decomp_2d_mpi
      use MPI
      use reference_data
#if defined(_GPU)
      use cudafor
      use cufft
      use openacc
#endif

      implicit none

      integer, intent(in) :: nt

      ! Local variables
      integer :: nx, ny, nz, st1, en1, st2, en2, st3, en3
      integer :: i, j, k, itest, ierror, format
      type(decomp_info), pointer :: ph => null(), sp => null()
      real(mytype), target, dimension(:), allocatable :: in1d
      real(mytype), pointer, contiguous, dimension(:, :, :) :: in
      complex(mytype), dimension(:, :, :), allocatable :: out
      double precision :: t1, t2, tt
      real(mytype) :: error

      ! Get the decomp_info object describing the arrays
      ph => decomp_2d_fft_get_ph()
      sp => decomp_2d_fft_get_sp()
      nx = ph%xsz(1); ny = ph%ysz(2); nz = ph%zsz(3)

      ! Get the format
      format = decomp_2d_fft_get_format()

      ! Allocate memory
      if (format == PHYSICAL_IN_X) then
         if (decomp_2d_fft_get_inplace_r2c() .or. decomp_2d_fft_get_inplace_c2r()) then
            allocate (in1d(max(product(ph%xsz), 2 * product(sp%xsz))))
            in(ph%xst(1):ph%xen(1), ph%xst(2):ph%xen(2), ph%xst(3):ph%xen(3)) => in1d
         else
            allocate (in(ph%xst(1):ph%xen(1), ph%xst(2):ph%xen(2), ph%xst(3):ph%xen(3)))
         end if
         call alloc_z(out, sp, opt_global=.true.)
      else
         if (decomp_2d_fft_get_inplace_r2c() .or. decomp_2d_fft_get_inplace_c2r()) then
            allocate (in1d(max(product(ph%zsz), 2 * product(sp%zsz))))
            in(ph%zst(1):ph%zen(1), ph%zst(2):ph%zen(2), ph%zst(3):ph%zen(3)) => in1d
         else
            allocate (in(ph%zst(1):ph%zen(1), ph%zst(2):ph%zen(2), ph%zst(3):ph%zen(3)))
         end if
         call alloc_x(out, sp, opt_global=.true.)
      end if

      ! Initialise input
      if (format == PHYSICAL_IN_X) then
         st1 = ph%xst(1); en1 = ph%xen(1)
         st2 = ph%xst(2); en2 = ph%xen(2)
         st3 = ph%xst(3); en3 = ph%xen(3)
      else
         st1 = ph%zst(1); en1 = ph%zen(1)
         st2 = ph%zst(2); en2 = ph%zen(2)
         st3 = ph%zst(3); en3 = ph%zen(3)
      end if
      do k = st3, en3
         do j = st2, en2
            do i = st1, en1
               in(i, j, k) = reference_rdata(i, j, k, nx, ny, nz)
            end do
         end do
      end do

      !$acc data copyin(in), copyout(out)

      ! First iteration is out of the loop with a dedicated timer
      t1 = MPI_WTIME()
      call decomp_2d_fft_3d(in, out)
      t1 = MPI_WTIME() - t1
      ! Reset the initial array
      !$acc kernels
      in = 0._mytype
      !$acc end kernels
      t2 = MPI_WTIME()
      call decomp_2d_fft_3d(out, in)
      t2 = MPI_WTIME() - t2
      ! Rescale input
      !$acc kernels
      in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
      !$acc end kernels
      ! MPI reduction of timers
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      if (nrank == 0) then
         write (*, *) '===== r2c interface ====='
         write (*, *) 'Grid ', nx, ny, nz
         if (format == PHYSICAL_IN_X) then
            write (*, *) 'Format PHYSICAL_IN_X'
         else
            write (*, *) 'Format PHYSICAL_IN_Z'
         end if
         if (decomp_2d_fft_get_inplace()) then
            write (*, *) 'In-place c2c'
            if (decomp_2d_fft_get_inplace_r2c()) write (*, *) 'In-place r2c'
            if (decomp_2d_fft_get_inplace_c2r()) write (*, *) 'In-place c2r'
         end if
         write (*, *) ''
         write (*, *) 'First iteration with dedicated timer'
         write (*, *) '     time (sec): ', t1 / dble(nproc), t2 / dble(nproc)
         write (*, *) ''
      end if

      ! Init the timer and run the tests
      t1 = 0.d0
      t2 = 0.d0
      do itest = 1, nt

         ! Forward FFT
         tt = MPI_WTIME()
         call decomp_2d_fft_3d(in, out)
         tt = MPI_WTIME() - tt
         t1 = t1 + tt

         ! Reset the initial array
         !$acc kernels
         in = 0._mytype
         !$acc end kernels

         ! Inverse FFT
         tt = MPI_WTIME()
         call decomp_2d_fft_3d(out, in)
         tt = MPI_WTIME() - tt
         t2 = t2 + tt

         ! Rescale input
         !$acc kernels
         in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
         !$acc end kernels

      end do

#if defined(_GPU)
      ierror = cudaDeviceSynchronize()
#endif

      ! MPI reduction of timers
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      t1 = t1 / dble(nproc) / dble(nt)
      t2 = t2 / dble(nproc) / dble(nt)

      ! Compute the sqrt of the L2 error
      error = 0._mytype
      if (format == PHYSICAL_IN_X) then
         st1 = ph%xst(1); en1 = ph%xen(1)
         st2 = ph%xst(2); en2 = ph%xen(2)
         st3 = ph%xst(3); en3 = ph%xen(3)
      else
         st1 = ph%zst(1); en1 = ph%zen(1)
         st2 = ph%zst(2); en2 = ph%zen(2)
         st3 = ph%zst(3); en3 = ph%zen(3)
      end if
      !$acc parallel loop default(present) reduction(+:error)
      do k = st3, en3
         do j = st2, en2
            do i = st1, en1
               error = error + (in(i, j, k) - reference_rdata(i, j, k, nx, ny, nz))**2
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(MPI_IN_PLACE, error, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      error = sqrt(error) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      ! Print some stuff
      if (nrank == 0) then
         write (*, *) 'error / mesh point: ', error
         write (*, *) 'Avg time (sec): ', t1, t2
         write (*, *) ''
      end if

      !$acc end data

      ! Free memory, nullify pointers
      if (allocated(in1d)) then
         deallocate (in1d)
      else
         deallocate (in)
      end if
      nullify (in)
      deallocate (out)
      nullify (sp)
      nullify (ph)

      ! Abort if the error is too high
      if (error > epsilon(error) * 5 * nt) then
         call decomp_2d_abort(__FILE__, __LINE__, int(log10(error)), "c2c test")
      end if

   end subroutine ntest_r2c

end module mod_ntest
