!! SPDX-License-Identifier: BSD-3-Clause

submodule(decomp_2d_fft) d2d_fft_log

   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d

   implicit none

contains

   !
   ! Log subroutine for the decomp_2d_fft module
   !
   module subroutine decomp_2d_fft_log(backend)

      implicit none

      ! Argument
      character(len=*), intent(in) :: backend

      ! Local variable
      integer :: io_unit

      if (d2d_log_is_active()) then
         io_unit = d2d_log_get_unit()
         write (io_unit, *) ''
         write (io_unit, *) '***** Using the '//trim(backend)//' FFT engine *****'
         write (io_unit, *) ''
         write (io_unit, *) 'Id of the backend : ', D2D_FFT_BACKEND
         if (format == PHYSICAL_IN_X) then
            write (io_unit, *) 'Format : Physical in x'
         else
            write (io_unit, *) 'Format : Physical in z'
         end if
         write (io_unit, *) ''
         if (nx_fft == nx_global .and. ny_fft == ny_global .and. nz_fft == nz_global) then
            write (io_unit, *) 'decomp_info object ph is a pointer to decomp_main'
         end if
         call decomp_info_print(decomp_2d_fft_get_ph(), io_unit, "ph")
         call decomp_info_print(decomp_2d_fft_get_sp(), io_unit, "sp")
         write (io_unit, *) ''
         if (inplace) then
            if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_GENERIC) then
               write (io_unit, *) 'OVERWRITE is supported but transforms are not performed in-place'
               write (io_unit, *) ''
            else if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_CUFFT .or. &
                     D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3 .or. &
                     D2D_FFT_BACKEND == D2D_FFT_BACKEND_MKL) then
               write (io_unit, *) 'OVERWRITE is supported but in-place FFT is limited to complex FFT'
               write (io_unit, *) ''
            else if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3_F03) then
               write (io_unit, *) 'OVERWRITE is supported for all transforms'
               if (inplace_r2c) write (io_unit, *) '   in-place r2c transform'
               if (inplace_c2r) write (io_unit, *) '   in-place c2r transform'
               write (io_unit, *) ''
            end if
         end if
         if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3) then
            call decomp_2d_warning(D2D_FFT_BACKEND, "FFTW3 backend will be removed")
         end if
         if (skip_x_c2c) write (io_unit, *) 'Skip X c2c transforms'
         if (skip_y_c2c) write (io_unit, *) 'Skip Y c2c transforms'
         if (skip_z_c2c) write (io_unit, *) 'Skip Z c2c transforms'
         call d2d_log_close_unit(io_unit)
      end if

   end subroutine decomp_2d_fft_log

end submodule d2d_fft_log
