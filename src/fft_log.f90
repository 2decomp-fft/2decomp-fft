!=======================================================================
! This is part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
! Copyright (C) 2021               the University of Edinburgh (UoE)
!
!=======================================================================

submodule(decomp_2d_fft) d2d_fft_log

   use decomp_2d_constants
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

      if ((decomp_log == D2D_LOG_STDOUT .and. nrank == 0) .or. &
          (decomp_log == D2D_LOG_TOFILE .and. nrank == 0) .or. &
          (decomp_log == D2D_LOG_TOFILE_FULL)) then
         io_unit = d2d_listing_get_unit()
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
         call decomp_info_print(ph, io_unit, "ph")
         call decomp_info_print(sp, io_unit, "sp")
         write (io_unit, *) ''
#ifdef OVERWRITE
         if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_GENERIC .or. &
             D2D_FFT_BACKEND == D2D_FFT_BACKEND_CUFFT) then
            call decomp_2d_warning("Selected FFT backend does not support overwrite")
         else if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3 .or. &
                  D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3_F03) then
            write (io_unit, *) 'OVERWRITE support limited to c2C and c2r transforms'
            write (io_unit, *) ''
         else if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_MKL) then
            write (io_unit, *) 'OVERWRITE support limited to c2c and c2r transforms'
            write (io_unit, *) ''
         endif
#else
         if (D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3 .or. &
             D2D_FFT_BACKEND == D2D_FFT_BACKEND_FFTW3_F03) then
            write (io_unit, *) 'Warning : c2c and c2r transforms do not preserve their complex input'
            write (io_unit, *) ''
         endif
#endif
         call d2d_listing_close_unit(io_unit)
      end if

   end subroutine decomp_2d_fft_log

end submodule d2d_fft_log
