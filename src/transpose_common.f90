!! SPDX-License-Identifier: BSD-3-Clause

module transpose_common

  use decomp_2d_constants, only: mytype

  implicit none

  private
  public transpose_real_fast
  
contains

  ! Fast implementation of transpose for real numbers avoiding communication
  subroutine transpose_real_fast(src, dst, nsize)

    implicit none

    real(mytype), dimension(:, :, :), intent(IN) :: src
    real(mytype), dimension(:, :, :), intent(OUT) :: dst
    integer, intent(in) :: nsize
#if defined(_GPU)
    integer :: istat
#endif

#if defined(_GPU)
    !$acc host_data use_device(src,dst)
    istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
    !$acc end host_data
    if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
    associate(foo => nsize)
    end associate
    dst = src
#endif

  end subroutine transpose_real_fast
  
end module transpose_common
