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

submodule (decomp_2d) d2d_log

  implicit none

  contains

  !
  ! Print some information about decomp_2d
  !
  module subroutine d2d_listing(given_io_unit)

    use iso_fortran_env, only : output_unit, compiler_version, compiler_options

    implicit none

    ! Argument
    integer, intent(in), optional :: given_io_unit

    ! Local variable
    integer :: io_unit
    integer :: version, subversion, ierror
#ifdef DEBUG
    character(len=64) :: fname
#endif

    ! Output log if needed
    if (decomp_log == D2D_LOG_NO) return
    if (decomp_log == D2D_LOG_STDOUT .and. nrank /= 0) return
    if (decomp_log == D2D_LOG_FILE .and. nrank /= 0) return

    ! If no IO unit provided, use stdout
    if (present(given_io_unit)) then
       io_unit = given_io_unit
    else
       io_unit = output_unit
    endif

    ! Header
    write (io_unit, *) '==========================================================='
    write (io_unit, *) '=================== Decomp2D - log ========================'
    write (io_unit, *) '==========================================================='

    ! Git hash if available
#if defined(VERSION)
    write (io_unit, *) 'Git version        : ', VERSION
#else
    write (io_unit, *) 'Git version        : unknown'
#endif

    ! Basic info
#ifdef DEBUG
    if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) &
       write (io_unit, *) 'I am mpi rank ', nrank
#endif
    write (io_unit, *) 'Total ranks ', nproc
    write (io_unit, *) 'Global data size : ', nx_global, ny_global, nz_global
    write (io_unit, *) 'p_row, p_col : ', dims(1), dims(2)
    write (io_unit, *) 'Periodicity : ', periodic_x, periodic_y, periodic_z
    write (io_unit, *) 'Number of bytes / float number : ', mytype_bytes
    write (io_unit, *) '==========================================================='

    ! Show detected flags, compiler options, version of the MPI library
    write (io_unit, *) 'Compile flags detected :'
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
    write (io_unit, *) 'Numerical precision: Double, saving in single'
#else
    write (io_unit, *) 'Numerical precision: Double'
#endif
#else
    write (io_unit, *) 'Numerical precision: Single'
#endif
    write (io_unit, *) 'Compiled with ', compiler_version()
    write (io_unit, *) 'Compiler options : ', compiler_options()
    call MPI_Get_version(version, subversion, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_Get_version")
    write (io_unit, '(" Version of the MPI library : ",I0,".",I0)') version, subversion
#ifdef DEBUG
    write (io_unit, *) 'Compile flag DEBUG detected'
    write (io_unit, *) '   debug level : ', decomp_debug
#endif
#ifdef PROFILER
    write (io_unit, *) 'Compile flag PROFILER detected'
#endif
#ifdef SHM
    write (io_unit, *) 'Compile flag SHM detected'
#endif
#ifdef EVEN
    write (io_unit, *) 'Compile flag EVEN detected'
#endif
#ifdef OCC
    write (io_unit, *) 'Compile flag OCC detected'
#endif
#ifdef OVERWRITE
    write (io_unit, *) 'Compile flag OVERWRITE detected'
#endif
#ifdef HALO_DEBUG
    write (io_unit, *) 'Compile flag HALO_DEBUG detected'
#endif
#ifdef SHM_DEBUG
    write (io_unit, *) 'Compile flag SHM_DEBUG detected'
#endif
#ifdef _GPU
    write (io_unit, *) 'Compile flag _GPU detected'
#endif
#ifdef _NCCL
    write (io_unit, *) 'Compile flag _NCCL detected'
#endif
    write (io_unit, *) '==========================================================='
    write (io_unit, *) 'Profiler id : ', decomp_profiler
#ifdef PROFILER
    call decomp_profiler_log(io_unit)
    write(io_unit, *) "   Profiling transpose : ", decomp_profiler_transpose
    write(io_unit, *) "   Profiling IO : ", decomp_profiler_io
    write(io_unit, *) "   Profiling FFT : ", decomp_profiler_fft
    write(io_unit, *) "   Profiling decomp_2d : ", decomp_profiler_d2d
#endif
    write (io_unit, *) '==========================================================='
    ! Info about each decomp_info object
    call decomp_info_print(decomp_main, io_unit, "decomp_main")
    call decomp_info_print(phG, io_unit, "phG")
    call decomp_info_print(ph1, io_unit, "ph1")
    call decomp_info_print(ph2, io_unit, "ph2")
    call decomp_info_print(ph3, io_unit, "ph3")
    call decomp_info_print(ph4, io_unit, "ph4")
#ifdef SHM_DEBUG
    write (io_unit, *) '==========================================================='
    call print_smp(io_unit)
#endif
    write (io_unit, *) '==========================================================='
    write (io_unit, *) '==========================================================='
#ifdef DEBUG
    !
    ! In DEBUG mode, rank 0 will also print environment variables
    !
    ! At high debug level, all ranks will print env. variables
    !
    ! The system call, if writing to a file, is not blocking if supported
    !
    if (nrank == 0 .or. decomp_debug >= D2D_DEBUG_LEVEL_INFO) then
       write (io_unit, *) '============== Environment variables ======================'
       write (io_unit, *) '==========================================================='
       write (io_unit, *) '==========================================================='
       if (io_unit == output_unit ) then
          call execute_command_line("env", wait = .true.)
       else
          inquire(unit = io_unit, name = fname, iostat = ierror)
          if (ierror /= 0) call decomp_2d_abort(__FILE__, &
                                                __LINE__, &
                                                ierror, &
                                                "No name for the log file")
          call execute_command_line("env >> "//trim(fname), wait = .false.)
       endif
    endif
#endif

  end subroutine d2d_listing

  !
  ! Print some information about given decomp_info object
  !
  module subroutine decomp_info_print(d2d, io_unit, d2dname)

    implicit none

    ! Arguments
    type(decomp_info), intent(in) :: d2d
    integer, intent(in) :: io_unit
    character(len=*), intent(in) :: d2dname

    ! Nothing to print if not initialized
    if (.not.allocated(d2d%x1dist)) then
      write (io_unit, *) 'Uninitialized decomp_info ', d2dname
      return
    endif

    !
    ! If DEBUG mode, print everything
    ! Otherwise, print only global size
    !
    write (io_unit, *) 'Decomp_info : ', d2dname
    write (io_unit, *) '   Global size : ', d2d%xsz(1), d2d%ysz(2), d2d%zsz(3)
#ifdef DEBUG
    write (io_unit, *) '   xsz, xst, xen : ', d2d%xsz, d2d%xst, d2d%xen
    write (io_unit, *) '   ysz, yst, yen : ', d2d%ysz, d2d%yst, d2d%yen
    write (io_unit, *) '   zsz, zst, zen : ', d2d%zsz, d2d%zst, d2d%zen
    write (io_unit, *) '   x1dist : ', d2d%x1dist
    write (io_unit, *) '   y1dist : ', d2d%y1dist
    write (io_unit, *) '   y2dist : ', d2d%y2dist
    write (io_unit, *) '   z2dist : ', d2d%z2dist
    write (io_unit, *) '   x1cnts : ', d2d%x1cnts
    write (io_unit, *) '   y1cnts : ', d2d%y1cnts
    write (io_unit, *) '   y2cnts : ', d2d%y2cnts
    write (io_unit, *) '   z2cnts : ', d2d%z2cnts
    write (io_unit, *) '   x1disp : ', d2d%x1disp
    write (io_unit, *) '   y1disp : ', d2d%y1disp
    write (io_unit, *) '   y2disp : ', d2d%y2disp
    write (io_unit, *) '   z2disp : ', d2d%z2disp
    write (io_unit, *) '   x1count : ', d2d%x1count
    write (io_unit, *) '   y1count : ', d2d%y1count
    write (io_unit, *) '   y2count : ', d2d%y2count
    write (io_unit, *) '   z2count : ', d2d%z2count
    write (io_unit, *) '   even : ', d2d%even
#ifdef SHM
    write (io_unit, *) '   listing of the SHM part is not yet implemented'
#endif
#endif

  end subroutine decomp_info_print

#ifdef SHM_DEBUG

  ! For debugging, print the shared-memory structure
  module subroutine print_smp(io_unit)

    implicit none
    
    ! Argument
    integer, intent(in) :: io_unit

    ! print out shared-memory information
    write(io_unit,*)'I am mpi rank ', nrank, 'Total ranks ', nproc
    write(io_unit,*)' '
    write(io_unit,*)'Global data size:'
    write(io_unit,*)'nx*ny*nz', nx,ny,nz
    write(io_unit,*)' '
    write(io_unit,*)'2D processor grid:'
    write(io_unit,*)'p_row*p_col:', dims(1), dims(2)
    write(io_unit,*)' '
    write(io_unit,*)'Portion of global data held locally:'
    write(io_unit,*)'xsize:',xsize
    write(io_unit,*)'ysize:',ysize
    write(io_unit,*)'zsize:',zsize
    write(io_unit,*)' '
    write(io_unit,*)'How pensils are to be divided and sent in alltoallv:'
    write(io_unit,*)'x1dist:',decomp_main%x1dist
    write(io_unit,*)'y1dist:',decomp_main%y1dist
    write(io_unit,*)'y2dist:',decomp_main%y2dist
    write(io_unit,*)'z2dist:',decomp_main%z2dist
    write(io_unit,*)' '
    write(io_unit,*)'######Shared buffer set up after this point######'
    write(io_unit,*)' '
    write(io_unit,*) 'col communicator detais:'
    call print_smp_info(decomp_main%COL_INFO, io_unit)
    write(io_unit,*)' '
    write(io_unit,*) 'row communicator detais:'
    call print_smp_info(decomp_main%ROW_INFO; io_unit)
    write(io_unit,*)' '
    write(io_unit,*)'Buffer count and displacement of per-core buffers'
    write(io_unit,*)'x1cnts:',decomp_main%x1cnts
    write(io_unit,*)'y1cnts:',decomp_main%y1cnts
    write(io_unit,*)'y2cnts:',decomp_main%y2cnts
    write(io_unit,*)'z2cnts:',decomp_main%z2cnts
    write(io_unit,*)'x1disp:',decomp_main%x1disp
    write(io_unit,*)'y1disp:',decomp_main%y1disp
    write(io_unit,*)'y2disp:',decomp_main%y2disp
    write(io_unit,*)'z2disp:',decomp_main%z2disp
    write(io_unit,*)' '
    write(io_unit,*)'Buffer count and displacement of shared buffers'
    write(io_unit,*)'x1cnts:',decomp_main%x1cnts_s
    write(io_unit,*)'y1cnts:',decomp_main%y1cnts_s
    write(io_unit,*)'y2cnts:',decomp_main%y2cnts_s
    write(io_unit,*)'z2cnts:',decomp_main%z2cnts_s
    write(io_unit,*)'x1disp:',decomp_main%x1disp_s
    write(io_unit,*)'y1disp:',decomp_main%y1disp_s
    write(io_unit,*)'y2disp:',decomp_main%y2disp_s
    write(io_unit,*)'z2disp:',decomp_main%z2disp_s

  end subroutine print_smp

  ! For debugging, print the shared-memory structure
  module subroutine print_smp_info(s, io_unit)

    implicit none

    ! Argument
    TYPE(SMP_INFO), intent(in) :: s
    integer, intent(in) :: io_unit

    write(io_unit,*) 'size of current communicator:', s%NCPU
    write(io_unit,*) 'rank in current communicator:', s%NODE_ME
    write(io_unit,*) 'number of SMP-nodes in this communicator:', s%NSMP
    write(io_unit,*) 'SMP-node id (1 ~ NSMP):', s%SMP_ME
    write(io_unit,*) 'NCORE - number of cores on this SMP-node', s%NCORE
    write(io_unit,*) 'core id (1 ~ NCORE):', s%CORE_ME
    write(io_unit,*) 'maximum no. cores on any SMP-node:', s%MAXCORE
    write(io_unit,*) 'size of SMP shared memory SND buffer:', s%N_SND
    write(io_unit,*) 'size of SMP shared memory RCV buffer:', s%N_RCV

  end subroutine print_smp_info
#endif

end submodule d2d_log
