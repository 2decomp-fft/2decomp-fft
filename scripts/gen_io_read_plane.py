#!/usr/bin/env python3

# Name of the output file
outfile = "read_plane.f90"

# Number of variables
nformat = 6
# Suffix for the subroutine name
ext = ["freal", "fcplx", "dreal", "dcplx", "ints", "logs"]


#
# Open the file
#
f = open(outfile, "w")

#
# Write the interface
#
f.write("\n")
f.write("   interface decomp_2d_read_plane\n")
for i in range(nformat):
    f.write("      module procedure read_plane_"+ext[i]+"\n")

f.write("   end interface decomp_2d_read_plane\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine read_plane_"+ext[i]+"(ipencil, var, varname, nplanes, &\n")
    f.write("                               opt_dirname, &\n")
    f.write("                               opt_mpi_file_open_info, &\n")
    f.write("                               opt_mpi_file_set_view_info, &\n")
    f.write("                               opt_reduce_prec, &\n")
    f.write("                               opt_decomp, &\n")
    f.write("                               opt_nb_req, &\n")
    f.write("                               opt_io)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    f.write("      integer, intent(IN) :: ipencil\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==4):
        f.write("      integer, contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==5):
        f.write("      logical, contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    f.write("      character(len=*), intent(in) :: varname\n")
    f.write("      integer, intent(in) :: nplanes\n")
    f.write("      character(len=*), intent(in), optional :: opt_dirname\n")
    f.write("      integer, intent(in), optional :: opt_mpi_file_open_info\n")
    f.write("      integer, intent(in), optional :: opt_mpi_file_set_view_info\n")
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp\n")
    f.write("      integer, intent(inout), optional :: opt_nb_req\n")
    f.write("      type(d2d_io_mpi), intent(inout), optional :: opt_io")
    f.write("\n")
    #
    # Local variables
    #
    f.write("      ! Local variables\n")
    if (i==2 or i==3):
        f.write("      logical :: reduce\n")
    f.write("      TYPE(DECOMP_INFO), pointer :: decomp\n")
    if (i==2):
        f.write("      real(real32), allocatable, dimension(:, :, :) :: tmp\n")
    elif (i==3):
        f.write("      complex(real32), allocatable, dimension(:, :, :) :: tmp\n")
    f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"io_read_plane\")\n")
    f.write("\n")
    #
    # Deal with optional arguments
    #
    f.write("      ! Use the provided decomp_info or the default one\n")
    f.write("      if (present(opt_decomp)) then\n")
    f.write("         decomp => opt_decomp\n")
    f.write("      else\n")
    f.write("         decomp => decomp_main\n")
    f.write("      end if\n")
    f.write("\n")
    if (i==2 or i==3):
        f.write("      ! One can write to single precision using opt_reduce_prec\n")
        f.write("      if (present(opt_reduce_prec)) then\n")
        f.write("         reduce = opt_reduce_prec\n")
        f.write("      else\n")
        f.write("         reduce = default_opt_reduce_prec\n")
        f.write("      end if\n")
        f.write("\n")
    #
    # Call the lower level IO subroutine
    #
    if (i==0 or i==1 or i==4 or i==5):
        f.write("      call read_plane(ipencil, varname, decomp, nplanes, &\n")
        f.write("                      opt_dirname=opt_dirname, &\n")
        f.write("                      opt_mpi_file_open_info=opt_mpi_file_open_info, &\n")
        f.write("                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &\n")
        f.write("                      opt_nb_req=opt_nb_req, &\n")
        f.write("                      opt_io=opt_io, &\n")
        if (i==0):
            f.write("                      freal=var)\n")
        elif (i==1):
            f.write("                      fcplx=var)\n")
        elif (i==4):
            f.write("                      ints=var)\n")
        elif (i==5):
            f.write("                      logs=var)\n")
    #
    if (i==2 or i==3):
        f.write("      if (reduce) then\n")
        f.write("\n")
        f.write("         allocate (tmp(size(var, 1), &\n")
        f.write("                       size(var, 2), &\n")
        f.write("                       size(var, 3)))\n")
        f.write("         call read_plane(ipencil, varname, decomp, nplanes, &\n")
        f.write("                         opt_dirname=opt_dirname, &\n")
        f.write("                         opt_mpi_file_open_info=opt_mpi_file_open_info, &\n")
        f.write("                         opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &\n")
        if (i==2):
            f.write("                         freal=tmp)\n")
            f.write("         var = real(tmp, kind=real64)\n")
        elif (i==3):
            f.write("                         fcplx=tmp)\n")
            f.write("         var = cmplx(tmp, kind=real64)\n")
        f.write("         deallocate (tmp)\n")
        f.write("      else\n")
        f.write("         call read_plane(ipencil, varname, decomp, nplanes, &\n")
        f.write("                         opt_dirname=opt_dirname, &\n")
        f.write("                         opt_mpi_file_open_info=opt_mpi_file_open_info, &\n")
        f.write("                         opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &\n")
        f.write("                         opt_nb_req=opt_nb_req, &\n")
        f.write("                         opt_io=opt_io, &\n")
        if (i==2):
            f.write("                         dreal=var)\n")
        elif (i==3):
            f.write("                         dcplx=var)\n")
        f.write("      end if\n")
    #
    # Cleaning + unused variables
    #
    f.write("\n")
    f.write("      nullify (decomp)\n")
    f.write("\n")
    if (i==0 or i==1 or i==4 or i==5):
        f.write("      associate (p => opt_reduce_prec)\n")
        f.write("      end associate\n")
        f.write("\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"io_read_plane\")\n")
    f.write("\n")
    #
    # Footer
    #
    f.write("   end subroutine read_plane_"+ext[i]+"\n")
    f.write("   !\n")

#
# Closing the file
#
f.close()
