#!/usr/bin/env python3

# Name of the output file
outfile = "adios_read_plane.f90"

# Number of variables
nformat = 4
# Suffix for the subroutine name
ext = ["freal", "fcplx", "dreal", "dcplx"]


#
# Open the file
#
f = open(outfile, "w")

#
# Write the interface
#
f.write("\n")
f.write("   interface decomp_2d_adios_read_plane\n")
for i in range(nformat):
    f.write("      module procedure read_plane_"+ext[i]+"\n")

f.write("   end interface decomp_2d_adios_read_plane\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine read_plane_"+ext[i]+"(io, var, varname, &\n")
    f.write("                               opt_family, &\n")
    f.write("                               opt_reduce_prec)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    f.write("      type(d2d_io_adios), intent(inout) :: io\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    f.write("      character(len=*), intent(in) :: varname\n")
    f.write("      type(d2d_io_family), intent(inout), optional :: opt_family\n")
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("\n")
    #
    # Local variables
    #
    if (i==2 or i==3):
        f.write("      ! Local variables\n")
        f.write("      logical :: reduce\n")
        if (i==2):
            f.write("      real(real32), allocatable, dimension(:, :, :) :: tmp\n")
        elif (i==3):
            f.write("      complex(real32), allocatable, dimension(:, :, :) :: tmp\n")
        f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"adios_read_plane\")\n")
    f.write("\n")
    #
    # Deal with optional arguments
    #
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
    if (i==0 or i==1):
        f.write("      call adios_read(io, varname, &\n")
        f.write("                      opt_family=opt_family, &\n")
        if (i==0):
            f.write("                      freal=var)\n")
        elif (i==1):
            f.write("                      fcplx=var)\n")
    #
    if (i==2 or i==3):
        f.write("      if (reduce) then\n")
        f.write("\n")
        f.write("         allocate (tmp(size(var, 1), &\n")
        f.write("                       size(var, 2), &\n")
        f.write("                       size(var, 3)))\n")
        f.write("         call adios_read(io, varname, &\n")
        f.write("                         opt_family=opt_family, &\n")
        if (i==2):
            f.write("                         freal=tmp)\n")
            f.write("         var = real(tmp, kind=real64)\n")
        elif (i==3):
            f.write("                         fcplx=tmp)\n")
            f.write("         var = cmplx(tmp, kind=real64)\n")
        f.write("         deallocate (tmp)\n")
        f.write("      else\n")
        f.write("         call adios_read(io, varname, &\n")
        f.write("                         opt_family=opt_family, &\n")
        if (i==2):
            f.write("                         dreal=var)\n")
        elif (i==3):
            f.write("                         dcplx=var)\n")
        f.write("      end if\n")
    #
    # Cleaning + unused variables
    #
    if (i==0 or i==1):
        f.write("      associate (p => opt_reduce_prec)\n")
        f.write("      end associate\n")
        f.write("\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"adios_read_plane\")\n")
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
