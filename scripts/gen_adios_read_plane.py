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
    f.write("   subroutine read_plane_"+ext[i]+"(io, var, varname, ipencil, &\n")
    f.write("                               opt_decomp, &\n")
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
    f.write("      integer, intent(in) :: ipencil\n")
    f.write("      class(info), intent(in), optional :: opt_decomp\n")
    f.write("      type(d2d_io_family), intent(inout), optional :: opt_family\n")
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("\n")
    #
    # Local variables
    #
    f.write("      ! Local variable(s)\n")
    f.write("      integer, dimension(3) :: sel_start, sel_count\n")
    if (i==2 or i==3):
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
    f.write("      if (present(opt_decomp)) then\n")
    f.write("         if (ipencil == 1) then\n")
    f.write("            sel_start = opt_decomp%xst\n")
    f.write("            sel_count = opt_decomp%xsz\n")
    f.write("         else if (ipencil == 2) then\n")
    f.write("            sel_start = opt_decomp%yst\n")                                                    
    f.write("            sel_count = opt_decomp%ysz\n")
    f.write("         else if (ipencil == 3) then\n")
    f.write("            sel_start = opt_decomp%zst\n")
    f.write("            sel_count = opt_decomp%zsz\n")
    f.write("         else\n")
    f.write("            call decomp_2d_abort(__FILE__, __LINE__, ipencil, \"Invalid value\")\n")
    f.write("         end if\n")
    f.write("      else\n")
    f.write("         if (ipencil == 1) then\n")
    f.write("            sel_start = decomp_main%xst\n")          
    f.write("            sel_count = decomp_main%xsz\n")
    f.write("         else if (ipencil == 2) then\n")
    f.write("            sel_start = decomp_main%yst\n")
    f.write("            sel_count = decomp_main%ysz\n")
    f.write("         else if (ipencil == 3) then\n")
    f.write("            sel_start = decomp_main%zst\n")
    f.write("            sel_count = decomp_main%zsz\n")
    f.write("         else\n")
    f.write("            call decomp_2d_abort(__FILE__, __LINE__, ipencil, \"Invalid value\")\n")
    f.write("         end if\n")
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
    if (i==0 or i==1):
        f.write("      call adios_read(io, varname, sel_start, sel_count, &\n")
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
        f.write("         call adios_read(io, varname, sel_start, sel_count, &\n")
        f.write("                         opt_family=opt_family, &\n")
        if (i==2):
            f.write("                         freal=tmp)\n")
            f.write("         var = real(tmp, kind=real64)\n")
        elif (i==3):
            f.write("                         fcplx=tmp)\n")
            f.write("         var = cmplx(tmp, kind=real64)\n")
        f.write("         deallocate (tmp)\n")
        f.write("      else\n")
        f.write("         call adios_read(io, varname, sel_start, sel_count, &\n")
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
