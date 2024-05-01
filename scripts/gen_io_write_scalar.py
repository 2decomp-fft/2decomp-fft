#!/usr/bin/env python3

# Name of the output file
outfile = "write_scalar.f90"

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
f.write("   interface decomp_2d_write_scalar\n")
for i in range(nformat):
    f.write("      module procedure write_scalar_"+ext[i]+"\n")

f.write("   end interface decomp_2d_write_scalar\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine write_scalar_"+ext[i]+"(io, n, var)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    f.write("      type(d2d_io_mpi), intent(INOUT) :: io\n")
    f.write("      integer, intent(IN) :: n\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:), intent(IN) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:), intent(IN) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:), intent(IN) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:), intent(IN) :: var\n")
    elif (i==4):
        f.write("      integer, contiguous, dimension(:), intent(IN) :: var\n")
    elif (i==5):
        f.write("      logical, contiguous, dimension(:), intent(IN) :: var\n")
    #
    f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"io_write_scalar\")\n")
    f.write("\n")
    #
    # Call the lower level IO subroutine
    #
    if (i==0 or i==1 or i==2 or i==3 or i==4 or i==5):
        f.write("      call write_scalar(io, n, &\n")
        if (i==0):
            f.write("                        freal=var)\n")
        elif (i==1):
            f.write("                        fcplx=var)\n")
        elif (i==2):
            f.write("                        dreal=var)\n")
        elif (i==3):
            f.write("                        dcplx=var)\n")
        elif (i==4):
            f.write("                        ints=var)\n")
        elif (i==5):
            f.write("                        logs=var)\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"io_write_scalar\")\n")
    f.write("\n")
    #
    # Footer
    #
    f.write("   end subroutine write_scalar_"+ext[i]+"\n")
    f.write("   !\n")

#
# Closing the file
#
f.close()
