# IO tests

These examples demostrate the use of the different I/O features of DECOMP2D&FFT
using both MPI I/O and ADIOS2. The tests performed are the following: 
- [Write Test](io_test.f90): write the files u1.dat, u2.dat and u3.dat defined in X, Y and Z pencils respectively
                             and write these files in a combined single file (used for checkpointing)
                             The program also checks that the written files are correct;
- [Read Test](io_read.f90): read files u1.dat, u2.dat and u3.dat which were written by io_test and check that they are the expected ones;
- [Write Real/Complex variables](io_var_test.f90): test the writing of real and complex data for scalar and 3D arrays 
                                                   and for different resolutions;
- [Write 2D planes files](io_plane_test.f90): test the writing of 2D planes in different direction;  
- [Timing IO](io_bench.f90): timing the writing on disk of a 3D array.


What to input: The program takes max 5 inputs as : 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]

If the decomposition is imposed both (1) and (2) are necessary. 
If the resolution is imposed (1-5) are necessary.

What to expect: 
- All programs print out a success message otherwise an error message with location of the error;
- io_bench gives also the timing of the writing on disk.
