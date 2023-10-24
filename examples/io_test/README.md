# IO tests

These examples demostrate the use of the different I/O features of DECOMP2D&FFT
using both MPI I/O and ADIOS2. The test performed are the following: 
- [Write Test](io_test.f90): write of the files u1.dat, u2.dat and u3.dat allocated in the X, Y and Z pencil decomposition respectively            and the write of all 3 files in a combined single file (used for checkpointing)
           The program checks that the written files are the excpeted.
- [Read Test](io_read.f90): read the files u1.dat, u2.dat and u3.dat written by io_test and check that they are the expected one. 
- [Write Real/Complex variables](io_var_test.f90): test the writing of real and complex data for scalar and 3D arrays 
               and for different resolutions 
- [Write 2D planes files](io_plane_test.f90): test the writing of 2D planes in different direction  
- [Timing IO](io_bench.f90): timing the writing on disk of a 3D array


What to input: The program takes max 5 inputs as : 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]

In case the decomposition is imposed both (1) and (2) are necessary. 
If also the resolution is imposed (1-5) are necessary

What to expect: 
- All programs print out a success message otherwise an error message with location of the error is reported
- io_bench gives also the timing of the writing on disk
