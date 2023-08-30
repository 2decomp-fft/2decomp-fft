io_bench, io_plane_test, io_read, io_test, io_var_test
---------

These examples demostrate the use of the different I/O features of DECOMP2D&FFT
using both MPI I/O and ADIOS2. The test performed are the following: 
- io_test: write of the files u1.dat, u2.dat and u3.dat allocated in the X, Y and Z pencil decomposition respectively            and the write of all 3 files in a combined single file (used for checkpointing)
           The program checks that the written files are the excpeted.
- io_read: read the files u1.dat, u2.dat and u3.dat written by io_test and check that they are the expected one. 
- io_var_test: test the writing of writing real and complex data for scalar and 3D arrays 
               and for different resolutions 
- io_plane_test: test the writing of 2D planes in different direction  
- io_bench: timing the writing on disk of a 3D array

What to input: The program takes max 6 inputs as : 
               1) p_row [optional]
	       2) p_col [optional] 
	       3) nx    [optional]
	       4) ny    [optional]
	       5) nz    [optional]
	       In case the decomposition is imposed both (1) and (2) are necessary. 
	       If also the resolution is imposed (1-5) are necessary

What to expect: 
- All programs print out a success message otherwise an error message with location of the error is reported
- io_bench gives also the timing of the writing on disk
