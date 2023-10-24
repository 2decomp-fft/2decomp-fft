# Test FFT for Z pencil decomposition 

List of the tests:
- [fft_c2c_z](fft_c2c_z.f90): Test Complex to Complex FFT transform; 
- [fft_r2c_x](fft_r2c_z.f90): Test Real to Complex FFT transform;


These programs can be used to test the fft trasform using Z as starting physical dimension. 
Both c2c (fft_c2c_z) and r2c/c2r (fft_r2c_z) transforms are tested.
The results should recover the input data up to machine accuracy
after a forward and a backward transform and proper normalisation.
The test automatically resize the problem depending on the number of MPI processes in use

What to input: The program takes max 6 inputs as : 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]
1. nt    [optional]

In case the decomposition is imposed both (1) and (2) are necessary. 
If also the resolution is imposed (1-5) are necessary

What to expect:
- The timing results 
- The error reported should be around machine accuracy (~ 10^-6 for single
  precision and 10^-15 for double)
- In case of the GENERIC FFT expect an increase in the order of the error
