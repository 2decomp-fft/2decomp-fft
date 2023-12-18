# Test FFT for multiple grids

List of the tests:
- [fft_multigrid](fft_multigrid.f90): Test FFT on multiple grids using the dedicated object;
- [fft_multigrid_no_obj](fft_multigrid_no_obj.f90): Test FFT on multiple grids without using the dedicated object;
- [fft_multigrid_inplace](fft_multigrid_inplace.f90): Test FFT on one grid with various options (physical in X or Z, in-place transforms)


Both c2c and r2c/c2r transforms are tested.
The results should recover the input data up to machine accuracy
after a forward and a backward transform and appropriate normalisation.
The test automatically resize the problem depending on the number of MPI processes in use

What to input: The program takes max 6 inputs as : 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]
1. nt    [optional]

If the decomposition is imposed both (1) and (2) are necessary. 
If the resolution is imposed (1-5) are necessary

What to expect:
- The timing results 
- The error reported should be around machine accuracy (~ 10^-6 for single
  precision and 10^-15 for double)
- In case of the GENERIC FFT expect an increase in the order of the error
