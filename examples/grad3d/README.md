# Halo test

List of the tests:
- [grad3d](grad3d.f90): Example to compute the gradient of a field. 

This example demonstrates the use 2DECOMP&FFT library to compute the gradient
of a field using an explicit second order finite difference scheme. 
The purpose is to show how to use the transpose operations to allow explicit calculation
of the gradient in all 3 directions. The results are written on a file and the function 
is periodic over the interval [0-1] 

What to input: The program takes max 5 inputs as: 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]

If the decomposition is imposed both (1) and (2) are necessary. 
If the resolution is imposed (1-5) are necessary.

What to expect: the output is the original function and the gradient in the 3 direction. 
                The program will also give the total error in L2 norm comparing with the anytical solution. 
