# Initialization test

List of the tests:
- [Initialization test](init_test.f90): Test for initialization process. 

This example demonstrates the initialisation of DECOMP2D&FFT library and tests 
that the size of the mesh in every direction is as expected


What to input: The program takes max 5 inputs as: 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]

If the decomposition is imposed both (1) and (2) are necessary. 
If the resolution is imposed (1-5) are necessary

What to expect: Success/error message if initialization process has been complted or 
                errors have been encountered.
