# Test and timing of the transpose operations 

List of the tests:
- [Test transpose operations](test2d.f90): Test the transpose operation for a 3D real array; 
- [Time transpose real](timing2d_real.f90): Time to solution for the transpose operation for a 3D real array; 
- [Time transpose complex](timing2d_complex.f90): Time to solution for the transpose operation for a 3D complex array.

The program test2d is to validate the 2D pencil decomposition library. It transposes
a set of data into different storage formats so that all communication routines
are tested.
The programs timing2d_real and timing2d_complex perform the same tests several times 
and are used to report timing and performances for real and complex transforms. 

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
- For test2d the output is a success message or an error message with the direction of swapping;
- For the timing, beside the error/success message timing for all transpostions
  (X->Y, Y->Z, Z->Y, Y->X) together with the sum are reported. 
