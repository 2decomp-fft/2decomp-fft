1. You want to contribute but have no idea how ? Please refer to the [Get started](#get-started) section.
2. You have identified a bug ? Please refer to the [Bug](#bug) section.
3. You have improved a part of the code or have developed a new functionality ? Please refer to the [Advanced](#advanced-contribution) section.

## Get started

## Bug

It appears that you have identified a bug in the 2decomp library.
If you are not sure this is really a bug in the library, you should go to the [discussions](https://github.com/2decomp-fft/2decomp-fft/discussions) section and open a new discussion.
Otherwise, follow the steps below.

Firstly, try to reproduce the error with a debug build of the library, a small problem size and a small number of MPI ranks.
It makes bug-hunting much easier.
Unfortunately, it is not always possible.
At least, please try to reproduce the bug on another machine with another compiler.

Secondly, if you have modified the source code of the 2decomp library, you must reproduce the bug without the modifications in 2decomp.
The development team will only provide support for sections of code available in the present repository.

Thirdly, you must provide a minimal working example.
The program using 2decomp and exposing the bug should be relatively small.
The development team will not provide support if the program exposing the bug is very long.
The programs available in the examples section are a good starting point for a minimal working example.

At this stage, you probably did your best to simplify the problem at hand.
Open a issue and select the bug report template.
Provide a meaningful title and do your best to complete all the sections of the template and provide the version of the library, the version of the compiler, the version of the MPI / FFT library, ...
If you think you have a fix for the bug, please expose it inside the issue.
It is recommended to wait for feedback before opening a pull-request.

## Advanced contribution
