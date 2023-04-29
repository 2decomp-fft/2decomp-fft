1. You want to contribute but have no idea how ? Please refer to the [Get started](#get-started) section.
2. You have identified a bug ? Please refer to the [Bug](#bug) section.
3. You have improved a part of the code or have developed a new functionality ? Please refer to the [advanced](#advanced-contribution) section.

## Get started

The recommended strategy to contribute is to start with a [discussion](https://github.com/2decomp-fft/2decomp-fft/discussions) or to pick a existing issue.
To modify or experiment with the code, fork the 2decomp github repository and commit changes in a dedicated branch of your fork.
When the modification is ready for review, one can open a pull request as exposed in the [advanced](#advanced-contribution) section below.

## Bug

It appears that you have identified a bug in the 2decomp library.
If you are not sure this is really a bug in the library, you should go to the [discussions](https://github.com/2decomp-fft/2decomp-fft/discussions) section and open a new discussion.
Otherwise, follow the steps below.

Firstly, try to reproduce the error with a debug build of the library, a small problem size and a small number of MPI ranks.
It makes bug-hunting much easier.
Unfortunately, it is not always possible.
Please note that for a debug build, the log contains all the environment variables.
Use it to hunt the bug but think twice before sharing it as it can expose sensitive and personal information.
At least, please try to reproduce the bug on another machine with another compiler.

Secondly, if you have modified the source code of the 2decomp library, you must reproduce the bug without the modifications in 2decomp.
The development team will only provide support for sections of code available in the present repository.

Thirdly, you must provide a minimal working example.
The program using 2decomp and exposing the bug should be relatively small.
The development team will not provide support if the program exposing the bug is very long.
The programs available in the examples section are a good starting point for a minimal working example.

At this stage, you probably did your best to simplify the problem at hand.
Open a issue and select the bug report template.
Provide a meaningful title, do your best to complete all the sections of the template and provide the version of the compiler, the version of the MPI / FFT library, ...
If you think you have a fix for the bug, please expose it inside the issue.
It is recommended to wait for feedback before opening a pull request.

## Advanced contribution

One should read this section before opening a pull request.
To fix a bug, please open a issue and use the bug report template first.
To improve a part of the code or develop a new functionality, please open a issue and use the feature request template first.
If you are not sure about your contribution, open a [discussion](https://github.com/2decomp-fft/2decomp-fft/discussions) first.

The code in a pull-request should be formatted using the `fprettify` program.
See the code sample in the `scripts` folder.

Pull requests must be focused, small, coherent and have a detailed description.
The longer the pull request, the harder the review.
Please empathise with your fellow contributors who are going to spend time reviewing your code.

As long as the pull request is open for discussion and not ready for merging, convert it to draft.
Whenever it is ready for merging, convert it to a regular pull request.
Please note that reviewers might push modifications directly to your branch.
