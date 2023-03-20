#!/usr/bin/env sh
#
# format.sh
#
# Runs the prettifier over 2decomp sources.

# fprettify options
# --enable-replacements : Replaces relational operators, e.g. .lt. -> <
# --c-relations : Use C-style relational operators, e.g. <=
# -w 3 : Whitespace preset - operators, print/read, *, /, +, -
fprettify --enable-replacements --c-relations -w 3 src/*f90 examples/*/*f90
