#!/bin/sh

gcc -g -O0 Trilaterate.c matutil.c  -o trilaterate -l lapack -l blas -l f2c -lm
