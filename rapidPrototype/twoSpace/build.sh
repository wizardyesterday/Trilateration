#! /bin/sh

gcc -g -O0 $1 -l lapack -l blas -l f2c -lm
