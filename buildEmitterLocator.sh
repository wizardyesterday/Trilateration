#!/bin/sh

gcc -c -g -O0 -I include src/matutil.c

g++ -g -O0 src/main.cc src/EmitterLocator.cc matutil.o  -o emitterLocator -I include -l lapack -l blas -l f2c -lm
