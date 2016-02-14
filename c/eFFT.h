//eFFT.h
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

double complex * fft(double complex * input,int length);
double complex * ifft(double complex * input,int length);