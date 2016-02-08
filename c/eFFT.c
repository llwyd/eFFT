//eFFT.c
//T.Lloyd 
//2015-2016

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#define M_PI 3.141592653589793

char * binaryadd(char *bits, unsigned int position, unsigned int totalbits){
	if (((bits[position] + 1)) > 1){
		bits[position] = 0;
		if (position == 0){
			bits = binaryadd(bits, totalbits - 1, totalbits);
		}
		else{
			bits = binaryadd(bits, (position - 1) % totalbits, totalbits);
		}
	}
	else{
		bits[position] += 1;
		return bits;
	}
}
int * bitrev(int * order, unsigned int bits, int N, int div){
	char * indiv = calloc(bits,sizeof(char));
	int temp = 0;
	order[0] = 0;
	for (int i = 0; i < N - 1; i++){
		for (int k = 0; k < bits; k++){
			temp = ((order[i] >> k) & 1);
			indiv[k] = (char)temp;
		}
		indiv = binaryadd(indiv, bits - div, bits);
		order[i + 1] = 0;
		for (int k = 0; k < bits; k++){
			order[i + 1] += (indiv[k] * pow(2, k));

		}
	}
	return order;
}

double complex * butterfly(complex double * input, int siglen, int norm, int sign){
	complex double normalise = (complex double)norm;
	unsigned int poww = ceil(log2(siglen));
	int N = pow(2, poww);
	unsigned int bits = log2(N);
	int * order = calloc(N,sizeof(int));
	complex double * output;
	output = (complex double*)realloc(input, N * sizeof(complex double));
	for (int i = siglen; i < N; i++){
		//Zero Padding if required
		output[i] = 0;
	}
	int inc = 0;
	complex double tempup;
	complex double tempdown;
	complex double outputup;
	complex double outputdown;
	complex double twid;
	for (int i = 1; i < bits + 1; i++){
		inc = N / pow(2, i);
		order = bitrev(order, bits, N, i);
		int * w = calloc(N/2,sizeof(int));
		int incre = 0;
		for (int j = 0; j < (N / 2); j++){
			if ((j%inc == 0) && (j != 0)){
				incre = (incre + inc) % (N / 2);
			}
			w[j] = incre;
		}
		incre = 0;
		int wc = 0;
		for (int k = 0,wc=0; k < N; k += 2,wc++){
			tempup = output[order[k]];
			tempdown = output[order[k + 1]];
			twid = (complex double)(cos(2 * M_PI*w[wc] / N))+(I*sign*sin(2 * M_PI*w[wc] / N));
			outputup = tempup + (tempdown*twid);
			outputdown = tempup -(tempdown*twid);
			output[order[k]] = outputup;
			output[order[k + 1]] = outputdown;
		}
	}
	//Re-Order
	complex double * spectra = calloc(N,sizeof(complex double));
	order = bitrev(order, bits, N, 1);
	for (int p = 0; p < N; p++){
		spectra[p] = output[order[p]]/normalise;
	}

	return spectra;
}


double complex * fft(double complex * input,int length){
	double complex * spectra = butterfly(input, length, 1, -1);
	return spectra;
}
double complex * ifft(double complex * input,int length){
	double complex * temporal = butterfly(input, length, length, 1);
	return temporal;
}