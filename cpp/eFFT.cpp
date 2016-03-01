#include <cmath>
#include <complex>
#include <stdlib.h>
#include "eFFT.h"

#define M_PI 3.141592653589793

namespace eFFT
{
int * eFFT::rrotate(int * input,unsigned int bits, int n){
	int lsb=0;
	for(int i=0;i<n;i++){	
		lsb=(input[i])&1;
		lsb<<=(bits-1);
		input[i]=lsb|((input[i])>>1);
	}
	return input;
}
int * eFFT::bitrev(int * order, unsigned int bits, int n){
	int r;
	int p;
	for(int k=0;k<n;k++){
		r=0;
		p=0;
		for(int i=bits-1;i>=0;i--){
			int x=pow(2,i);
			x&=k;
			x>>=i;
			p=x<<((bits-1)-(i));
			r|=p;
		}
		order[k]=r;
	}
	return order;
}
std::complex<double> * eFFT::butterfly(std::complex<double> * input, int siglen, int norm, int sign){
	std::complex<double> normalise = std::complex<double>(norm, 0);
	unsigned int poww = ceil(log2(siglen));
	int N = pow(2, poww);
	unsigned int bits = log2(N);
	int * order = new int[N];
	std::complex<double> * output;
	output = (std::complex<double>*)realloc(input, N * sizeof(std::complex<double>));
	for (int i = siglen; i < N; i++){
		//Zero Padding if required
		output[i] = std::complex<double>(0, 0);
	}
	int inc = 0;
	std::complex<double> tempup;
	std::complex<double> tempdown;
	std::complex<double> outputup;
	std::complex<double> outputdown;
	std::complex<double> twid;
	order = bitrev(order, bits, N);
	for (int i = 1; i < bits + 1; i++){
		inc = N / pow(2, i);
		int *w = new int[N / 2];
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
			twid = std::complex<double>(cos(2 * M_PI*w[wc] / N), sign*sin(2 * M_PI*w[wc] / N));
			outputup = tempup + (tempdown*twid);
			outputdown = tempup -(tempdown*twid);
			output[order[k]] = outputup;
			output[order[k + 1]] = outputdown;
		}
		order=rrotate(order,bits,N);
	}
	//Re-Order
	std::complex<double> * spectra = new std::complex<double>[N];
	for (int p = 0; p < N; p++){
		spectra[p] = output[order[p]]/normalise;
	}

	return spectra;
}
std::complex<double> * eFFT::fft(std::complex<double> * input,int length){
	std::complex<double> * spectra = butterfly(input, length, 1, -1);
	return spectra;
}
std::complex<double> * eFFT::ifft(std::complex<double> * input,int length){
	std::complex<double> * temporal = butterfly(input, length, length, 1);
	return temporal;
}
}
