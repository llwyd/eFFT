#ifndef EFFT_H
#define EFFT_H

//eFFT.h
//cpp header
//T.Lloyd & S.Ismail
//2015-2016

#include <cmath>
#include <complex>

namespace eFFT
{
class eFFT
{
private:
	static int * rrotate(int * input,unsigned int bits, int n);
	static int * bitrev(int * order, unsigned int bits, int n);
	static std::complex<double> * butterfly(std::complex<double> * input, int siglen, int norm, int sign);
public:
	static std::complex<double> * fft(std::complex<double> * input,int length);
	static std::complex<double> * ifft(std::complex<double> * input,int length);
};
}

#endif // EFFT_H
