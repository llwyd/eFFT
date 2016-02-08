//eFFT.h
//cpp header
//T.Lloyd & S.Ismail
//2015-2016


#include <cmath>
#include <complex>
class eFFT
{
private:
		char * binaryadd(char *bits, unsigned int position, unsigned int totalbits);
		int * bitrev(int * order, unsigned int bits, int N, int div);
		std::complex<double> * butterfly(std::complex<double> * input, int siglen, int norm, int sign);
public:
	std::complex<double> * fft(std::complex<double> * input,int length);
	std::complex<double> * ifft(std::complex<double> * input,int length);
}