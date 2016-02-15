//eFFT.h
//cpp header
//T.Lloyd & S.Ismail
//2015-2016


#include <cmath>
#include <complex>
class eFFT
{
private:
	int * rrotate(int * input,unsigned int bits, int n);
	int * bitrev(int * order, unsigned int bits, int n);
	std::complex<double> * butterfly(std::complex<double> * input, int siglen, int norm, int sign);
public:
	std::complex<double> * fft(std::complex<double> * input,int length);
	std::complex<double> * ifft(std::complex<double> * input,int length);
	eFFT();
	~eFFT();
};