# eFFT
The 'Easy Fast Fourier Transform' (eFFT) is a simple to use 1D radix-2 FFT library programmed in C.  It was developed late 2015/early 2016
##Authors
T.Lloyd and S.Ismail
##Background
The library was originally developed as a learning exercise to understand and implement the famous algorithm developed by Cooley and Tukey. Although many existing and faster implementations of the FFT exist the authors wanted to develop their own version to better understand the algorithm and to use it in projects.  Further information about how the FFT works can be found in the bibliography at the bottom.

##How To Use
The library contains two simple functions which are detailed below where the `input` is the data to under go transformation and 'length' is the number of data points.

###FFT
`double complex * fft(double complex * input,int length);`

###IFFT
'double complex * ifft(double complex * input,int length);'

##Bibliography
[1] J. W. Cooley, J. W. Tukey. An Algorithm for the Machine Calculation of Complex Fourier Series. (1965) Math. Comput.
[2] P. Gaydecki. Foundations of Digital Signal Processing: theory, algorithms and hardware design. (2004) IET Publishing.
[3] J. Wu. Fast Fourier Transform(FFT). Available: http://www.cmlab.csie.ntu.edu.tw/cml/dsp/training/coding/transform/fft.html. (N/A)