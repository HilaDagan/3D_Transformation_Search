#ifndef P_FFTCOMMON_H
#define P_FFTCOMMON_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include <cassert>
#include <macros.h>

/**
 * A class with common static functions for performing FFT.
 */
class FFTCommon {
public:

    /**
     * Cast the given input to complex and store the result in 'out'
     */
    template<typename T1, typename T2>
    static void cast_values(const std::vector<T1> &in, std::vector<T2> & out) {
        for (unsigned long i = 0; i < out.size(); i++) {
            out[i] = in[i];
        }
    }

    static void complex_to_double(const std::vector<std::complex<double>>& in, std::vector<double>& out);

    static void dft(int n0, int n1, int n2, std::complex<double> *in,
                    std::complex<double> *out, int sign);

    static void fftw_forward(std::vector<std::complex<double>>& in,
                              std::vector<std::complex<double>>& out,
                              int nx, int ny, int nz);

    static void fftw_backward(std::vector<std::complex<double> >& in,
                              std::vector<std::complex<double> >& out,
                              int nx, int ny, int nz);
};

#endif //P_FFTCOMMON_H
