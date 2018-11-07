
#include "FFTCommon.h"


/**
 * Cast the given input to double and store the result in 'out'.
 */
void FFTCommon::complex_to_double(const std::vector<std::complex<double>>& in,
                                             std::vector<double>& out) {
    for(unsigned long i=0; i<out.size(); i++) {
        out[i] = std::norm(in[i]);
    }
}

/**
 * Perform DFT.
 * @param n0, n1, n2 - gives array (density map) dimensions.
 * @param in, out - input/output arrays
 * @param sign - the direction of the transform FFTW_FORWARD (-1) or FFTW_BACKWARD (+1)
 */
void FFTCommon::dft(int n0, int n1, int n2, std::complex<double> *in,
         std::complex<double> *out, int sign){
    fftw_plan plan = fftw_plan_dft_3d(n0, n1, n2, (fftw_complex *) in,
                                      (fftw_complex *) out, sign,
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    assert(plan != nullptr);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * Transform to Fourier space.
 */
void FFTCommon::fftw_forward(std::vector<std::complex<double>>& in,
                                        std::vector<std::complex<double>>& out,
                                        int nx, int ny, int nz)
{
    dft(nz, ny, nx, &in[0], &out[0], FFTW_FORWARD);
}

/**
 * Transform from Fourier space.
 */
void FFTCommon::fftw_backward(std::vector<std::complex<double> >& in,
                               std::vector<std::complex<double> >& out,
                               int nx, int ny, int nz)
{
    dft(nz, ny, nx, &in[0], &out[0], FFTW_BACKWARD);
}
