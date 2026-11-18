#if __has_include("stdafx.h")
#include "stdafx.h"
#endif
#include "CFFT-fftw3.h"
#include <mutex>
static std::mutex fft_initGuard;

void CFFT_fftw3::init_r2c(int size, scalar_t* in, complex_t* out) {
    std::scoped_lock lock(fft_initGuard);
    _destroy();
    plan = fftw_plan_dft_r2c_1d(size, in, reinterpret_cast<fftw_complex*>(out), FFTW_ESTIMATE);
    if (!plan) throw std::bad_alloc();
}
void CFFT_fftw3::init_c2r(int size, complex_t* in, scalar_t* out) {
    std::scoped_lock lock(fft_initGuard);
    _destroy();
    plan = fftw_plan_dft_c2r_1d(size, reinterpret_cast<fftw_complex*>(in), out, FFTW_ESTIMATE);
    if (!plan) throw std::bad_alloc();
}
void CFFT_fftw3::destroy() {
    std::scoped_lock lock(fft_initGuard);
    _destroy();
}
