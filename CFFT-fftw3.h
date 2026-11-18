#pragma once
#if __has_include(<fftw/api/fftw3.h>)
#include <fftw/api/fftw3.h>
#else
#include <fftw3.h>
#endif
#include <complex>

class CFFT_fftw3 {
public:
    typedef double scalar_t;

    // These types have identical binary layout, we use them interchangably, to avoid issues with fftw_complex in a vector.
    typedef std::complex<scalar_t> complex_t;
    static_assert(sizeof(fftw_complex) == sizeof(complex_t));

    void init_r2c(int size, scalar_t* in, complex_t* out);
    void init_c2r(int size, complex_t* in, scalar_t* out);
    void destroy();
    ~CFFT_fftw3() { destroy(); }
    void execute() { fftw_execute(plan); }
    
    // FFT library quirks
    // Does FFT alter input buffer? If not, some buffer nulling can be skipped
    static constexpr bool alters_input_buffer = false;
    // Does FFT library output entire mirrored range, or only first half of it?
    static int r2c_size(int inSize) { return inSize / 2 + 1; }

private:
    void _destroy() {
        fftw_destroy_plan(plan);
    }
    fftw_plan plan = nullptr;
};

#ifdef _MSC_VER
#pragma comment(lib, "fftw3.lib")
#endif