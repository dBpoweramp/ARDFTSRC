#pragma once
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#include "pocketfft_hdronly.h"
#include <functional>
#include <complex>

template<typename _scalar_t = double>
class CFFT_pocketfft {
public:
    using scalar_t = _scalar_t;
    using complex_t = std::complex<scalar_t>;
    
    void init_forward(int size, const complex_t * in, complex_t * out) {
        using namespace pocketfft;
        m_execute = [this, in, out, shape = shape_t{ (size_t)size }, stride_in = stride_t{ sizeof(*in) }, stride_out = stride_t{sizeof(*out)}] {
            c2c(shape, stride_in, stride_out, {0}, FORWARD, in, out, (scalar_t)1.0);
        };
    }
    void init_backward(int size, const complex_t * in, complex_t * out) {
        using namespace pocketfft;
        m_execute = [this, in, out, shape = shape_t{ (size_t)size }, stride_in = stride_t{ sizeof(*in) }, stride_out = stride_t{sizeof(*out)}] {
            c2c(shape, stride_in, stride_out, {0}, BACKWARD, in, out, (scalar_t)1.0);
        };
    }
    
    void init_r2c(int size, const scalar_t* in, complex_t* out) {
        using namespace pocketfft;
        m_execute = [this, in, out, shape = shape_t{ (size_t)size }, stride_in = stride_t{ sizeof(*in) }, stride_out = stride_t{sizeof(*out)}] {
            r2c(shape, stride_in, stride_out, 0, FORWARD, in, out, (scalar_t)1.0);
        };
    }
    void init_c2r(int size, const complex_t* in, scalar_t* out) {
        using namespace pocketfft;
        m_execute = [this, in, out, shape = shape_t{ (size_t)size }, stride_in = stride_t{ sizeof(*in) }, stride_out = stride_t{ sizeof(*out) }] {
            c2r(shape, stride_in, stride_out, 0, BACKWARD, in, out, (scalar_t)1.0);
        };
    }
    void execute() { 
        m_execute();
    }
    // FFT library quirks
    // Does FFT alter input buffer? If not, some buffer nulling can be skipped
    static constexpr bool alters_input_buffer = false;
    // Does FFT library output entire mirrored range, or only first half of it?
    static constexpr bool full_output = false;
    static int r2c_size(int inSize) { return full_output ? inSize : inSize / 2 + 1; }
private:
    std::function<void()> m_execute;
};
