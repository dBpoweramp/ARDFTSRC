#pragma once
//---------------------------------------------------------------------------
// pocketfft backend for CARDFTSRC (transform_t), for any scalar type pocketfft_hdronly.h supports.
//
//  - Each object builds its plan (factorisation and twiddle factors) once, in init_*, and keeps it.  pocketfft's own
//    r2c() / c2r() build a new plan on every call unless a plan cache is configured, which used to happen for every
//    block.  The arithmetic is unchanged: the data goes through the same plan->exec as in r2c() / c2r().
//  - pocketfft_hdronly.h here is a patched copy.  For double it computes correctly rounded twiddle factors (in
//    double-double), which lowers the FFT rounding noise by 1-2 dB rms (up to 4 dB on the worst bin) and the
//    resampler's by about 1 dB.  Define ARDFTSRC_FFT_DOUBLE_TWIDDLES to get the stock double-precision twiddles
//    and output identical to earlier versions.
//  - An object can be moved but not copied.
//---------------------------------------------------------------------------
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#include "pocketfft_hdronly.h"
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>

template<typename _scalar_t = double>
class CFFT_pocketfft {
public:
    using scalar_t = _scalar_t;
    using complex_t = std::complex<scalar_t>;

    void init_forward(int size, const complex_t* in, complex_t* out)  { init_c2c(size, in, out, Mode::c2c_forward); }
    void init_backward(int size, const complex_t* in, complex_t* out) { init_c2c(size, in, out, Mode::c2c_backward); }

    void init_r2c(int size, const scalar_t* in, complex_t* out) {
        n = static_cast<size_t>(size);
        rplan = std::make_unique<pocketfft::detail::pocketfft_r<scalar_t>>(n);
        cplan.reset();
        rwork.assign(n, scalar_t(0));
        mode = Mode::r2c; in_r = in; out_c = out;
    }
    void init_c2r(int size, const complex_t* in, scalar_t* out) {
        n = static_cast<size_t>(size);
        rplan = std::make_unique<pocketfft::detail::pocketfft_r<scalar_t>>(n);
        cplan.reset();
        rwork.assign(n, scalar_t(0));
        mode = Mode::c2r; in_c = in; out_r = out;
    }

    void execute() {
        using cmplx = pocketfft::detail::cmplx<scalar_t>;
        switch (mode) {
        case Mode::r2c: {
            std::copy(in_r, in_r + n, rwork.data());
            rplan->exec(rwork.data(), scalar_t(1), true);
            // halfcomplex [r0, r1, i1, r2, i2, ...] to complex, as pocketfft's r2c() does
            out_c[0] = complex_t(rwork[0], scalar_t(0));
            size_t i = 1, k = 1;
            for (; i + 1 < n; i += 2, ++k) out_c[k] = complex_t(rwork[i], rwork[i + 1]);
            if (i < n) out_c[k] = complex_t(rwork[i], scalar_t(0));
            break;
        }
        case Mode::c2r: {
            // complex to halfcomplex, as pocketfft's c2r() does for a backward transform
            rwork[0] = in_c[0].real();
            size_t i = 1, k = 1;
            for (; i + 1 < n; i += 2, ++k) { rwork[i] = in_c[k].real(); rwork[i + 1] = in_c[k].imag(); }
            if (i < n) rwork[i] = in_c[k].real();
            rplan->exec(rwork.data(), scalar_t(1), false);
            std::copy(rwork.begin(), rwork.end(), out_r);
            break;
        }
        case Mode::c2c_forward:
        case Mode::c2c_backward:
            std::copy(in_c, in_c + n, cwork.data());
            cplan->exec(reinterpret_cast<cmplx*>(cwork.data()), scalar_t(1), mode == Mode::c2c_forward);
            std::copy(cwork.begin(), cwork.end(), out_c);
            break;
        case Mode::none:
            break;
        }
    }

    // FFT library quirks
    // Does FFT alter input buffer? If not, some buffer nulling can be skipped
    static constexpr bool alters_input_buffer = false;
    // Does FFT library output entire mirrored range, or only first half of it?
    static constexpr bool full_output = false;
    static int r2c_size(int inSize) { return full_output ? inSize : inSize / 2 + 1; }

private:
    enum class Mode { none, r2c, c2r, c2c_forward, c2c_backward };

    void init_c2c(int size, const complex_t* in, complex_t* out, Mode m) {
        n = static_cast<size_t>(size);
        cplan = std::make_unique<pocketfft::detail::pocketfft_c<scalar_t>>(n);
        rplan.reset();
        cwork.assign(n, complex_t(scalar_t(0), scalar_t(0)));
        mode = m; in_c = in; out_c = out;
    }

    Mode mode = Mode::none;
    size_t n = 0;
    const scalar_t* in_r = nullptr;
    const complex_t* in_c = nullptr;
    scalar_t* out_r = nullptr;
    complex_t* out_c = nullptr;
    std::unique_ptr<pocketfft::detail::pocketfft_r<scalar_t>> rplan;
    std::unique_ptr<pocketfft::detail::pocketfft_c<scalar_t>> cplan;
    std::vector<scalar_t> rwork;
    std::vector<complex_t> cwork;
};
