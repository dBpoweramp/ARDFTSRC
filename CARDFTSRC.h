#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include <numbers>
#include <limits>
#include <algorithm>
#include <complex>
#include <cstdint>
#include <type_traits>

#include "FFT-pocketfft.h"

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

template<typename transform_t = CFFT_pocketfft<> >
class CARDFTSRC {
public:
    static constexpr double default_bandwidth = 0.97;
    static constexpr unsigned default_quality = 8192;
    using complex_t = typename transform_t::complex_t;
    using scalar_t = typename transform_t::scalar_t;

    // Fused multiply-add for the built-in types; plain x*y+z for class types (e.g. boost cpp_double_double).
    // (if constexpr instead of in-class explicit specialisations, which only MSVC accepts)
    static scalar_t my_fma(scalar_t x, scalar_t y, scalar_t z) {
        if constexpr (std::is_floating_point_v<scalar_t>) return std::fma(x, y, z);
        else return x * y + z;
    }
    
    //----------------------------------------------------------------------------------------------------------
  /*  static double planck_taper_down(double u, double alpha = 2.0) {
        if (u <= 0.0) return 1.0;
        if (u >= 1.0) return 0.0;

        // Descending logistic/Planck-like transition.
        // alpha = 2.0 is equivalent to the current tanh shape.
        const double x = alpha * (1.0 / (1.0 - u) - 1.0 / u);

        // Stable computation of 1 / (1 + exp(x)).
        if (x >= 0.0) {
            const double e = std::exp(-x);
            return e / (1.0 + e);
        }
        else {
            const double e = std::exp(x);
            return 1.0 / (1.0 + e);
        }
    }*/

    //----------------------------------------------------------------------------------------------------------
    // smooths the FFT bin selection to find a macthing even factor
    //----------------------------------------------------------------------------------------------------------
    static bool is_fft_smooth_factor(std::uint64_t n) {
        for (std::uint64_t p : {2ull, 3ull, 5ull, 7ull, 11ull}) {
            while (n % p == 0)
                n /= p;
        }

        return n == 1;
    }

    //----------------------------------------------------------------------------------------------------------
    static std::uint64_t next_fft_smooth_even_factor(std::uint64_t factor) {
        if (factor < 1)
            factor = 1;

        if (factor & 1u)
            ++factor;

        while (!is_fft_smooth_factor(factor))
            factor += 2;

        return factor;
    }

    //----------------------------------------------------------------------------------------------------------
    // rational_cosine_taper_down 
    //----------------------------------------------------------------------------------------------------------
    static double choose_taper_p(int transition_bins) {
        if (transition_bins < 6000)
            return 3.25;

        if (transition_bins < 10300)
            return 3.5;

        if (transition_bins < 11200)
            return 3.625;

        return 3.75;
    }

    //----------------------------------------------------------------------------------------------------------
    template<class T>
    static T rational_cosine_taper_down(T x, T p) {
        x = std::clamp(x, T(0), T(1));

        if (x <= T(0))
            return T(1);

        if (x >= T(1))
            return T(0);

        const T y = T(1) - x;

        // w = x^p / (x^p + y^p), computed stably.
        const T r = p * (std::log(x) - std::log(y));

        T w;

        if (r >= T(0)) {
            const T e = std::exp(-r);
            w = T(1) / (T(1) + e);
        }
        else {
            const T e = std::exp(r);
            w = e / (T(1) + e);
        }

        const T q = T(1) - w;

        const T half_pi = T(0.5) * std::numbers::pi_v<T>;

        T v;

        if (w <= q) {
            const T c = std::cos(half_pi * w);
            v = c * c;
        }
        else {
            const T s = std::sin(half_pi * q);
            v = s * s;
        }

        if (!std::isfinite(v))
            return T(0);

        return std::clamp(v, T(0), T(1));
    }

    //----------------------------------------------------------------------------------------------------------
    void build_rational_cosine_taper(void)
    {
        taper.resize(transition_bins);
        double taper_p = choose_taper_p(transition_bins);
        for (int i = 0; i < transition_bins; ++i) {
            const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(transition_bins);
            const double v = rational_cosine_taper_down<double>(x, static_cast<double>(taper_p));

            taper[i] = static_cast<scalar_t>(std::clamp(v, 0.0, 1.0));
        }
    }

    //----------------------------------------------------------------------------------------------------------
    // beta cdf taper
    //----------------------------------------------------------------------------------------------------------
    template<class T>
    static T clamp01_finite(T v) {
        if (!std::isfinite(v))
            return T(0);

        if (v < T(0))
            return T(0);

        if (v > T(1))
            return T(1);

        return v;
    }

    //----------------------------------------------------------------------------------------------------------
    // Computes I_x(a, a) for integer a.
    // Intended for x in [0, 0.5]. Uses double internally.
    //----------------------------------------------------------------------------------------------------------
    static double beta_cdf_equal_integer_low(double x, int order) {
        if (x <= 0.0)
            return 0.0;

        if (x >= 0.5)
            return 0.5;

        if (order <= 0)
            return x;

        const int n = 2 * order - 1;
        const double y = 1.0 - x;

        // term for j = 0:
        // C(n,0) * x^0 * y^n
        double term = std::pow(y, n);
        double sum = 0.0;

        for (int j = 0; j <= n; ++j) {
            if (j >= order)
                sum += term;

            if (j == n)
                break;

            // Recurrence:
            // term[j + 1] = term[j] * ((n - j) / (j + 1)) * (x / y)
            term *= (static_cast<double>(n - j) /
                static_cast<double>(j + 1)) *
                (x / y);
        }

        return std::clamp(sum, 0.0, 0.5);
    }

    //----------------------------------------------------------------------------------------------------------
    void build_beta_cdf_taper(int order) 
    {
        taper.resize(static_cast<std::size_t>(transition_bins));

        if (transition_bins <= 0)
            return;

        //constexpr double snap = 1.0e-18L; // about -360 dB amplitude
        constexpr double snap = 1.0e-16; // about -320 dB amplitude

        const int M = transition_bins;
        const int half = M / 2;

        for (int i = 0; i < half; ++i) {
            const double x =
                (static_cast<double>(i) + 0.5) /
                static_cast<double>(M);

            // lo = I_x(order, order), near zero at the beginning.
            double lo = beta_cdf_equal_integer_low(x, order);
            double hi = 1.0 - lo;

            if (lo < snap)
                lo = 0.0;

            if (1.0 - hi < snap)
                hi = 1.0;

            taper[static_cast<std::size_t>(i)] =
                static_cast<scalar_t>(hi);

            taper[static_cast<std::size_t>(M - 1 - i)] =
                static_cast<scalar_t>(lo);
        }

        if (M & 1) {
            taper[static_cast<std::size_t>(half)] = scalar_t{ 0.5 };
        }

        // Monotonicity repair for rounding noise.
        for (int i = 1; i < M; ++i) {
            auto& cur = taper[static_cast<std::size_t>(i)];
            const auto& prev = taper[static_cast<std::size_t>(i - 1)];

            if (cur > prev)
                cur = prev;
        }
    }

    //----------------------------------------------------------------------------------------------------------
    // Kaiser-kernel taper
    //
    //   taper(x) = ( integral_{x}^{1} g(v) dv ) / ( integral_{0}^{1} g(v) dv ),   x = (i + 0.5) / M
    //   g(v)     = I0( beta * sqrt(1 - (2v - 1)^2) )               (Kaiser window across the transition band)
    //
    // The impulse response of the block filter is  h(t) = G(t) * sin(2*pi*fc*t) / (pi*t)  where G is the
    // Fourier transform of the kernel g (the derivative of the taper).  For the beta-CDF / rational-cosine
    // tapers G(t) only decays algebraically (~t^-(order+1)); for the Kaiser kernel it is ~exp(-beta) beyond
    // its main lobe (which ends at t = beta / (pi * W), W = transition width in Hz).  For the same transition
    // width this lowers the residual impulse response at +/-0.3 s (the window impulse.m truncates to) by
    // >100 dB, or equivalently allows a ~40% narrower transition band for the same residual.
    //
    // Series used:  I0(b*sqrt(1-v^2)) = sum_k (b/2)^(2k)/(k!)^2 * (1-v^2)^k
    //               J_k(u) = integral_0^u (1-v^2)^k dv ,  J_0 = u ,  J_k = (2k*J_{k-1} + u*(1-u^2)^k) / (2k+1)
    // All terms are positive, so the recursion is stable (integrating outward from the centre of the kernel).
    //----------------------------------------------------------------------------------------------------------
    void build_kaiser_taper(double beta)
    {
        using T = scalar_t;         

        // Normal mode: M = transition_bins points at x = (i + 0.5) / M on bins [tr_nb - M, tr_nb).
        // Straddle mode: M + 1 points at x = i / M on bins [tr_nb - M/2, tr_nb + M/2], so that
        //   taper[M/2 - a] + taper[M/2 + a] == 1  (bin tr_nb - a and its alias bin tr_nb + a are complementary).
        const int M = transition_bins;
        const int npts = straddle_nyquist ? M + 1 : M;
        taper.resize(static_cast<std::size_t>(npts));
        if (M <= 0)
            return;

        if (beta <= 0.0)
            beta = 1.0e-6;                      // degenerates to a linear taper

        // c_k = (beta/2)^(2k) / (k!)^2
        std::vector<T> coef;
        {
            const T q = T(beta) * T(beta) / T(4);
            const T tiny = std::numeric_limits<T>::epsilon() * T(1.0e-4);
            T c = T(1), cmax = T(1);
            for (int k = 1; k < 100000; ++k) {
                coef.push_back(c);
                if (c > cmax) cmax = c;
                c *= q / (T(k) * T(k));
                if (k > 8 && c < cmax * tiny)
                    break;
            }
        }

        // integral_0^u I0(beta*sqrt(1-v^2)) dv  for u in [0,1]
        auto half_integral = [&coef](T u) -> T {
            const T w = T(1) - u * u;
            T Jk = u, wk = T(1), sum = coef[0] * Jk;
            for (std::size_t k = 1; k < coef.size(); ++k) {
                wk *= w;
                Jk = (T(2 * k) * Jk + u * wk) / T(2 * k + 1);
                sum += coef[k] * Jk;
            }
            return sum;
        };

        const T half_area = half_integral(T(1));

        const int half = npts / 2;

        for (int i = 0; i < half; ++i) {
            const T x = straddle_nyquist ? T(i) / T(M)                    // 0 .. 0.5, passband side of the transition
                                         : (T(i) + T(0.5)) / T(M);
            const T u = T(1) - T(2) * x;                               // 1 .. 0,   distance from centre, kernel units
            T lo = T(0.5) * (T(1) - half_integral(u) / half_area);     // ~0 at the Nyquist end
            if (lo < T(0)) lo = T(0);
            if (lo > T(0.5)) lo = T(0.5);

            taper[static_cast<std::size_t>(i)] = static_cast<scalar_t>(T(1) - lo);
            taper[static_cast<std::size_t>(npts - 1 - i)] = static_cast<scalar_t>(lo);
        }

        if (npts & 1)
            taper[static_cast<std::size_t>(half)] = scalar_t{ 0.5 };

        // Monotonicity repair for rounding noise.
        for (int i = 1; i < npts; ++i) {
            auto& cur = taper[static_cast<std::size_t>(i)];
            const auto& prev = taper[static_cast<std::size_t>(i - 1)];
            if (cur > prev)
                cur = prev;
        }
    }
    

    CARDFTSRC(int input_samplerate_, int output_samplerate_, int channels_, unsigned quality = default_quality, double bandwidth = default_bandwidth, double phaset = 0.0, int Betacdftaper = 0, double kaiser_beta = 0.0, bool straddle_nyquist_ = false) :
        input_samplerate(input_samplerate_), output_samplerate(output_samplerate_), num_channels(channels_) {
        int common_divisor = std::gcd(input_samplerate, output_samplerate);
        in_nb_samples = input_samplerate / common_divisor;
        out_nb_samples = output_samplerate / common_divisor;

        int factor = (int)next_fft_smooth_even_factor(((int64_t)(quality) + (int64_t)(std::min(in_nb_samples, out_nb_samples)) - 1) / (int64_t)(std::min(in_nb_samples, out_nb_samples)));

        in_nb_samples *= factor;
        out_nb_samples *= factor;

        in_rdft_size = in_nb_samples * 2;
        out_rdft_size = out_nb_samples * 2;
        in_offset = (in_rdft_size - in_nb_samples) / 2;
        out_offset = (out_rdft_size - out_nb_samples) / 2; // This is not used in the original logic but kept for consistency
        tr_nb_samples = std::min(in_nb_samples, out_nb_samples);
        transition_bins = static_cast<int>(std::ceil(tr_nb_samples * (1.0 - bandwidth)));

        prev_chunk.resize(out_nb_samples * num_channels, scalar_t{});

        transition_bins = std::clamp(transition_bins, 1, (int)tr_nb_samples);

        // Straddle mode: the transition band is centred on the output Nyquist frequency instead of ending there.
        // The half above Nyquist is folded (aliased) back onto the half below, exactly as sampling at the output
        // rate would do, so the effective filter is a smooth linear-phase transition straddling Nyquist:
        // the -3 dB point sits at Nyquist and the noise-floor PSD stays flat to within ~W/4 of it, at the cost of
        // aliasing the source band [Nyquist, Nyquist + W/2] into [Nyquist - W/2, Nyquist]. 
        // Only meaningful when down-sampling and with the Kaiser taper.
        straddle_nyquist = straddle_nyquist_ && kaiser_beta > 0.0 && in_nb_samples > out_nb_samples;
        if (straddle_nyquist) {
            transition_bins = std::clamp(transition_bins + (transition_bins & 1), 2, (int)tr_nb_samples);   // even
            if (tr_nb_samples + transition_bins / 2 >= in_nb_samples)
                straddle_nyquist = false;        // no room above Nyquist in the input spectrum
        }
        passband_bins = straddle_nyquist ? tr_nb_samples - transition_bins / 2 : tr_nb_samples - transition_bins;

        if (kaiser_beta > 0.0)
            build_kaiser_taper(kaiser_beta);     // preferred: fastest impulse-response decay for a given transition width
        else if (Betacdftaper == 0)
            build_rational_cosine_taper();
        else
            build_beta_cdf_taper(Betacdftaper);

        if (phaset != 0) {
            phase.resize(tr_nb_samples);
            for (int n = 0; n < tr_nb_samples; ++n) {
                const double x = (double)n / tr_nb_samples;
                constexpr double shift = -0.4;
                using namespace std::numbers;
                const double z = (sqrt2 / 2) * (exp((x + shift) * ln10 * pi) - exp(shift * ln10 * pi));
                const double w = -phaset * z;

                phase[n] = complex_t{ cos(w), sin(w) };
            }
        }

        // FFT setup
        rdftc.resize(std::max(in_rdft_size, out_rdft_size) / 2 + 1, (scalar_t)0.0f);
        fftw_in.resize(in_rdft_size);
        if constexpr (!transform_t::alters_input_buffer) std::fill(fftw_in.begin(), fftw_in.end(), 0.0f);
        rfft.init_r2c(in_rdft_size, fftw_in.data(), rdftc.data());
        ifftw_out.resize(out_rdft_size);
        irfft.init_c2r(out_rdft_size, rdftc.data(), ifftw_out.data());
        scale = (double)out_nb_samples / (double)in_nb_samples;
    }
    ~CARDFTSRC() {
    }

    //----------------------------------------------------------------------------------------------------------
    template<typename sample_t> // float or double
    void process(const sample_t* inAudio, sample_t* outAudio) {
        process_ex([inAudio](size_t idx) { return inAudio[idx]; }, [outAudio](size_t idx, double arg) { outAudio[idx] = (sample_t)arg; });
    }
    //----------------------------------------------------------------------------------------------------------
    template<typename sample_t> // float or double
    void process(const sample_t* inAudio, std::nullptr_t) {
        process_ex([inAudio](size_t idx) { return inAudio[idx]; }, [](size_t, double) {});
    }
    //----------------------------------------------------------------------------------------------------------
    template<typename get_t, typename put_t>
    void process_ex(get_t && get, put_t && put) {
        for (int ch = 0; ch < num_channels; ++ch) {
            int inPtr = ch;
            // clear zero padding
            std::fill(fftw_in.begin(),
                fftw_in.begin() + in_offset,
                scalar_t{});

            std::fill(fftw_in.begin() + in_offset + in_nb_samples,
                fftw_in.end(),
                scalar_t{});

            for (int k = 0; k < in_nb_samples; ++k) {
                fftw_in[k + in_offset] = (scalar_t)get(inPtr); inPtr += num_channels;
            }
            if constexpr (transform_t::alters_input_buffer) std::fill(fftw_in.begin() + in_offset + in_nb_samples, fftw_in.end(), 0.0f);
            // RFFT
            rfft.execute();

            if (straddle_nyquist) {
                // Phase first: bin tr+a carries its own phase into the fold, where conjugation mirrors it,
                // which is what sampling at the output rate does to it.
                for (size_t n = 0; n < phase.size(); ++n) rdftc[n] *= phase[n];

                // Output bin tr-a = direct term of input bin tr-a  +  alias of input bin tr+a (conjugated, its own weight).
                const std::size_t tr = static_cast<std::size_t>(tr_nb_samples);
                const int half = transition_bins / 2;
                for (int a = 1; a <= half; ++a) {
                    complex_t& lo = rdftc[tr - static_cast<std::size_t>(a)];
                    const complex_t& hi = rdftc[tr + static_cast<std::size_t>(a)];
                    const complex_t hi_conj{ hi.real(), -hi.imag() };
                    lo = lo * taper[static_cast<std::size_t>(half - a)] + hi_conj * taper[static_cast<std::size_t>(half + a)];
                }
                // The c2r counts the Nyquist bin once: both +/-Nyquist input components land on it.
                rdftc[tr] = complex_t(scalar_t(2) * std::real(rdftc[tr]) * taper[static_cast<std::size_t>(half)], scalar_t(0));
                for (size_t n = tr + 1; n < rdftc.size(); ++n) rdftc[n] = 0;
            }
            else {
                for (size_t n = (size_t)tr_nb_samples; n < rdftc.size(); ++n) rdftc[n] = 0;

                for (size_t n = 0; n < phase.size(); ++n) rdftc[n] *= phase[n];

                for (int i = 0; i < transition_bins; ++i) {
                    rdftc[static_cast<std::size_t>(passband_bins + i)] *= taper[i];
                }
            }

            // IRFFT
            irfft.execute();

            // Normalize and store
            // A wide class type (double-double) computes the gain in its own precision rather than inheriting a
            // double-rounded one; the built-in types keep the original double computation.
            const scalar_t factor = std::is_floating_point_v<scalar_t>
                ? static_cast<scalar_t>(scale / static_cast<double>(out_rdft_size))
                : scalar_t(out_nb_samples) / (scalar_t(in_nb_samples) * scalar_t(out_rdft_size));

            {
                auto outPtr = ch;
                auto inPrev = prev_chunk.data() + ch;

                for (int k = 0; k < out_nb_samples; ++k) {
                    const scalar_t y = my_fma( ifftw_out[static_cast<std::size_t>(k)], factor, *inPrev);

                    put(outPtr, static_cast<double>(y));

                    outPtr += num_channels;
                    inPrev += num_channels;
                }
            }

            {
                auto outPtr = prev_chunk.data() + ch;

                for (int k = 0; k < out_nb_samples; ++k) {
                    *outPtr = ifftw_out[static_cast<std::size_t>(k + out_nb_samples)] * factor;

                    outPtr += num_channels;
                }
            }
        }
    }

    int InChunkSize() const { return(in_nb_samples); }
    int OutChunkSize() const { return out_nb_samples; }
    int FirstBlocksSampleDelay() const { return(out_offset); }
    int FinalBlocksZeroPadSample() const { return(in_offset); }

    const int input_samplerate, output_samplerate, num_channels;

    int in_nb_samples = 0;
    int out_nb_samples = 0;

    int in_rdft_size;
    int out_rdft_size;
    int in_offset;
    int out_offset;
    int tr_nb_samples;

    std::vector<scalar_t> prev_chunk;
    
    int transition_bins = 0;
    int passband_bins = 0;
    bool straddle_nyquist = false;
    std::vector<scalar_t> taper; 
    

    std::vector<scalar_t> fftw_in;
    std::vector<complex_t> rdftc;
    std::vector<scalar_t> ifftw_out;
    transform_t rfft, irfft;
    double scale = 0;

    std::vector<complex_t> phase;
};
