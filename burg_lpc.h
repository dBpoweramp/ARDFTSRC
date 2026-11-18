#pragma once
// Burg-method linear prediction, used to extrapolate audio past the start or end of a file.
//
// Compared with the autocorrelation method in vorbis_lpc, Burg's method estimates the predictor from the data alone
// (no implied zeros outside the analysis window), which resolves tonal components far more sharply, and it always
// gives reflection coefficients |k| < 1, i.e. a stable predictor.  The extrapolation runs through the all-pole
// lattice filter of those coefficients rather than the equivalent direct-form recursion: with poles close to the
// unit circle (tones, band-limited material) the direct form can diverge in double precision while the lattice
// stays stable.
#include <vector>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <algorithm>

namespace burg_lpc {

    // Reflection coefficients of an order 'order' predictor for x[0..n) (order is reduced to n / 2 for short input;
    // it stops early on silence).
    inline std::vector<double> reflection_coefficients(const double* x, size_t n, size_t order)
    {
        std::vector<double> k;
        if (order > n / 2) order = n / 2;
        if (order == 0) return k;
        k.reserve(order);
        std::vector<double> f(x, x + n), b(x, x + n);     // forward / backward prediction errors
        for (size_t m = 1; m <= order; ++m) {
            double num = 0, den = 0;
            for (size_t i = m; i < n; ++i) {
                num += f[i] * b[i - 1];
                den += f[i] * f[i] + b[i - 1] * b[i - 1];
            }
            if (!(den > 0)) break;
            const double km = -2.0 * num / den;
            k.push_back(km);
            for (size_t i = n - 1; i >= m; --i) {           // descending, so b[i - 1] is still the previous order's
                const double fi = f[i], bi = b[i - 1];
                f[i] = fi + km * bi;
                b[i] = bi + km * fi;
            }
        }
        return k;
    }

    // Continues x[0..n) by 'len' samples (out[0] follows x[n - 1]) with zero excitation of the lattice filter of k.
    inline void extrapolate(const std::vector<double>& k, const double* x, size_t n, double* out, size_t len)
    {
        const size_t p = k.size();
        if (p == 0) {
            for (size_t t = 0; t < len; ++t) out[t] = 0;
            return;
        }
        std::vector<double> bs(p, 0.0), f(p);       // bs[m] = b_m[t - 1], m = 0 .. p-1

        // Analysis lattice over the history; it is FIR of length p + 1, so the state is exact once p samples are in.
        for (size_t t = 0; t < n; ++t) {
            double fm = x[t], bcur = x[t];             // f_0[t], b_0[t]
            for (size_t m = 1; m <= p; ++m) {
                const double bprev = bs[m - 1];        // b_{m-1}[t-1]
                const double fnext = fm + k[m - 1] * bprev;
                const double bnext = bprev + k[m - 1] * fm;
                bs[m - 1] = bcur;                      // b_{m-1}[t]
                bcur = bnext;
                fm = fnext;
            }
        }

        // Synthesis lattice with f_p = 0.
        for (size_t t = 0; t < len; ++t) {
            double fm = 0;
            for (size_t m = p; m >= 1; --m) {
                fm -= k[m - 1] * bs[m - 1];            // f_{m-1}[t]
                f[m - 1] = fm;
            }
            const double y = std::isfinite(f[0]) ? f[0] : 0.0;
            for (size_t m = p - 1; m >= 1; --m)
                bs[m] = bs[m - 1] + k[m - 1] * f[m - 1];
            bs[0] = y;
            out[t] = y;
        }
    }

    //--------------------------------------------------------------------
    // Edge extrapolation for a resampler: extends audio before its start or after its end so the first and last
    // blocks see a continuation of the signal instead of a step from silence.  Shared by Ardftsrc.cpp and dBpoweramp.

    struct EdgeParams {
        size_t analysis;    // frames next to the edge the predictor is fitted to (~0.19 s)
        size_t order;       // predictor order (256 at 44.1 / 48 kHz)
        size_t predict;     // frames predicted (0.5 s); beyond this the extension is zero
    };

    inline EdgeParams edge_params(int sampleRate)
    {
        const double scale = sampleRate / 44100.0;
        EdgeParams p;
        p.analysis = static_cast<size_t>(std::clamp<long long>(std::llround(8192 * scale), 1024, 65536));
        p.order = static_cast<size_t>(std::clamp<long long>(std::llround(256 * scale), 32, 1024));
        p.predict = static_cast<size_t>((std::max)(sampleRate / 2, 1));
        return p;
    }

    // in: interleaved audio, 'frames' x 'channels', of which the analysis length nearest the edge is used
    //     (the first frames for pre = true, the last frames for pre = false).
    // Returns 'len' interleaved frames: for pre the last frame immediately precedes in[0], for post the first frame
    // immediately follows the last frame of in.  The prediction is faded to zero with a half-Hann window over
    // min(len, predict) frames (1 at the boundary) and is zero beyond; a channel whose prediction exceeds 4x the
    // peak of its analysed audio (or is not finite) is left at zero.
    inline std::vector<double> extrapolate_edge(const double* in, size_t frames, size_t channels, size_t len,
                                                bool pre, int sampleRate)
    {
        std::vector<double> ext(len * channels, 0.0);
        if (frames == 0 || channels == 0 || len == 0)
            return ext;

        const EdgeParams params = edge_params(sampleRate);
        const size_t n = (std::min)(frames, params.analysis);
        const size_t predictLen = (std::min)(len, params.predict);
        const double pi = 3.14159265358979323846;

        std::vector<double> ana(n), pred(predictLen);
        for (size_t c = 0; c < channels; ++c) {
            // analysis data in time order ending next to the gap: reversed for the start
            for (size_t i = 0; i < n; ++i)
                ana[i] = pre ? in[(n - 1 - i) * channels + c] : in[(frames - n + i) * channels + c];

            const std::vector<double> k = reflection_coefficients(ana.data(), n, params.order);
            extrapolate(k, ana.data(), n, pred.data(), predictLen);

            double peakIn = 0, peakOut = 0;
            for (double v : ana) peakIn = (std::max)(peakIn, std::abs(v));
            for (double v : pred) peakOut = (std::max)(peakOut, std::abs(v));
            if (!(peakOut <= 4.0 * peakIn))
                continue;

            for (size_t d = 0; d < predictLen; ++d) {       // d = distance from the boundary
                const double w = 0.5 * (1.0 + std::cos(pi * (static_cast<double>(d) + 0.5) / static_cast<double>(predictLen)));
                const size_t pos = pre ? len - 1 - d : d;
                ext[pos * channels + c] = pred[d] * w;
            }
        }
        return ext;
    }
}
