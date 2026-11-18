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
#include <cmath>

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
}
