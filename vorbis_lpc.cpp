#include "vorbis_lpc.h"
#include <vector>
#include <cmath>
#include <cstring>

namespace vorbis_lpc {

    static void from_data(sample_t* data, sample_t* lpci, size_t n, size_t m)
    {
        std::vector<double> aut(m + 1), lpc(m);
        double error;
        double epsilon;
        size_t i, j;

        /* autocorrelation, p+1 lag coefficients */
        j = m + 1;
        while (j--) {
            double d = 0; /* double needed for accumulator depth */
            for (i = j; i < n; i++)d += (double)data[i] * data[i - j];
            aut[j] = d;
        }

        /* Generate lpc coefficients from autocorr values */

        /* set our noise floor to about -100dB */
        error = aut[0] * (1. + 1e-10);
        epsilon = 1e-9 * aut[0] + 1e-10;

        for (i = 0; i < m; i++) {
            double r = -aut[i + 1];

            if (error < epsilon) {
                memset(lpc.data() + i, 0, (m - i) * sizeof(*lpc.data()));
                break;
            }

            /* Sum up this iteration's reflection coefficient; note that in
            Vorbis we don't save it.  If anyone wants to recycle this code
            and needs reflection coefficients, save the results of 'r' from
            each iteration. */

            for (j = 0; j < i; j++)r -= lpc[j] * aut[i - j];
            r /= error;

            /* Update LPC coefficients and total error */

            lpc[i] = r;
            for (j = 0; j < i / 2; j++) {
                double tmp = lpc[j];

                lpc[j] += r * lpc[i - 1 - j];
                lpc[i - 1 - j] += r * tmp;
            }
            if (i & 1)lpc[j] += lpc[j] * r;

            error *= 1. - r * r;

        }

#if 0 // disabled - degrades intermodulation test score
        /* slightly damp the filter */
        {
            double g = .99;
            double damp = g;
            for (j = 0; j < m; j++) {
                lpc[j] *= damp;
                damp *= g;
            }
        }
#endif
        for (j = 0; j < m; j++)lpci[j] = (float)lpc[j];
    }

    static void predict(sample_t* coeff, sample_t* prime, int m, sample_t* data, long n)
    {
        /* in: coeff[0...m-1] LPC coefficients
        prime[0...m-1] initial values (allocated size of n+m-1)
        out: data[0...n-1] data samples */

        long i, j, o, p;
        std::vector<sample_t> work(m + n, 0.f);

        if ( prime )
            for (i = 0; i < m; i++)
                work[i] = prime[i];

        for (i = 0; i < n; i++) {
            sample_t y = 0;
            o = i;
            p = m;
            for (j = 0; j < m; j++)
                y -= work[o++] * coeff[--p];

            if (!std::isnormal(y)) y = 0;

            data[i] = work[o] = y;
        }
    }


    void postextrapolate_chunk(sample_t* outpcm, const sample_t* inpcm, size_t ext_samples, size_t nch)
    {
        constexpr size_t order = 32;
        const size_t sc = ext_samples;

        if (sc <= order * 2) {
            for (size_t walk = 0; walk < sc * nch; ++walk) outpcm[walk] = 0;
            return;
        }

        std::vector<sample_t> work(sc * 2), lpc(order);
        for (size_t c = 0; c < nch; c++) {
            for (size_t i = 0; i < sc; i++) {
                work[i] = inpcm[i * nch + c];
            }
            vorbis_lpc::from_data(work.data(), lpc.data(), sc, order);
            vorbis_lpc::predict(lpc.data(), work.data() + sc - order, order, work.data() + sc, (long)sc);
            for (size_t i = 0; i < sc; i++) {
                outpcm[i * nch + c] = work[i + sc];
            }
        }
    }

    void preextrapolate_chunk(sample_t* outpcm, const sample_t* inpcm, size_t ext_samples, size_t nch)
    {
        constexpr size_t order = 32;
        const size_t sc = ext_samples;

        if (sc <= order * 2) { // safety
            for (size_t walk = 0; walk < sc * nch; ++walk) outpcm[walk] = 0;
            return;
        }

        std::vector<sample_t> work(sc * 2), lpc(order);

        for (size_t c = 0; c < nch; c++) {
            for (size_t i = 0; i < sc; i++) { // need to run the extrapolation in reverse!
                work[i] = inpcm[(sc - i - 1) * nch + c];
            }
            vorbis_lpc::from_data(work.data(), lpc.data(), sc, order);
            vorbis_lpc::predict(lpc.data(), work.data() + sc - order, order, work.data() + sc, (long)sc);
            for (size_t i = 0; i < sc; i++) {
                outpcm[(sc - i - 1) * nch + c] = work[i + sc];
            }
        }
    }
    std::vector<sample_t> postextrapolate_chunk(const sample_t* inpcm, size_t ext_samples, size_t nch) {
        std::vector<sample_t> ret; ret.resize(ext_samples * nch);
        postextrapolate_chunk(ret.data(), inpcm, ext_samples, nch);
        return ret;
    }
    std::vector<sample_t> preextrapolate_chunk(const sample_t* inpcm, size_t ext_samples, size_t nch) {
        std::vector<sample_t> ret; ret.resize(ext_samples * nch);
        preextrapolate_chunk(ret.data(), inpcm, ext_samples, nch);
        return ret;
    }

}
