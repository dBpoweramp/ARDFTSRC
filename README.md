# ARDFTSRC
ARDFTSRC Audio Resampler, by mycroft, port to C and Windows by dBpoweramp.com

The resampler itself is the header-only class `CARDFTSRC` (`CARDFTSRC.h`); `Ardftsrc.cpp` is a command line front end for it.
FFTs are done by the bundled, header-only [pocketfft](https://github.com/mreineck/pocketfft), and an optional double-double
precision path uses Boost.Multiprecision.

## Usage

    ardftsrc.exe "c:\in.wav" "c:\out.wav" 48000 24 [mode] [--f128] [--no-lpc | --vorbis-lpc]
    ardftsrc.exe "c:\in.wav" "c:\out.wav" 48000 24 8192 0.99        (numeric quality [bandwidth], as in earlier versions)

`48000` is the output frequency and `24` the output bit depth: 8, 16, 24, 32 (32-bit float) or 64 (64-bit float).
64-bit float output keeps the resampler's result unrounded, which is the only output format in which the `_f128`
precision can make a measurable difference. Input can be WAV or AIFF at 8, 16, 24 or 32-bit, or 64-bit float
(WAV, including WAVE_FORMAT_EXTENSIBLE, or AIFF-C `fl64`).

| Mode | Quality | Bandwidth | Phase | Taper |
|---|---|---|---|---|
| `fast` | 2048 | 0.96 | linear | rational cosine |
| `standard` | 8192 | 0.98 | linear | rational cosine |
| *(none / `""`)* — default | 65536 | 0.9964 | linear | Kaiser, beta 32 |
| `high` | 262144 | 0.9967 | linear | Kaiser, beta 32 |
| `insane` | 262144 | 0.997 | linear | Kaiser, beta 32 |
| `max bandwidth` | 262144 | 0.998 | linear | Kaiser, beta 40 |
| `min phase` | 65536 | 0.99 | 0.3 | Kaiser, beta 32 |
| `max bandwidth, min phase` | 262144 | 0.999 | 0.3 | Kaiser, beta 32 |

Mode names are case-insensitive; names containing spaces need quotes, e.g. `"max bandwidth, min phase"`.

**`_f128`** on the end of any mode name (`high_f128`, `"min phase_f128"`, or just `_f128` for the default mode), or the
`--f128` option, runs the whole resampler — taper, FFTs, overlap-add — in double-double precision
(`boost::multiprecision::cpp_double_double`, 106-bit significand, about 32 significant digits; not IEEE binary128).
This takes the numerical error from about -310 dB (double) to about -620 dB; it is roughly 20x slower, and the
difference is far below anything a 24-bit or 32-bit float output file can hold (see 64-bit float output above).

**Edges**: the start and end of the file are extended by linear prediction before resampling, so the first and last
blocks see a continuation of the signal rather than a step from silence. The predictor (`burg_lpc.h`) is fitted with
Burg's method to the ~0.19 s next to each edge, order 256 at 44.1/48 kHz (scaled with the sample rate), and run as a
lattice filter for up to 0.5 s with a half-Hann fade. In tests that cut segments out of music and compared the result
with resampling the uncut file, Burg lowered the error in the first and last 100 ms by typically 10-15 dB versus the
earlier `vorbis_lpc` order-32 method (zero padding is around -55 dB), and on steady tones it reaches the numerical
floor. Unpredictable, noise-like material gains nothing from any method.
`--vorbis-lpc` selects the earlier method (output identical to the previous version), `--no-lpc` zero-pads.

The output is delay-compensated and exactly `round(input length * out rate / in rate)` samples long.

## Requirements (add to project folder)

    \boost\         Boost 1.89 or later, headers only: copy the 'boost' folder of a Boost distribution here
                    (so that \boost\multiprecision\cpp_double_fp.hpp exists), or set the BOOST_ROOT environment variable

C++20 is required (Visual Studio 2022, toolset v143). AudioFile is now part of the source tree (a modified copy, see
below), and Eigen and FFTW are no longer needed.

Do not build with fast floating-point math (`/fp:fast`, `-ffast-math`): double-double arithmetic depends on exact IEEE
rounding and silently degrades to double precision. The project sets `/fp:precise`.

## Files

| File | |
|---|---|
| `Ardftsrc.cpp` | command line front end, mode table |
| `AudioFile.h` | WAV / AIFF reading and writing: modified copy of [AudioFile](https://github.com/adamstark/AudioFile) (MIT licence, notice in the file) |
| `CARDFTSRC.h` | the resampler class, templated on the FFT backend and its scalar type |
| `FFT-pocketfft.h`, `pocketfft_hdronly.h` | pocketfft backend (the default) |
| `burg_lpc.h` | Burg / lattice LPC extrapolation for the file edges (default) |
| `vorbis_lpc.h`, `vorbis_lpc.cpp` | earlier autocorrelation LPC extrapolation (`--vorbis-lpc`) |
| `CFFT-fftw3.h`, `CFFT-fftw3.cpp` | optional FFTW3 backend (`CARDFTSRC<CFFT_fftw3>`, double only); excluded from the build by default |

`AudioFile.h` is based on upstream commit 4554279 (October 2025), with changes marked `ARDFTSRC:`. It reads and writes
64-bit float WAV and reads 64-bit float AIFF-C. It takes the sample format of WAVE_FORMAT_EXTENSIBLE files from their
SubFormat (upstream reads extensible 32-bit float as integer PCM). It handles files between 2 and 4 GB. It skips the pad
byte after odd-sized chunks (upstream fails on such files). It writes a `fact` chunk in float WAV files. AudioFile's own
test suite passes unchanged.

`pocketfft_hdronly.h` carries a small patch (marked `ARDFTSRC patch`) so that it works with class types such as
`cpp_double_double`: the twiddle factors use argument-dependent `sin`/`cos` and a full-precision pi, and the hard-coded
butterfly constants are parsed from their decimal text for such types instead of being rounded to `long double`.
Float, double and long double results are unchanged.

## Building with GCC / Clang

    g++ -std=c++20 -O2 -I. -I<path to boost> Ardftsrc.cpp vorbis_lpc.cpp -o ardftsrc
