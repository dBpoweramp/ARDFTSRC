// Ardftsrc.cpp : ARDFTSRC - Licensed under WTFPL
//   Created by mycroft @ hydrogenaudio.org
//   C port by the team @ dBpoweramp.com
//
// Command line front end for the CARDFTSRC class (CARDFTSRC.h), running on the header-only pocketfft transform
// (FFT-pocketfft.h).  Named modes select the presets below; "_f128" on the end of any mode name runs the whole
// resampler (taper, FFTs, overlap-add) in double-double precision via boost::multiprecision::cpp_double_double.

#ifndef _SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING
#define _SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING   // MSVC STL4037: std::complex<cpp_double_double> is used deliberately
#endif

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <optional>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cwctype>
#include <algorithm>

#include <boost/multiprecision/cpp_double_fp.hpp>   // Boost >= 1.89 (header only)
#include "AudioFile.h"                              // WAV / AIFF I/O (modified copy in this tree: adds 64-bit float)
#include "CARDFTSRC.h"
#include "vorbis_lpc.h"
#include "burg_lpc.h"

// "F128" precision: a double-double, i.e. an unevaluated sum of two doubles (106-bit significand, ~32 digits).
// Not IEEE binary128 (113 bits), but close, and it builds on every compiler including MSVC.
using f128_t = boost::multiprecision::cpp_double_double;

//--------------------------------------------------------------------
// Taper selection (see CARDFTSRC.h):  kaiserBeta > 0 -> Kaiser-kernel taper (recommended; its impulse response
// collapses ~exp(-beta) beyond t = beta / (pi * W), W = transition width = (1 - bandwidth) * output Nyquist),
// otherwise Betacdftaper > 0 -> beta-CDF taper, otherwise rational_cosine_taper.
// straddle = true centres the transition band on the output Nyquist frequency and folds the upper half back
// (aliasing the source band [Nyquist, Nyquist + W/2] into [Nyquist - W/2, Nyquist], as SoX does), which puts the
// -3 dB point at Nyquist itself.  Kaiser taper and down-sampling only; ignored otherwise.
struct STArdMode { const wchar_t* Name; int quality; double bandwidth; double phase; // phase: 0 = linear .. 1 = minimum phase
                    int Betacdftaper; double kaiserBeta; bool straddle; bool useF128; };

static const STArdMode ArdModes[] = {
    //  name                          quality   bandwidth phase  betaCDF kaiser straddle F128   ** F128 is set when _f128 is on the end of the name
    {L"fast",                          2048,    0.96,   0.0,   0,      0.0,  false, false},
    {L"standard",                      8192,    0.98,   0.0,   0,      0.0,  false, false},
    {L"",                              65536,   0.9964, 0.0,   0,     32.0,  false, false},       // default
    {L"high",                          262144,  0.9967, 0.0,   0,     32.0,  false, false},
    {L"insane",                        262144,  0.997,  0.0,   0,     32.0,  false, false},       // W = 35 Hz, -3 dB 19 Hz below Nyquist
    {L"max bandwidth",                 262144,  0.998,  0.0,   0,     40.0,  false, false},       // W = 44 Hz; noise-floor bandwidth 100 %, impulse frequency 100 %
    {L"min phase",                     65536,   0.99,   0.3,   0,     32.0,  false, false},
    {L"max bandwidth, min phase",      262144,  0.999,  0.3,   0,     32.0,  false, false},
};

//--------------------------------------------------------------------
static bool EqualsNoCase(const std::wstring& a, const std::wstring& b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
        [](wchar_t x, wchar_t y) { return std::towlower(x) == std::towlower(y); });
}

//--------------------------------------------------------------------
// Returns a copy of the matching mode, with useF128 set when the name ends in "_f128" (any name can carry it,
// "_f128" on its own is the default mode).  A copy rather than a pointer into the table, so concurrent lookups
// cannot see each other's F128 flag.
static std::optional<STArdMode> FindArdMode(const wchar_t* Mode) {
    std::wstring MatchName(Mode ? Mode : L"");
    const std::wstring Suffix = L"_f128";
    bool IsF128 = false;
    if (MatchName.size() >= Suffix.size() && EqualsNoCase(MatchName.substr(MatchName.size() - Suffix.size()), Suffix)) {
        IsF128 = true;
        MatchName.erase(MatchName.size() - Suffix.size());
    }
    for (const auto& walk : ArdModes) {
        if (EqualsNoCase(MatchName, walk.Name)) {
            STArdMode found = walk;
            found.useF128 = IsF128;
            return found;
        }
    }
    return std::nullopt;
}

//--------------------------------------------------------------------
static std::wstring Widen(const std::string& s) { return std::wstring(s.begin(), s.end()); }   // mode names are ASCII
static std::string Narrow(const wchar_t* s) { std::string r; for (; s && *s; ++s) r += (char)*s; return r; }

static bool ParseLong(const char* s, long& out) {
    char* end = nullptr; out = std::strtol(s, &end, 10);
    return end != s && *end == 0;
}
static bool ParseDouble(const char* s, double& out) {
    char* end = nullptr; out = std::strtod(s, &end);
    return end != s && *end == 0;
}
static int64_t CeilDiv(int64_t a, int64_t b) { return (a + b - 1) / b; }

//--------------------------------------------------------------------
static void PrintUsage() {
    std::cout <<
        "\nUsage:  ardftsrc <in.wav> <out.wav> <out frequency> <out bitdepth> [mode] [options]\n"
        "        (out bitdepth: 8, 16, 24, 32 = 32-bit float, 64 = 64-bit float)\n"
        "        ardftsrc <in.wav> <out.wav> <out frequency> <out bitdepth> <quality> [bandwidth] [options]\n"
        "\nModes (append _f128 to any mode for double-double precision, e.g. \"high_f128\"; \"_f128\" alone = default):\n";
    for (const auto& m : ArdModes) {
        const std::string name = *m.Name ? "\"" + Narrow(m.Name) + "\"" : "(none)  default";
        std::cout << "  " << std::left << std::setw(28) << name << std::right
                  << " quality " << std::setw(6) << m.quality << "  bandwidth " << std::setw(6) << m.bandwidth;
        if (m.phase != 0) std::cout << "  phase " << m.phase;
        std::cout << "\n";
    }
    std::cout <<
        "\nOptions:\n"
        "  --f128         double-double precision (same as the _f128 suffix)\n"
        "  --no-lpc       zero-pad the start and end of the file instead of LPC extrapolation\n"
        "  --vorbis-lpc   previous edge extrapolation (vorbis_lpc, order 32) instead of Burg\n"
        "\nExamples:  ardftsrc \"c:\\in.wav\" \"c:\\out.wav\" 48000 24 high\n"
        "           ardftsrc \"c:\\in.wav\" \"c:\\out.wav\" 44100 24 \"max bandwidth, min phase_f128\"\n"
        "           ardftsrc \"c:\\in.wav\" \"c:\\out.wav\" 48000 24 8192 0.99      (numeric quality / bandwidth)\n";
}

//--------------------------------------------------------------------
enum class EdgeMode { Burg, Vorbis, Zero };

// Builds 'len' samples per channel (interleaved) of extrapolation before the start (pre) or after the end (post) of
// the file, so the first and last blocks see a continuation of the signal instead of a step from silence.
//
// Burg (default): an order ~256 predictor (scaled with the sample rate) from the ~0.19 s of audio next to the edge,
// run as a lattice filter for up to 0.5 s and faded out with a half-Hann window (1 at the file boundary).
// Measured against resampling the uncut signal, the error in the first/last 100 ms of output drops by ~10-15 dB on
// music and to the numerical floor on steady tones, compared with the vorbis_lpc method.
// Vorbis: the previous method - vorbis_lpc order 32, analysis and prediction over the whole edge, faded over it.
// Either way a channel whose prediction exceeds 4x the peak of the analysed audio is left at zero.
static std::vector<double> EdgeExtrapolation(const AudioFile<double>::AudioBuffer& samples, int64_t numSamples,
                                             int64_t len, bool pre, int sampleRate, EdgeMode method) {
    const size_t nch = samples.size();
    std::vector<double> ext(static_cast<size_t>(len) * nch, 0.0);
    if (method == EdgeMode::Zero || numSamples <= 0 || nch == 0 || len <= 0)
        return ext;

    const double pi = 3.14159265358979323846;
    const double rateScale = sampleRate / 44100.0;
    const int64_t analysisLen = method == EdgeMode::Burg
        ? std::min(numSamples, std::clamp<int64_t>(std::llround(8192 * rateScale), 1024, 65536))
        : std::min(len, numSamples);
    const int64_t predictLen = method == EdgeMode::Burg ? std::min<int64_t>(len, std::max(sampleRate / 2, 1)) : analysisLen;
    const int64_t fadeLen = method == EdgeMode::Burg ? predictLen : len;
    const size_t order = static_cast<size_t>(std::clamp<int64_t>(std::llround(256 * rateScale), 32, 1024));

    // Analysis data in time order, ending next to the gap: reversed for the start of the file.
    std::vector<double> ana(static_cast<size_t>(analysisLen)), pred(static_cast<size_t>(predictLen));
    for (size_t c = 0; c < nch; ++c) {
        for (int64_t i = 0; i < analysisLen; ++i)
            ana[static_cast<size_t>(i)] = pre ? samples[c][static_cast<size_t>(analysisLen - 1 - i)]
                                              : samples[c][static_cast<size_t>(numSamples - analysisLen + i)];

        if (method == EdgeMode::Burg) {
            const auto k = burg_lpc::reflection_coefficients(ana.data(), ana.size(), order);
            burg_lpc::extrapolate(k, ana.data(), ana.size(), pred.data(), pred.size());
        }
        else    // vorbis_lpc::postextrapolate_chunk continues its input, which is already ordered towards the gap
            vorbis_lpc::postextrapolate_chunk(pred.data(), ana.data(), ana.size(), 1);

        double peakIn = 0, peakOut = 0;
        for (double v : ana) peakIn = std::max(peakIn, std::abs(v));
        for (double v : pred) peakOut = std::max(peakOut, std::abs(v));
        if (!(peakOut <= 4.0 * peakIn))                 // unstable prediction (or NaN): leave this channel zero
            continue;

        for (int64_t d = 0; d < predictLen; ++d) {      // d = distance from the file boundary
            const double w = 0.5 * (1.0 + std::cos(pi * (static_cast<double>(d) + 0.5) / static_cast<double>(fadeLen)));
            const int64_t pos = pre ? len - 1 - d : d;
            ext[static_cast<size_t>(pos) * nch + c] = pred[static_cast<size_t>(d)] * w;
        }
    }
    return ext;
}

//--------------------------------------------------------------------
// Resamples the whole file.  The class is fed a stream  [edge | audio | edge | zeros...]  in blocks of InChunkSize(),
// where 'edge' is half a block (FinalBlocksZeroPadSample()) of extrapolation (or zeros).  The output delay
// (FirstBlocksSampleDelay()) plus the resampled lead-in is dropped, so the output is aligned with the input and
// exactly round(numSamples * outRate / inRate) samples long.
template<typename scalar_t>
static bool Resample(const AudioFile<double>& in, int outRate, const STArdMode& mode, EdgeMode edgeMode,
                     AudioFile<double>::AudioBuffer& out) {
    const int inRate = static_cast<int>(in.getSampleRate());
    const int nch = in.getNumChannels();
    const int64_t numIn = in.getNumSamplesPerChannel();

    CARDFTSRC<CFFT_pocketfft<scalar_t>> src(inRate, outRate, nch, static_cast<unsigned>(mode.quality), mode.bandwidth,
                                            mode.phase, mode.Betacdftaper, mode.kaiserBeta, mode.straddle);

    const int64_t inChunk = src.InChunkSize(), outChunk = src.OutChunkSize();
    const int64_t edge = src.FinalBlocksZeroPadSample();                 // half an input block
    const int64_t edgeOut = edge * outChunk / inChunk;                    // the same span at the output rate (exact)
    const int64_t skip = src.FirstBlocksSampleDelay() + edgeOut;          // output samples dropped at the start
    const int64_t numOut = (numIn * outRate + inRate / 2) / inRate;
    const int64_t numChunks = std::max(CeilDiv(edge + numIn + edge, inChunk), CeilDiv(skip + numOut, outChunk));

    std::cout << "Block size: " << inChunk << " in / " << outChunk << " out samples   FFT: "
              << 2 * inChunk << " / " << 2 * outChunk << "\n\n";

    const std::vector<double> head = EdgeExtrapolation(in.samples, numIn, edge, true, inRate, edgeMode);
    const std::vector<double> tail = EdgeExtrapolation(in.samples, numIn, edge, false, inRate, edgeMode);

    // Stream sample s (per channel): [0, edge) head, [edge, edge + numIn) audio, then tail, then zeros.
    auto streamSample = [&](int64_t s, int c) -> double {
        if (s < edge) return head[static_cast<size_t>(s) * nch + c];
        s -= edge;
        if (s < numIn) return in.samples[c][static_cast<size_t>(s)];
        s -= numIn;
        if (s < edge) return tail[static_cast<size_t>(s) * nch + c];
        return 0.0;
    };

    out.assign(nch, std::vector<double>(static_cast<size_t>(numOut), 0.0));

    for (int64_t chunk = 0; chunk < numChunks; ++chunk) {
        std::cout << "Resampling " << std::fixed << std::setprecision(2) << (100.0 * chunk / numChunks) << " %\r" << std::flush;
        const int64_t inBase = chunk * inChunk;
        const int64_t outBase = chunk * outChunk - skip;
        src.process_ex(
            [&](size_t idx) { return streamSample(inBase + static_cast<int64_t>(idx / nch), static_cast<int>(idx % nch)); },
            [&](size_t idx, double v) {
                const int64_t o = outBase + static_cast<int64_t>(idx / nch);
                if (o >= 0 && o < numOut) out[idx % nch][static_cast<size_t>(o)] = v;
            });
    }
    std::cout << "Resampling 100.00 %" << std::endl;
    return numOut > 0;
}

//--------------------------------------------------------------------
int main(int argc, char* argv[])
{
    std::cout << ("ARDFTSRC - Licensed under WTFPL\n  Created by mycroft @ hydrogenaudio.org\n  C port by the team @ dBpoweramp.com\n");
    if (argc < 5) {
        PrintUsage();
        return 0;
    }

    const std::string input_file = argv[1];
    const std::string output_file = argv[2];
    long output_samplerate = 0, output_bitdepth = 0;
    if (!ParseLong(argv[3], output_samplerate) || output_samplerate <= 0) {
        std::cerr << "Error: invalid output frequency " << argv[3] << std::endl;
        return 1;
    }
    if (!ParseLong(argv[4], output_bitdepth) || (output_bitdepth != 8 && output_bitdepth != 16 && output_bitdepth != 24 && output_bitdepth != 32 && output_bitdepth != 64)) {
        std::cerr << "Error: output bitdepth must be 8, 16, 24, 32 (float) or 64 (float)" << std::endl;
        return 1;
    }

    // Remaining arguments: options, then either a mode name or numeric quality [bandwidth].
    bool forceF128 = false;
    EdgeMode edgeMode = EdgeMode::Burg;
    std::vector<std::string> rest;
    for (int i = 5; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--f128") forceF128 = true;
        else if (a == "--no-lpc") edgeMode = EdgeMode::Zero;
        else if (a == "--vorbis-lpc") edgeMode = EdgeMode::Vorbis;
        else rest.push_back(a);
    }

    STArdMode mode = *FindArdMode(L"");
    std::string modeLabel = "default";
    long quality = 0;
    if (!rest.empty() && ParseLong(rest[0].c_str(), quality)) {
        // Numeric form (as in earlier versions): quality [bandwidth], class default (rational cosine) taper.
        double bandwidth = 0.95;
        if (rest.size() >= 2 && !ParseDouble(rest[1].c_str(), bandwidth)) {
            std::cerr << "Error: invalid bandwidth " << rest[1] << std::endl;
            return 1;
        }
        if (quality <= 0 || !(bandwidth > 0.0 && bandwidth < 1.0)) {
            std::cerr << "Error: quality must be > 0 and bandwidth between 0 and 1" << std::endl;
            return 1;
        }
        mode = STArdMode{ L"custom", (int)quality, bandwidth, 0.0, 0, 0.0, false, false };
        modeLabel = "custom";
    }
    else if (!rest.empty()) {
        const auto found = FindArdMode(Widen(rest[0]).c_str());
        if (!found) {
            std::cerr << "Error: unknown mode \"" << rest[0] << "\"" << std::endl;
            PrintUsage();
            return 1;
        }
        mode = *found;
        modeLabel = *mode.Name ? Narrow(mode.Name) : "default";
    }
    mode.useF128 = mode.useF128 || forceF128;

    // ######################################

    AudioFile<double> inputAudioFile;
    inputAudioFile.shouldLogErrorsToConsole(true);
    if (!inputAudioFile.load(input_file)) {
        std::cerr << "Error: Could not load input file " << input_file << std::endl;
        return 1;
    }
    const int in_samplerate = (int)inputAudioFile.getSampleRate();
    const int num_channels = inputAudioFile.getNumChannels();
    if (in_samplerate <= 0 || num_channels <= 0 || inputAudioFile.getNumSamplesPerChannel() <= 0) {
        std::cerr << "Error: input file has no audio" << std::endl;
        return 1;
    }

    std::cout << "\nSource: " << input_file
              << "\n  Sample Rate:" << in_samplerate
              << "\n  Bit Depth:" << inputAudioFile.getBitDepth()
              << "\n  Channels:" << num_channels;

    std::cout << "\n\nOutput: " << output_file
              << "\n  Sample Rate:" << output_samplerate
              << "\n  Bit Depth:" << output_bitdepth
              << "\n  Channels:" << num_channels;

    std::cout << "\n\nParameters  Mode: " << modeLabel
              << "   Quality: " << mode.quality
              << "   Bandwidth: " << mode.bandwidth
              << "   Phase: " << mode.phase
              << "\n            Taper: ";
    if (mode.kaiserBeta > 0) std::cout << "Kaiser, beta " << mode.kaiserBeta << (mode.straddle ? ", straddling Nyquist" : "");
    else if (mode.Betacdftaper > 0) std::cout << "beta-CDF, order " << mode.Betacdftaper;
    else std::cout << "rational cosine";
    std::cout << "\n            Precision: " << (mode.useF128 ? "double-double (boost cpp_double_double)" : "double")
              << "   Edges: " << (edgeMode == EdgeMode::Burg ? "Burg LPC extrapolation" : edgeMode == EdgeMode::Vorbis ? "vorbis LPC extrapolation" : "zero padding") << "\n";

    // ######################################

    const auto t0 = std::chrono::steady_clock::now();
    AudioFile<double>::AudioBuffer outBuffer;
    const bool ok = mode.useF128
        ? Resample<f128_t>(inputAudioFile, (int)output_samplerate, mode, edgeMode, outBuffer)
        : Resample<double>(inputAudioFile, (int)output_samplerate, mode, edgeMode, outBuffer);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (!ok) {
        std::cerr << "No output was generated" << std::endl;
        return 1;
    }
    std::cout << "Time: " << std::fixed << std::setprecision(2) << seconds << " s\n";

    AudioFile<double> outputAudioFile;
    outputAudioFile.shouldLogErrorsToConsole(true);
    outputAudioFile.setAudioBuffer(outBuffer);
    outputAudioFile.setSampleRate((uint32_t)output_samplerate);
    outputAudioFile.setBitDepth((int)output_bitdepth);
    if (!outputAudioFile.save(output_file, AudioFileFormat::Wave)) {
        std::cerr << "Error: Could not save output file " << output_file << std::endl;
        return 1;
    }

    return 0;
}
