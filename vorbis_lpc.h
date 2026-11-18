#pragma once
#include <vector>
#include <cstddef>

namespace vorbis_lpc {
	typedef double sample_t;

	void postextrapolate_chunk(sample_t* outpcm, const sample_t* inpcm, size_t ext_samples, size_t nch);
	void preextrapolate_chunk(sample_t* outpcm, const sample_t* inpcm, size_t ext_samples, size_t nch);
	std::vector<sample_t> postextrapolate_chunk(const sample_t* inpcm, size_t ext_samples, size_t nch);
	std::vector<sample_t> preextrapolate_chunk(const sample_t* inpcm, size_t ext_samples, size_t nch);
}
