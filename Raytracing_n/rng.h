#ifndef RNG_H
#define RNG_H

#include <cmath>
#include <stdint.h>
static const float FloatOneMinusEpsilon = 0.99999994;

static const float OneMinusEpsilon = FloatOneMinusEpsilon;

#define PCG32_DEFAULT_STATE 0x853c49e6748fea9bULL
#define PCG32_DEFAULT_STREAM 0xda3e39cb94b95bdbULL
#define PCG32_MULT 0x5851f42d4c957f2dULL

class rng {
public:
	rng();

	uint32_t UniformUInt32();

	float UniformFloat() {
		return std::fmin(OneMinusEpsilon, float(UniformUInt32() * 2.3283064365386963e-10f));
	}
private:
	uint64_t state, inc;
};

inline rng::rng() : state(PCG32_DEFAULT_STATE), inc(PCG32_DEFAULT_STREAM) {}

inline uint32_t rng::UniformUInt32() {
	uint64_t oldstate = state;
	state = oldstate * PCG32_MULT + inc;
	uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
	uint32_t rot = (uint32_t)(oldstate >> 59u);
	return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}
#endif