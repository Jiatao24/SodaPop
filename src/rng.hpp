#ifndef RNG_HPP
#define RNG_HPP

#include "pcg/pcg_random.hpp"

extern pcg32 rng;

// Seed RNG using std::random_device.
void initialize_rng();

// Seed RNG using user-provided value.
//   pcg32::state_type is uint64_t.
void initialize_rng(pcg32::state_type seed);

#endif
