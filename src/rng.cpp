#include <random>
#include "rng.hpp"

pcg32 rng;

// Seed RNG using std::random_device.
void initialize_rng()
{
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    rng = pcg32(seed_source);
}

// Seed RNG using user-provided value.
//   pcg32::state_type is uint64_t.
void initialize_rng(pcg32::state_type seed)
{
    rng = pcg32(seed);
}
