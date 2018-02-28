#include <random>
#include "rng.hpp"

pcg32 g_rng;
static std::uniform_real_distribution<double> uniform_dist(0., 1.);

// Seed RNG using std::random_device.
void initialize_rng()
{
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    g_rng = pcg32(seed_source);
}

// Seed RNG using user-provided value.
//   pcg32::state_type is uint64_t.
void initialize_rng(pcg32::state_type seed)
{
    g_rng = pcg32(seed);
}

// Return random number in range [0, 1)
double random_number()
{
    return uniform_dist(g_rng);
}
