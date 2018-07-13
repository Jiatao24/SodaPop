#include <random>
#include "rng.h"

static pcg_extras::seed_seq_from<std::random_device> seedSource;
static std::uniform_real_distribution<double> uniformDist(0., 1.);


// Our uniform random bit generator is a global variable.
// By default (without setRngSeed), the seeding is random.
pcg32 g_rng(seedSource);


// Seed RNG using user-provided value.
//   pcg32::state_type is uint64_t.
void setRngSeed(pcg32::state_type seed)
{
    g_rng = pcg32(seed);
}


// Set RNG stream using user-provided value. Function has no effect if
//   setRngSeed is called after this function. pcg32::state_type is
//   uint64_t.
void setRngStream(pcg32::state_type stream)
{
    g_rng.set_stream(stream);
    g_rng();                    // burn a value
}


// Return random number in range [0, 1)
double randomNumber()
{
    return uniformDist(g_rng);
}
