#ifndef OCCJSDM_RNG_H
#define OCCJSDM_RNG_H

// Serial, R-seeded random number generation.
//
// Every RNG-bearing C++ sampling body runs directly on the main thread,
// including helpers with historical _parallel or _TS names. Deterministic
// probability workers may still use TBB. Never call this engine, R's RNG or
// Armadillo's R-backed RNG from a worker thread.
//
// TBB workers cannot be assigned independent streams using
// omp_get_thread_num(): outside an OpenMP region each reports zero. The old
// per-thread engines therefore duplicated streams as well as making draw
// order depend on scheduling. One advancing stream avoids both defects.
//
// runOccJSDM() draws one base seed from R per fit. Chains and sampling blocks
// continue that stream rather than re-seeding by species or iteration. The
// generation counter resets the engine and normal cache when a new fit is
// seeded, but is deliberately not itself seed material.

#include <random>

// Shared by all translation units and accessed only on the main thread.
inline unsigned int& occjsdm_rng_base_seed() {
  static unsigned int seed = 12345u;
  return seed;
}

inline unsigned int& occjsdm_rng_generation() {
  static unsigned int generation = 1u;   // first use seeds
  return generation;
}

// One advancing mt19937 for the serial sampler.
inline std::mt19937& get_rng() {
  static std::mt19937 rng;
  static unsigned int seeded_at = 0u;   // 0 == never seeded

  const unsigned int current = occjsdm_rng_generation();
  if (seeded_at != current) {
    // Retain the former main-thread seed mapping so single-thread draws
    // remain unchanged. A generation change only requests a reset.
    std::seed_seq seq{ occjsdm_rng_base_seed(), 0u };
    rng.seed(seq);
    seeded_at = current;
  }
  return rng;
}

inline double runif() {
  // std::uniform_real_distribution is stateless, so unlike rnorm() below it
  // needs no reset when the engine is re-seeded.
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(get_rng());
}

inline double rnorm() {
  static std::normal_distribution<double> dist(0.0, 1.0);
  static unsigned int dist_reset_at = 0u;

  std::mt19937& rng = get_rng();   // re-seeds if R has set a new base seed

  // Re-seeding the engine is not enough. std::normal_distribution generates
  // deviates in pairs (Box-Muller) and caches the second one; that cache
  // survives a call to rng.seed(). Without the reset below, the first rnorm()
  // of a fit can return a value left over from the *previous* fit, so whether
  // a fit reproduces depends on the parity of the normal draws that happened
  // to precede it in the session. That is exactly what made this reproducible
  // when run alone but not when run after other tests in the same process.
  const unsigned int current = occjsdm_rng_generation();
  if (dist_reset_at != current) {
    dist.reset();
    dist_reset_at = current;
  }
  return dist(rng);
}

#endif // OCCJSDM_RNG_H
