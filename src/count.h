#ifndef OCCJSDM_COUNT_H
#define OCCJSDM_COUNT_H

// Deterministic formulas for the count model (negative binomial fitted by
// Polya-Gamma augmentation). Design:
// dev/superpowers/specs/2026-11-12-count-model-design.md, "How it is fitted".
//
// Nothing here draws a random number, so these may be called from any thread.
// The sampler that uses them, samplePGcount_cpp(), is in jsdm.cpp.

#include <algorithm>
#include <cmath>
#include <R_ext/Arith.h>

// Mean of a Polya-Gamma variable with shape b at predictor z:
// b * tanh(|z| / 2) / (2 |z|). Near zero the ratio is 0 / 0, so the first two
// terms of its series are used instead.
inline double pg_mean(double b, double z) {
  const double a = std::fabs(z);
  if (a < 1e-3) {
    return b * (0.25 - a * a / 48.0);
  }
  return b * std::tanh(0.5 * a) / (2.0 * a);
}

// Variance of the same variable: b * (sinh|z| - |z|) / (4 |z|^3 cosh^2(|z|/2)),
// written with t = tanh(|z| / 2) so that nothing overflows at a large
// predictor. Near zero the numerator cancels, so the series is used.
inline double pg_var(double b, double z) {
  const double a = std::fabs(z);
  if (a < 1e-2) {
    return b * (1.0 / 24.0 - a * a / 120.0);
  }
  const double t = std::tanh(0.5 * a);
  return b * (2.0 * t - a * (1.0 - t * t)) / (4.0 * a * a * a);
}

// log(exp(a) + exp(b)) without overflow.
inline double count_logsumexp(double a, double b) {
  return std::max(a, b) + std::log1p(std::exp(-std::fabs(a - b)));
}

// Log probability of the count y under a negative binomial with mean
// exp(log_mean) and size r, so that the variance is mean + mean^2 / r. An
// effort offset is part of log_mean. A missing count gives NA; a count that
// is negative or not a whole number has probability zero.
inline double nb_logpmf(double y, double log_mean, double r) {
  if (ISNAN(y)) {
    return NA_REAL;
  }
  if (!R_FINITE(y) || y < 0.0 || y != std::floor(y)) {
    return R_NegInf;
  }
  const double log_r = std::log(r);
  return std::lgamma(y + r) - std::lgamma(r) - std::lgamma(y + 1.0) +
    r * log_r + y * log_mean -
    (y + r) * count_logsumexp(log_mean, log_r);
}

#endif // OCCJSDM_COUNT_H
