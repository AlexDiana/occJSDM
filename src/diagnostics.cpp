#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Per-chain moments and effective sample size for a block of posterior draws,
// one parameter at a time, for computeDiagnostics() (R/diagnostics.R).
//
// The statistics are the ones coda computes, so that the printed convergence
// summary is unchanged: the effective sample size of one chain is
// coda::effectiveSize(), that is n * var(x) / spectrum0.ar(x), where
// spectrum0.ar() fits an autoregression by Yule-Walker with the order chosen
// by AIC (stats::ar()) and returns its spectral density at frequency zero.
// Nothing here is shared between parameters, so time and memory grow with the
// number of parameters and not with its square.

namespace {

// Dot product with four accumulators. A single running sum is a dependency
// chain the compiler may not reorder, which makes the autocovariances (most of
// the work) several times slower than they need to be.
inline double dot(const double* a, const double* b, int n) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    s0 += a[i] * b[i];
    s1 += a[i + 1] * b[i + 1];
    s2 += a[i + 2] * b[i + 2];
    s3 += a[i + 3] * b[i + 3];
  }
  for (; i < n; ++i) s0 += a[i] * b[i];
  return (s0 + s1) + (s2 + s3);
}

struct ArWork {
  std::vector<double> r, v, phi, tmp, phisum;
  explicit ArWork(int order_max)
      : r(order_max + 1), v(order_max + 1), phi(order_max + 1),
        tmp(order_max + 1), phisum(order_max + 1) {}
};

// Effective sample size of one chain, given the centred series xc and its
// sample variance. Returns NA when the series is one that this code does not
// reproduce coda on with confidence (see below); the R caller recomputes
// those with coda itself.
double chain_ess(const double* xc, int n, double var, int order_max,
                 ArWork& w) {
  if (n < 2 || order_max < 1) return NA_REAL;

  // spectrum0.ar() first regresses the chain on the iteration index and
  // calls the spectral density zero when the residuals have no spread:
  //   identical(all.equal(sd(residuals(lm(x ~ z))), 0), TRUE)
  // which is sd <= sqrt(.Machine$double.eps). A chain that never moved is
  // the case that matters; an exactly linear one is covered too.
  const double centre = 0.5 * (n - 1);
  double sxz = 0.0, szz = 0.0;
  for (int i = 0; i < n; ++i) {
    const double z = i - centre;
    sxz += z * xc[i];
    szz += z * z;
  }
  const double slope = sxz / szz;
  double rsum = 0.0;
  for (int i = 0; i < n; ++i) rsum += xc[i] - slope * (i - centre);
  const double rmean = rsum / n;
  double rss = 0.0;
  for (int i = 0; i < n; ++i) {
    const double e = xc[i] - slope * (i - centre) - rmean;
    rss += e * e;
  }
  const double rsd = std::sqrt(rss / (n - 1));
  if (std::isnan(rsd)) return NA_REAL;
  if (rsd <= std::sqrt(std::numeric_limits<double>::epsilon())) return 0.0;

  // Autocovariances at lags 0..order_max with divisor n, as
  // acf(type = "covariance") gives them to ar.yw().
  for (int lag = 0; lag <= order_max; ++lag) {
    w.r[lag] = dot(xc, xc + lag, n - lag) / n;
  }
  if (!(w.r[0] > 0.0)) return NA_REAL;

  // Levinson-Durbin recursion: the Yule-Walker coefficients and the
  // innovations variance at every order up to order_max. Only the sum of the
  // coefficients is needed at each order.
  w.v[0] = w.r[0];
  w.phisum[0] = 0.0;
  for (int m = 1; m <= order_max; ++m) {
    double acc = w.r[m];
    for (int j = 1; j < m; ++j) acc -= w.phi[j] * w.r[m - j];
    const double k = acc / w.v[m - 1];
    for (int j = 1; j < m; ++j) w.tmp[j] = w.phi[j] - k * w.phi[m - j];
    for (int j = 1; j < m; ++j) w.phi[j] = w.tmp[j];
    w.phi[m] = k;
    w.v[m] = w.v[m - 1] * (1.0 - k * k);
    double s = 0.0;
    for (int j = 1; j <= m; ++j) s += w.phi[j];
    w.phisum[m] = s;
  }

  // Order by AIC, as ar.yw() does: n * log(var.pred) + 2 * order + 2. A
  // non-finite or tied minimum means a numerically singular fit, where
  // ar.yw() either stops or depends on rounding; leave those to coda.
  int order = 0;
  int n_min = 0;
  double best = std::numeric_limits<double>::infinity();
  for (int m = 0; m <= order_max; ++m) {
    const double aic = n * std::log(w.v[m]) + 2.0 * m + 2.0;
    if (!std::isfinite(aic)) return NA_REAL;
    if (aic < best) {
      best = aic;
      order = m;
      n_min = 1;
    } else if (aic == best) {
      ++n_min;
    }
  }
  if (n_min != 1) return NA_REAL;

  const double var_pred = w.v[order] * n / (n - (order + 1.0));
  const double gain = 1.0 - w.phisum[order];
  const double spec = var_pred / (gain * gain);
  if (std::isnan(spec)) return NA_REAL;
  return spec == 0.0 ? 0.0 : n * var / spec;
}

}  // namespace

// draws is an array whose last two dimensions are iteration and chain and
// whose leading dimensions, flattened, index n_params parameters. It is read
// in place. Returns the mean, sample variance and effective sample size of
// every chain of every parameter, each as an n_params by nchain matrix. It
// draws no random numbers, so it is exported without Rcpp's RNG scope and
// leaves .Random.seed untouched.
// [[Rcpp::export(rng = false)]]
Rcpp::List chain_moments_ess_cpp(const Rcpp::NumericVector& draws,
                                 int n_params, int niter, int nchain) {
  if (n_params < 1 || niter < 1 || nchain < 1) {
    Rcpp::stop("n_params, niter and nchain must all be at least 1");
  }
  const R_xlen_t per_chain = static_cast<R_xlen_t>(n_params) * niter;
  if (draws.size() != per_chain * nchain) {
    Rcpp::stop("draws does not have n_params * niter * nchain elements");
  }

  Rcpp::NumericMatrix mean(n_params, nchain);
  Rcpp::NumericMatrix var(n_params, nchain);
  Rcpp::NumericMatrix ess(n_params, nchain);

  // The largest AR order ar.yw() considers by default.
  const int order_max = std::min(
      niter - 1, static_cast<int>(std::floor(10.0 * std::log10(niter))));
  ArWork work(std::max(order_max, 0));

  // One parameter's draws are n_params apart in memory. Reading them one
  // parameter at a time would miss the cache on every draw, so a group of
  // parameters is gathered into contiguous rows first, in the order the
  // array is laid out.
  const int group = 32;
  std::vector<double> buf(static_cast<std::size_t>(group) * niter);
  const double* x = draws.begin();

  for (int p0 = 0; p0 < n_params; p0 += group) {
    const int n_group = std::min(group, n_params - p0);
    for (int ch = 0; ch < nchain; ++ch) {
      const double* base = x + ch * per_chain + p0;
      for (int it = 0; it < niter; ++it) {
        const double* src = base + static_cast<R_xlen_t>(it) * n_params;
        for (int k = 0; k < n_group; ++k) {
          buf[static_cast<std::size_t>(k) * niter + it] = src[k];
        }
      }
      for (int k = 0; k < n_group; ++k) {
        double* xc = buf.data() + static_cast<std::size_t>(k) * niter;
        // The mean in two passes, as R computes it, so that a chain that
        // never moved has a mean equal to its value and a variance of
        // exactly zero.
        double s = 0.0;
        for (int i = 0; i < niter; ++i) s += xc[i];
        double m = s / niter;
        double t = 0.0;
        for (int i = 0; i < niter; ++i) t += xc[i] - m;
        m += t / niter;
        for (int i = 0; i < niter; ++i) xc[i] -= m;
        const double v = niter > 1
            ? dot(xc, xc, niter) / (niter - 1)
            : NA_REAL;
        mean(p0 + k, ch) = m;
        var(p0 + k, ch) = v;
        ess(p0 + k, ch) = chain_ess(xc, niter, v, order_max, work);
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(Rcpp::Named("mean") = mean,
                            Rcpp::Named("var") = var,
                            Rcpp::Named("ess") = ess);
}
