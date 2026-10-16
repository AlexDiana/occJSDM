get_param <- function(params, key, default = 0) {
  if (!is.null(params[[key]])) {
    return(params[[key]])
  } else {
    return(default)
  }
}

process_covariates <- function(data_info, covariates, group_by_col, n_obs,
                               remove_intercept = FALSE,
                               spline_vars = F) {

  if (length(covariates) > 0) {

    # NULL means the caller has already supplied one ordered row per unit.
    # Keep identifier columns available when requested as covariates.
    df <- data_info
    if (!is.null(group_by_col)) {
      df <- df %>%
        dplyr::group_by(!!rlang::sym(group_by_col)) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(covariates), ~ dplyr::first(.x)),
                         .groups = "drop")
    }
    df <- df %>%
      dplyr::select(dplyr::all_of(covariates)) %>%

      dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor))

    if (any(is.infinite(as.matrix(df)))) stop("Infinite values (Inf or -Inf) detected in covariates.")
    if (any(is.nan(as.matrix(df)))) stop("NaN values detected in covariates.")

    list_X <- create_covariates_matrix(
      df,
      spline_vars = spline_vars,
      remove_intercept = remove_intercept)
    out_matrix <- list_X$X
    list_matrix <- list_X$list_matrix
    X0 <- list_X$X0

  } else {

    if (!remove_intercept) {
      out_matrix <- (matrix(1, n_obs, 1))
    } else {
      out_matrix <- (matrix(0, n_obs, 0))
    }

    list_matrix <- NULL
    X0 <- NULL
  }

  list("df" = out_matrix,
       "list_matrix" = list_matrix,
       "X0" = X0)

}

createDataIdx <- function(n, M, P, K, twostage, primerId_p = NULL){
  # P is now a vector with one entry per sample (length N = sum(M)), giving
  # the number of primers actually used for that sample -- it no longer has
  # to be constant across samples. primerId_p (length N2 = sum(P)) gives the
  # *global* primer identity (an index into 1:length(primerNames)) for each
  # (sample, primer) block, in the same order samples/primers are visited
  # below; this is what lets p/q be looked up by actual primer identity
  # rather than by a primer's position within a given sample's block.

  N <- sum(M)
  sumM <- c(0, cumsum(M)[-n])
  idx_z_w <- rep(NA, N)

  if(twostage){
    N2 <- sum(P)
    N3 <- sum(K)
    sumP <- c(0, cumsum(P)[-N])
    sumK <- c(0, cumsum(K)[-N2])

    idx_z_p <- rep(NA, N2)
    idx_w_p <- rep(NA, N2)
    idx_z_k <- rep(NA, N3)
    idx_p_k <- rep(NA, N3)
    idx_w_k <- rep(NA, N3)
  } else {
    idx_z_p <- NULL
    idx_w_p <- NULL
    idx_z_k <- NULL
    idx_p_k <- NULL
    idx_w_k <- NULL
  }

  idx_z <- 1
  idx_w <- 1
  idx_p <- 1
  idx_k <- 1

  for (i in 1:n) {
    for (m in 1:M[i]) {
      idx_z_w[idx_w] <- i

      if(twostage){

        sample_idx <- sumM[i] + m

        for (pp in seq_len(P[sample_idx])) {
          idx_z_p[idx_p] <- i

          primer_id <- if (!is.null(primerId_p)) primerId_p[idx_p] else pp

          for (k in 1:K[sumP[sample_idx] + pp]) {

            idx_z_k[idx_k] <- idx_z
            idx_w_k[idx_k] <- idx_w
            idx_p_k[idx_k] <- primer_id

            idx_k <- idx_k + 1
          }
          idx_z_p[idx_p] <- idx_z
          idx_w_p[idx_p] <- idx_w

          idx_p <- idx_p + 1
        }

      }


      idx_w <- idx_w + 1
    }
    idx_z <- idx_z + 1
  }

  list("idx_z_w" = idx_z_w,
       "idx_z_p" = idx_z_p,
       "idx_w_p" = idx_w_p,
       "idx_p_k" = idx_p_k,
       "idx_w_k" = idx_w_k,
       "idx_z_k" = idx_z_k)

}

inferDataModel <- function(data){

  stages_type <- NULL

  data_info <- data$info
  OTU <- data$OTU

  site_present <- !is.null(data_info$Site)
  sample_present <- !is.null(data_info$Sample)
  primer_present <- !is.null(data_info$Primer)

  if(site_present){ # site column present

    if(all(!duplicated(data_info$Site))){ # no repeated observation at each site

      stages_type <- "no_stage"
      message("No repeated observation for any site")

    } else { # repeated observation at each site

      if(!sample_present){ # sample column missing

        message("Sample column missing, but detected repeated observations for some sites")
        stages_type <- "one_stage"

      } else { # sample column present

        if(all(!duplicated(data_info$Sample))){ # no repeated observation in each sample

          stages_type <- "one_stage"
          message("No repeated observation for any sample")

        } else {

          stages_type <- "two_stage"
          message("Detected repeated observation for some sample")

        }

      }

    }
  } else { # site column missing
    stages_type <- "no_stage"
    message("Site column missing")
  }

  # infer the data type
  {
    # if(stages_t)
    if(mode(OTU) != "numeric"){
      stop("Data in OTU need to be numeric")
    } else if(all(round(OTU) == OTU, na.rm = T)) { ## TODO: write that the data are integer
      if(all(unique(OTU) %in% c(0,1,NA), na.rm = T)){
        message("Only 0 and 1 detected")
        data_type <- "binary"
      } else {
        message("Count data detected")
        data_type <- "counts"
      }
    } else {
      message("Continuous observations detected")
      data_type <- "continuous"
    }
  }

  if(stages_type == "no_stage"){
    if(data_type == "binary"){
      model <- "binary"
    } else if(data_type == "continuous"){
      model <- "continuous"
    } else if(data_type == "counts"){
      model <- "counts"
    }
  } else if(stages_type == "one_stage"){
    model <- "occupancy"
  } else if(stages_type == "two_stage"){
    model <- "two_stage"
  } else {
    stop("No model recognised")
  }

  if(model == "binary") message(paste("occJSDM has inferred one-stage binary data"))
  if(model == "continuous") message(paste("occJSDM has inferred one-stage continuous data"))
  if(model == "counts") message(paste("occJSDM has inferred one-stage counts data"))
  if(model == "occupancy") message(paste("occJSDM has inferred occupancy data"))
  if(model == "two_stage") message(paste("occJSDM has inferred two stage (eDNA style) data"))


  message("Check your data if this is not what you were expecting!")

  if(model == "counts") stop("Counts model not supported yet :(")

  model

}

create_waic_quantities <- function(n_obs){

  M2 <- rep(0, n_obs)
  mean_log <- rep(0, n_obs)
  mean_lik <- rep(0, n_obs)

  list("M2" = M2,
       "mean_log" = mean_log,
       "mean_lik" = mean_lik)
}

#' runOccJSDM
#'
#' Fit an occupancy joint species distribution model (occJSDM), optionally
#' accounting for environmental and detection covariates, species traits,
#' spatial autocorrelation, and for eDNA-style data, a two-stage observation
#' process (false-negative and false-positive detection errors in the field
#' and in the lab).
#'
#' @details
#' The model actually fit is inferred automatically from the shape of
#' \code{data} (via the internal \code{inferDataModel()}): whether
#' \code{data$info} has repeated observations per \code{Site} (and, within a
#' \code{Site}, per \code{Sample}) determines whether the data are treated
#' as single-visit JSDM data, classical (one-stage) occupancy data, or
#' two-stage (eDNA-style) occupancy/detection data. A message is printed
#' describing which model was inferred -- check that it matches expectations
#' for your data.
#'
#' @param data A list with two elements:
#' \describe{
#'   \item{info}{A data.frame with one row per observation (i.e. per PCR
#'   replicate for two-stage eDNA data), containing \code{Site}, \code{Sample},
#'   and \code{Primer} id columns (as applicable to the inferred model) plus
#'   the occupancy/collection/spatial covariate columns named in
#'   \code{occCovariates}/\code{collCovariates}/\code{spatCovariates}. See the
#'   vignette for the full expected shape.}
#'   \item{OTU}{A matrix of dimension (N x S), where N is
#'   \code{nrow(data$info)} and S is the number of species, containing the
#'   number of reads of each species in each observation.}
#'   \item{traits}{(Optional) species (rows) by trait (columns) matrix or
#'   data.frame of species traits (numeric and/or categorical), matched to
#'   \code{colnames(data$OTU)} by row name.}
#' }
#' @param listParams (Optional) list of model-size hyperparameters:
#' \describe{
#'   \item{n_factors}{Number of latent JSDM factors used for the residual
#'   species covariance. Capped to the number of species if larger.}
#'   \item{n_lattrait}{Number of latent trait factors (\code{gt}). Default
#'   \code{floor(sqrt(min(S, ncov_psi)))}, where \code{S} is the number of
#'   species and \code{ncov_psi} is the number of columns in the processed
#'   occupancy-covariate design matrix. Set to \code{0} to disable latent
#'   trait factors.}
#'   \item{n_supportpoints}{Number of spatial support points used to
#'   approximate the Gaussian process over site coordinates when
#'   \code{spatCovariates} is non-empty. Defaults to
#'   20 percent of the number of unique coordinate locations, rounded down.
#'   Capped to the number of unique locations. Requesting that number uses
#'   every location as a support point. Short-range spatial patterns may
#'   require more support points; check that increasing their number does
#'   not materially change the estimated field or occupancy probabilities.}
#' }
#' @param threshold Threshold used to truncate the reads to binary detections
#' for occupancy/two-stage models. Reads greater than or equal to the
#' threshold are considered a detection (default \code{1}). Must be
#' a single positive numeric value. Values between \code{0} and \code{1}
#' are accepted, but observations are still converted to binary detections.
#' For integer read counts, positive thresholds up to \code{1} give the
#' same detections as the default. Continuous detection modeling
#' (\code{threshold = 0}) is not supported.
#' @param occCovariates (Optional) vector of the names of the occupancy
#' covariates (site-level). Names should match column names in
#' \code{data$info}.
#' @param collCovariates (Optional) vector of the names of the collection
#' covariates (sample-level, relevant for occupancy/two-stage models only).
#' Names should match column names in \code{data$info}.
#' @param spatCovariates (Optional) vector of the names of the site
#' coordinate columns used to fit a spatially autocorrelated random field
#' over occupancy. Names should match column names in \code{data$info}. If
#' omitted, no spatial field is fit.
#' @param MCMCparams (Optional) list of MCMC settings: \code{nchain} (number
#' of chains), \code{nburn} (number of burn-in iterations to discard),
#' \code{niter} (number of post-burn-in iterations to keep), \code{nthin}
#' (thinning interval, default \code{1}) and \code{ncores} (the total core
#' budget, default \code{1}; it must be a single whole number of at least 1,
#' and a value above the number of cores detected on the machine gives a
#' warning and the fit proceeds). With \code{ncores = 1} the chains run one
#' after another in this process. With \code{ncores > 1} they run on
#' \code{min(nchain, ncores)} worker processes, each with
#' \code{floor(ncores / workers)} RcppParallel threads; for example four chains
#' with \code{ncores = 4} run on four single-threaded workers and with
#' \code{ncores = 8} on four workers of two threads. Any remainder of the
#' budget is idle: \code{ncores = 3} with two chains gives two single-threaded
#' workers. A single chain never
#' leaves this process: with \code{nchain = 1} and \code{ncores > 1} the fit
#' runs here with its RcppParallel thread count set to \code{ncores} for the
#' duration and restored afterwards, and no progress file is written. Defaults
#' to \code{list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1, ncores = 1)}.
#' @param summarisedLatentPresences (Optional) logical, default \code{TRUE}.
#' If \code{TRUE}, the returned latent-presence/probability arrays
#' (\code{z_output}, \code{w_output}, \code{psi_output}, \code{theta_output})
#' are collapsed to posterior means, which keeps the fitted object small. If
#' \code{FALSE}, the full posterior samples are kept instead, which allows
#' computing credible intervals for conditional/predictive occupancy
#' probabilities but produces a much larger object.
#' @param listPriors (Optional) list of prior hyperparameters. Currently
#' supported: \code{a_p}/\code{b_p} (Beta prior shape parameters for the true
#' positive lab stage rate,
#' default \code{5}/\code{1}), \code{a_q}/\code{b_q} (Beta prior shape parameters
#' for the false positive
#' lab stage rate, default \code{1}/\code{20})  and
#' \code{a_theta0}/\code{b_theta0} (Beta prior shape parameters for the
#' false positive field stage rate, default \code{1}/\code{20}), and
#' \code{b_betatheta_slope_var} (prior variance on the collection-covariate
#' slopes, default \code{2}). The last of these is provisional: it exists so
#' the slope prior can be varied without editing the source, and the default
#' reproduces the previous hard-coded behaviour exactly. Its value is still
#' under review, so treat a non-default setting as a diagnostic rather than a
#' recommended configuration.
#'
#' For continuous responses, \code{tau_prior} selects the noise prior:
#' \code{"inverse_gamma"} (current default) places an inverse-gamma prior on
#' each species' noise variance, with \code{a_tau} (shape, default \code{5})
#' and \code{b_tau} (rate for the reciprocal variance, default \code{5}).
#' \code{"half_cauchy"} places a half-Cauchy prior on the noise standard
#' deviation, with \code{tau_scale} (default \code{1}) in response units.
#' The half-Cauchy option allows noise close to zero; the inverse-gamma
#' default strongly discourages it. Each scale, shape and rate must be a
#' finite positive number. These settings do not change the binary,
#' occupancy or two-stage detection models. The chosen continuous noise
#' prior is saved in \code{infos$noise_prior}.
#'
#' For spatial fits, \code{sigma_bs_prior = "half_cauchy"} enables an
#' experimental half-Cauchy prior on the shared residual spatial-coefficient
#' standard deviation, with \code{sigma_bs_scale} (default \code{1}) on the
#' linear-predictor scale. This amplitude is not necessarily the realised
#' spatial field's standard deviation, especially with reduced support.
#' The default \code{sigma_bs_prior = "inverse_gamma"} retains the existing
#' inverse-gamma prior on its variance (shape \code{10}, scale \code{1}).
#' These settings require spatial covariates; \code{sigma_bs_scale} requires
#' the half-Cauchy choice and must be a finite positive number. The applied
#' choice is saved in \code{infos$spatial_sd_prior}. Half-Cauchy chains start
#' at the supplied scale (the prior median); default chains retain their
#' existing starting value. This is separate from
#' continuous observation noise and leaves all other occupancy and detection
#' priors unchanged.
#' In full-rank binary spatial fits, the unbounded half-Cauchy can leave the
#' spatial-amplitude posterior mean infinite and spatial-field means
#' non-integrable. Use posterior medians and quantiles for spatial summaries
#' with this experimental option, and check rank-based chain convergence.
#' Bounded occupancy-probability means remain meaningful. This option does
#' not establish improved spatial recovery.
#'
#' \code{sigma_b0} (experimental) sets how far each species' baseline
#' occupancy \code{B0} may range. \code{B0} is the species' occupancy at a
#' site with average covariate values, on the logit scale, and has a Normal(0,
#' \code{sigma_b0}^2) prior. The default \code{1} is the previous fixed prior
#' and leaves default fits unchanged. It puts about 95\% of the prior on
#' baseline occupancies between 12\% and 88\%, so very rare or very common
#' species can be pulled towards the middle; \code{sigma_b0 = 2} widens that
#' range to about 2\% to 98\%. The value must be a finite positive number. It
#' applies to spatial and non-spatial fits and does not affect the collection
#' (detection) intercept. The applied prior is saved in
#' \code{infos$intercept_prior}. For continuous data \code{B0} is on the
#' response scale, so the occupancy interpretation above does not apply.
#'
#' A simulation study compared \code{sigma_b0} = 1, 2, 3 and 5
#' (\code{dev/simstudy/occupancy-intercept-prior/} in the source repository).
#' In binary fits, wider values reduced the overestimation of low occupancy
#' probabilities, clearly in spatial fits and only slightly in non-spatial
#' fits. In two-stage (eDNA) fits the data only weakly separate how often a
#' species is present from how often a sample catches it when present, and the
#' default prior probably helps by keeping \code{B0} from drifting along that
#' trade-off: wider values made the MCMC mix worse, markedly at 3 and 5, and 3
#' and 5 also increased the bias of probabilities between 0.2 and 0.8.
#' Occupancy fits share that collection stage and were not tested, and neither
#' were continuous fits. The default therefore stays at \code{1}. Keep it for
#' occupancy and two-stage data; if you try a larger value for binary data,
#' check chain convergence.
#'
#' \code{b0_prior = "community"} (experimental, off by default) estimates the
#' centre of that prior instead of fixing it at 0: each species' \code{B0} has
#' a Normal(\code{mu_b0}, \code{sigma_b0}^2) prior, and the community mean
#' \code{mu_b0} has a Normal(0, \code{mu_b0_sd}^2) prior with \code{mu_b0_sd}
#' defaulting to \code{2}, which puts about 95\% of the prior on community
#' mean baseline occupancies between 2\% and 98\%. \code{mu_b0} is drawn from
#' its conjugate conditional distribution each iteration, so the species are
#' pulled towards their own community's average rather than towards 50\%,
#' which matters when most species are rare. Its draws are returned in
#' \code{results_output$jsdm_output$mu_b0_output}, a \code{niter} by
#' \code{nchain} array, and the applied prior is saved in
#' \code{infos$intercept_prior}. The default, \code{b0_prior = "fixed"},
#' leaves fits unchanged; \code{mu_b0_sd} requires \code{"community"} and must
#' be a finite positive number. With few species the community mean is weakly
#' determined, and in occupancy and two-stage fits it can trade off against
#' the collection intercepts, so check the convergence of \code{mu_b0} and
#' \code{B0}. This option is under review and is not a recommended setting.
#'
#' \code{theta_intercept_prior = "community"} (experimental, off by default)
#' does the same for the collection intercepts, \code{beta_theta[1, ]}, of
#' occupancy and two-stage fits. Each has a fixed Normal(0, 1) prior by
#' default; with this option it has a Normal(\code{mu_theta}, 1) prior, and
#' the community mean \code{mu_theta} has a Normal(0, \code{mu_theta_sd}^2)
#' prior with \code{mu_theta_sd} defaulting to \code{2}. \code{mu_theta} is
#' drawn from its conjugate conditional distribution each iteration, after the
#' collection coefficients, and a species with no occupied site in an
#' iteration has its collection coefficients drawn from their prior. Its
#' draws are returned in \code{results_output$mu_theta_output}, a
#' \code{niter} by \code{nchain} array, and the applied prior is saved in
#' \code{infos$collection_intercept_prior}. It can be used alone or with
#' \code{b0_prior = "community"}, so that neither the occupancy nor the
#' collection side has a fixed prior centre. The default,
#' \code{theta_intercept_prior = "fixed"}, leaves fits unchanged;
#' \code{mu_theta_sd} requires \code{"community"} and must be a finite
#' positive number, and binary and continuous fits, which have no collection
#' stage, reject the option. In simulations of rare-dominated communities
#' with 2 and 3 field samples per site, the two options together placed the
#' community's occupancy close to the truth, but the rarest species were
#' still estimated too common and their collection intercepts too low; at 3
#' field samples per site \code{b0_prior = "community"} alone did better.
#' \code{mu_b0} mixed slowly in both: use chains of at least the default
#' length and check the convergence of \code{mu_theta}, \code{mu_b0} and the
#' intercepts. This option is under review and is not a recommended setting.
#'
#' @section Field-contamination uncertainty:
#' \code{theta0} is the probability that DNA enters a field sample when
#' the species is absent from the site. Its default is Beta(1,20), set by
#' \code{a_theta0 = 1} and \code{b_theta0 = 20} in \code{listPriors}.
#' In targeted low-contamination, non-spatial simulations with ten
#' communities per design, its 95\% intervals contained the truth 97-98\%
#' of the time. General interval calibration remains unestablished.
#'
#' Compared with historical Beta(1,30), the current prior produced intervals
#' about 10-19\% wider and reduced downward mean error. Average absolute
#' error changed little, with no clear advantage for either prior, and some
#' validation fits retained numerical diagnostic flags. The default remains
#' Beta(1,20). Check chain diagnostics and prior sensitivity; tightening the
#' prior solely to reach 95\% coverage in selected simulations is not
#' supported by this study. This concerns field contamination, separately
#' from the laboratory false-positive rate \code{q}. See the
#' \href{https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/theta0-prior-mean/REPORT.md}{prior-sensitivity report}
#' and the repository README's Known limitations.
#'
#' @return A list with:
#' \describe{
#'   \item{results_output}{Posterior samples/summaries, including
#'   \code{jsdm_output} (JSDM coefficient posteriors), \code{beta_theta_output},
#'   \code{p_output}, \code{q_output}, \code{theta0_output} (detection-process
#'   posteriors), the legacy complete-data score \code{WAIC} (use
#'   \code{computeSiteWAIC()} for observed-data model comparison), and
#'   (per \code{summarisedLatentPresences})
#'   either posterior means or full samples of \code{z_output} (latent
#'   occupancy), \code{w_output} (latent collection), \code{psi_output}
#'   (occupancy probabilities), and \code{theta_output} (collection
#'   probabilities).}
#'   \item{infos}{Metadata about the fitted model (e.g. \code{S}, \code{n},
#'   \code{M}, \code{K}, \code{speciesNames}, \code{primerNames},
#'   \code{siteNames}, \code{ncov_psi}, \code{ncov_theta}, the \code{OTU}
#'   matrix used for fitting, and standardisation info for \code{X_psi}/
#'   \code{Xs}), used internally by the plotting/output helpers (e.g.
#'   \code{\link{plotOccupancyRates}}, \code{\link{plotDetectionRates}}). It
#'   also records \code{MCMCparams} (the resolved \code{nchain}, \code{nburn},
#'   \code{niter} and \code{nthin}), \code{chain_seeds} (a 2 by \code{nchain}
#'   matrix of the R and C++ seeds of each chain) and \code{execution} (how
#'   the fit was run: \code{ncores}, \code{workers}, \code{threads_per_worker},
#'   \code{thread_setting} (\code{"budget"} or \code{"ambient"}),
#'   \code{blas_threads_requested} (the BLAS thread settings requested of the
#'   workers), \code{worker_pids} (the worker process ids),
#'   \code{blas_declared} (whether the environment variable
#'   \code{OCCJSDM_REFERENCE_BLAS} was set), \code{blas} and \code{lapack}
#'   (the BLAS and LAPACK libraries R reports), \code{elapsed} (the elapsed
#'   seconds per chain) and \code{progress_file}). Only \code{execution}
#'   differs between two fits of the same seed at different \code{ncores} when
#'   the BLAS is single-threaded, as the Details explain.}
#'   \item{Tr}{The species trait matrix actually used to fit the model.}
#'   \item{X_theta}{The collection covariate design matrix (including
#'   intercept) used to fit the model.}
#'   \item{X_s}{The site coordinates used for the spatial field.}
#'   \item{X_psi}{The (standardised) occupancy covariate design matrix used
#'   to fit the model.}
#' }
#'
#' @details
#' Every chain is seeded on its own: \code{runOccJSDM()} draws two seeds per
#' chain from R's random number generator, one for the R-level draws and one
#' for the package's C++ sampler, so \code{set.seed()} before fitting
#' reproduces a fit on the same platform and build, and the fit is the same
#' whatever \code{ncores} says: a chain's draws depend on its seeds alone, not
#' on which chains ran before it or in which process. Consecutive fits without
#' resetting the seed consume new draws. Because of the per-chain seeding, a
#' seeded fit's numbers differ from those of the 0.1.1 beta, in which chains
#' continued one stream.
#'
#' Random sampling runs serially on R's main thread; deterministic probability
#' calculations may use several RcppParallel threads, and the thread count
#' does not change the sampling stream. The calling session's RcppParallel
#' thread setting is preserved. Identity across \code{ncores} also requires
#' that the BLAS does not change its arithmetic with the thread count. R's
#' reference BLAS, the one CRAN's macOS and Windows installers ship, is
#' single-threaded. Worker processes are asked to run any other BLAS
#' single-threaded through the usual thread variables, a request that has not
#' been verified against a threaded BLAS. The package cannot verify which BLAS
#' it is running on, and \code{extSoftVersion()} only indicates the library: a
#' name such as \code{libRblas} is an indication, not proof. With an unpinned
#' OpenBLAS, MKL or Accelerate an in-process fit can differ from worker fits.
#'
#' A fit advances the session's random number stream by its data preparation
#' and by the chain seeds it draws, and no further: running the chains, in
#' process or in workers, does not move it, so the session ends the fit in the
#' same state at every \code{ncores} and consecutive fits receive different
#' seeds. One caveat: under the Box-Muller normal kind R keeps a cached deviate
#' outside \code{.Random.seed}, which a fit discards, so the next normal draw
#' after a fit can differ from what an unfitted session would have given.
#'
#' A parallel fit prints one line naming its progress file, which the workers'
#' iteration messages go to; use the \code{tail -f} command it shows to watch
#' them. Warnings raised inside a chain in a parallel fit also go to that file,
#' where each worker prints them when it stops at the end of the fit, and they
#' are not raised again in this session.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' fitmodel <- runOccJSDM(data,
#' listParams = list(n_factors = 2),
#' occCovariates = c("X_psi.1", "X_psi.2"),
#' collCovariates = c("X_theta"))
#' }
#'
#' @export
#' @import dplyr
#'
runOccJSDM <- function(data,
                       listParams = list(),
                       threshold = 1,
                       occCovariates = c(),
                       collCovariates = c(),
                       spatCovariates = c(),
                       MCMCparams = list(nchain = 2,
                                         nburn = 5000,
                                         niter = 5000,
                                         nthin = 1),
                       summarisedLatentPresences = T,
                       listPriors = list()){

  {
    # Debugging scaffold: a set of arguments to paste in when stepping through
    # this function by hand. Everything here must stay commented out. Live
    # exactly once (8f9f315) breaks the function outright -- every call fails
    # with "object 'occ_data_effort' not found" -- and this is now the fourth
    # time it has happened (8f9f315, 46d8804, and again between e5d0105 and
    # 11981a1). 11981a1 re-commented it; 522b89e uncommented it again in the
    # same "some more parallelisation" commit, seemingly as a side effect of
    # editing nearby lines rather than a deliberate change. If you need this
    # for stepping through by hand, paste it into the console, not here.
    #
    # data = occ_data_effort
    # listParams = list(n_factors = 3)
    # threshold = 1
    # occCovariates = c("season_year","z_prop_closed")
    # collCovariates = c("predator_season_year", "z_log_effort_m_total")
    # spatCovariates = c()
    # MCMCparams = list(
    # nchain = 2,
    # nburn = 5000,
    # niter = 5000,
    # nthin = 1 )
    # listPriors = list(
    # a_theta0 = 1,
    # b_theta0 = 100)# listParams = list(n_factors = 3, splineVars = F)
    # summarisedLatentPresences = T
    # threshold = 1
    # occCovariates = c( "X_psi.EnvCov.1", "X_psi.EnvCov.2", "X_psi.EnvCov.3",
    #                    "X_psi.EnvCov.4", "X_psi.EnvCov.5" ,"X_psi.EnvCov.6")
    # collCovariates = c("X_theta.1","X_theta.2")
    # spatCovariates = c("Xs.1","Xs.2")
    # MCMCparams = list(
    #   nchain = 2,
    #   nburn = 200,
    #   niter = 200,
    #   nthin = 1 )
    # listPriors = list(
    #   a_theta0 = 1,
    #   b_theta0 = 100)
    # listPriors = list()
  }

  # data structure infer
  {
    if(is.null(data$info) || is.null(data$OTU)){
      stop("data_info or OTU missing")
    } else if (is.vector(data$OTU)) {
      stop("data$OTU cannot be a vector. Please provide it as a matrix or data frame.")
    }else {
      data_info <- as.data.frame(data$info)
      OTU <- data$OTU
    }

    if(nrow(data_info) != nrow(OTU)){
      stop("OTU and data_info cannot have different number of rows")
    }

    model <- inferDataModel(data)

    if(model %in% c("binary","occupancy","two_stage")){
      jsdmModel <- "binary"
    } else if(model == "continuous") {
      jsdmModel <- "continuous"
    }  else if(model == "counts") {
      jsdmModel <- "counts"
    }

    if(is.null(data_info$Site)){
      data_info$Site <- 1:nrow(data_info)
    }

    if(model == "occupancy" & is.null(data_info$Sample)){
      data_info$Sample <- 1:nrow(data_info)
    }

    data_info <- data_info %>%
      mutate(rownum = row_number())

    if(model == "occupancy"){
      data_info <- data_info %>%
        dplyr::arrange(Site, Sample)
    } else if(model == "two_stage"){
      data_info <- data_info %>%
        dplyr::arrange(Site, Sample, Primer)

    }

    OTU <- OTU[data_info$rownum, , drop = FALSE]

    data_info <- data_info %>% select(-rownum)
  }

  # read OTU data
  {
    y <- OTU
    if(is.null(dim(y))){
      S <- 1
    } else {
      S <- ncol(y)
    }

    if(model %in% c("occupancy","two_stage")){

      # truncate data

      if (!is.numeric(threshold) || length(threshold) != 1) {
        stop("'threshold' must be a single numeric value.")
      }
      if (threshold <= 0) { # You already note threshold > 0 in your code, enforce it cleanly!
        stop("'threshold' must be strictly greater than 0.")
      }

      readsAboveThreshold <- y >= threshold
      readsBelowThreshold <- y < threshold
      y[readsAboveThreshold] <- 1
      y[readsBelowThreshold] <- 0

    }

    if(model %in% c("binary","continuous","counts")){
      z <- y
    } else if(model == "occupancy"){
      w <- y
    }

    # check for nas
    {
      y_NA <- is.na(y)
      mode(y_NA) <- "integer"

      if(sum(y_NA) > 0 & model != "two_stage"){
        stop("NAs are allowed only in the two-stage model", .call = F)
      }
    }

    speciesNames <- colnames(data$OTU)
    if(is.null(speciesNames)){
      speciesNames <- 1:S
    }

    if (any(duplicated(speciesNames))) stop("Duplicated species names found in OTU matrix.")

    n_obs <- nrow(y)


  }

  # data checks
  {

    if(any(is.na(data_info$Site)) |
       any(is.na(data_info$Sample)) |
       any(is.na(data_info$Primer))){
      stop("NA in Site, Sample or Primer columns")
    }

    if(is.null(occCovariates)) occCovariates <- c()
    if(is.null(collCovariates)) collCovariates <- c()
    if(is.null(spatCovariates)) spatCovariates <- c()

    if(!(length(spatCovariates) %in% c(0,2))) stop("Two columns needed for the spatial covariates")

    if(!all(c(occCovariates, collCovariates, spatCovariates) %in% colnames(data$info))){
      stop("Covariate names provided not in data$info")
    }

    if(!all(sapply(data_info[,spatCovariates], is.numeric))) {
      stop("Non numeric columns in spatial covariates")
    }

  }

  # clean the data
  {

    # samples per site
    {
      if(model %in% c("occupancy","two_stage")){
        # Canonical sample rows use the same typed (Site, Sample) ordering
        # as the primer/PCR tables below and the latent-state indices.
        data_samples <- data_info %>%
          dplyr::group_by(Site, Sample) %>%
          dplyr::slice(1) %>%
          dplyr::ungroup()

        M_df <- data_samples %>%
          dplyr::group_by(Site) %>%
          dplyr::summarise(M = n(),
                           .groups = "keep")

        M <- M_df$M
        names(M) <- M_df$Site
        siteNames <- M_df$Site

        n <- length(M)
        N <- sum(M)

        sumM <- c(0, cumsum(M)[-n])

        # Retain the display label, but never use pasted keys for grouping:
        # numeric IDs sort lexically and distinct pairs can share a label.
        data_info$SiteSample <- paste(data_info$Site,data_info$Sample, sep = "-")

      } else {
        n <- nrow(data_info)
        siteNames <- 1:n
        M <- NULL
      }

    }

    # marker per samples
    {
      if(model == "two_stage"){
        # number of markers

        P_df <- data_info %>%
          dplyr::group_by(Site, Sample, Primer) %>%
          dplyr::slice(1) %>%
          dplyr::group_by(Site, Sample) %>%
          dplyr::summarise(P = n(),
                           .groups = "keep")

        # P is now allowed to differ across samples: each sample can use a
        # different subset (and number) of primers.
        P <- P_df$P
        names(P) <- P_df$Sample

        # primerNames is the *global* set of distinct primers across the
        # whole dataset; p/q are indexed by primer identity (row into this
        # vector), not by a primer's position within a given sample.
        primerNames <- sort(unique(data_info$Primer))
        maxP <- length(primerNames)

        N2 <- sum(P)
        sumP <- c(0, cumsum(P)[-N])

      } else {
        P <- NULL
        primerNames <- NULL
      }
    }

    # pcr per marker
    {
      if(model == "two_stage"){
        data_K <- data_info %>%
          dplyr::group_by(Site, Sample, Primer) %>%
          dplyr::summarise(K = n(),
                           dplyr::across(contains("Species"),function(x){sum(x > 0)}),
                           .groups = "keep"
          ) %>%
          dplyr::ungroup()

        K <- data_K$K
        sumK <- c(0, cumsum(K)[-N2])

        N3 <- sum(K)

        # global primer identity (index into primerNames) for each
        # (Site, Sample, Primer) block -- data_K is at the same
        # (Site, Sample, Primer) granularity and in the same order (both
        # derived from the already Site/Sample/Primer-sorted data_info), so
        # this aligns row-for-row with K/sumP.
        primerId_p <- match(data_K$Primer, primerNames)

      } else {
        K <- NULL
        primerId_p <- NULL
      }

    }

    if(model %in% c("occupancy","two_stage")){
      list_idx <- createDataIdx(n, M, P, K, model == "two_stage", primerId_p = primerId_p)
      idx_z_w <- list_idx$idx_z_w
      idx_z_k <- list_idx$idx_z_k
      idx_w_p <- list_idx$idx_w_p
      idx_z_p <- list_idx$idx_z_p
      idx_w_k <- list_idx$idx_w_k
      idx_p_k <- list_idx$idx_p_k
    } else {
      list_idx <- NULL
    }

  }

  # create covariates matrix
  {

    splineVars <- get_param(listParams, "splineVars", F)

    # For occupancy covariates (group by Site, includes intercept)
    {
      list_X_psi <- process_covariates(data_info, occCovariates, "Site", n,
                                  remove_intercept = TRUE,
                                  spline_vars = splineVars)
      X_psi <- list_X_psi$df
      list_Xpsi_mat <- list_X_psi$list_matrix
      X0_psi <- list_X_psi$X0
    }

    # For the spatial field
    {
      list_Xs <- process_covariates(data_info, spatCovariates, "Site", n,
                               remove_intercept = TRUE,
                               spline_vars = F)
      Xs <- list_Xs$df
      list_Xs_mat <- list_Xs$list_matrix
    }

    # For collection covariates (canonical sample rows, includes intercept)
    {
      if(model %in% c("occupancy","two_stage")){
        list_X_theta <- process_covariates(data_samples, collCovariates, NULL,
                                           N, remove_intercept = FALSE,
                                           spline_vars = F)
        X_theta <- list_X_theta$df
        list_X_theta_mat <- list_X_theta$list_matrix
      } else {
        X_theta <- NULL
        list_X_theta_mat <- NULL
      }
    }

    # For the traits
    {
      if (!is.null(data$traits) && any(duplicated(rownames(data$traits)))) {
        stop("Duplicated species names found in traits matrix rownames.")
      }

      if(!is.null(data$traits)){
        speciesNamesInTraitsMatrix <- rownames(data$traits)
        if(!all(speciesNames %in% speciesNamesInTraitsMatrix)){
          stop("Species names in OTU not present in traits matrix")
        }

        idx_speciesNames <- match(speciesNames, speciesNamesInTraitsMatrix)
        Tr <- data$traits
        Tr <- Tr[idx_speciesNames, , drop = FALSE]

        list_Tr <- create_covariates_matrix(
          Tr,
          spline_vars = FALSE,
          remove_intercept = TRUE
        )
        Tr <- list_Tr$X
        list_Tr_mat <- list_Tr$list_matrix
        rownames(Tr) <- speciesNames
        traitsNames <- colnames(Tr)
      } else {
        Tr <- matrix(NA, S, 0)
        list_Tr_mat <- NULL
      }

    }

    ncov_psi <- ncol(X_psi)
    ncov_theta <- ncol(X_theta)
    g <- ncol(Tr)

  }

  # set unknonwn parameters
  {
    d <- get_param(listParams, "n_factors")

    if(ncol(OTU) > 1){
      if(d > ncol(OTU)){
        message("More species than factors. The number of factors will be capped to the
            number of species")
        d <- ncol(OTU)
      }
    } else {
      if(d > 0){
        message("Cannot have a factor model with just one species")
        d <- 0
      }
    }


    gt_default <- floor(sqrt(min(S, ncov_psi)))
    gt <- get_param(listParams, "n_lattrait", gt_default)

    if(ncol(Xs) > 0){
      ns <- nrow(unique(Xs))
      ps <- get_param(listParams, "n_supportpoints", getDefaultSupportPoints(ns))

    } else {
      ps <- 0
    }


  }

  # priors
  {
    a_theta0 <- ifelse(is.null(listPriors$a_theta0), 1, listPriors$a_theta0)
    b_theta0 <- ifelse(is.null(listPriors$b_theta0), 20, listPriors$b_theta0)
    a_p <- ifelse(is.null(listPriors$a_p), 5, listPriors$a_p)
    b_p <- ifelse(is.null(listPriors$b_p), 1, listPriors$b_p)
    a_q <- ifelse(is.null(listPriors$a_q), 1, listPriors$a_q)
    b_q <- ifelse(is.null(listPriors$b_q), 20, listPriors$b_q)

    if(model %in% c("occupancy","two_stage")){
      b_betatheta <- rep(0, ncov_theta)

      b_betatheta_slope_var <- ifelse(is.null(listPriors$b_betatheta_slope_var),
                                      2, listPriors$b_betatheta_slope_var)
      B_betatheta <- diag(b_betatheta_slope_var, nrow = ncov_theta)

      prior_beta_theta <- 0
      prior_beta_theta_sd <- 1
      b_betatheta[1] <- prior_beta_theta
      B_betatheta[1,1] <- prior_beta_theta_sd^2
    }
    # NULL unless the experimental theta_intercept_prior = "community" is set.
    collection_intercept_prior <- read_collection_intercept_prior(listPriors, model)

  }

  # chain parameters
  {
    nchain <- get_param(MCMCparams, "nchain", 2)
    nburn <- get_param(MCMCparams, "nburn", 1000)
    niter <- get_param(MCMCparams, "niter", 1000)
    nthin <- get_param(MCMCparams, "nthin", 1)
    ncores <- resolve_ncores(MCMCparams)
  }

  # precompute spatial quantities
  {
    # Spatial covariates matrix
    list_Xs <- computeSpatialSummaries(Xs, ps)
    Xs_centers <- list_Xs$Xs_centers
    Xs_index <- list_Xs$Xs_index
    X_s_centers <- list_Xs$X_s_centers
    X_tilde <- list_Xs$X_tilde
    X_s <- list_Xs$X_s
    ps <- list_Xs$ps

    length_grid_ls <- 10
    l_s_grid <- seq(0.01, 0.3, length.out = length_grid_ls)
    list_SoRSummaries <- precomputeSORmatrices(l_s_grid, list_Xs)
  }

  # precompute jsdm params
  {
    list_data <- list(
      "X" = X_psi,
      "Tr" = Tr)

    a_sigmab <- 10; b_sigmab <- 1
    a_sigmabs <- 10; b_sigmabs <- 1
    spatial_sd_prior <- read_spatial_sd_prior(listPriors,ps>0)
    intercept_prior <- read_intercept_prior(listPriors)
    a_sigmah <- 10; b_sigmah <- 1
    a_tau <- 5; b_tau <- 5
    noise_prior <- if (model == "continuous") read_noise_prior(listPriors) else NULL
    if (identical(noise_prior$type,"inverse_gamma")) {
      a_tau <- noise_prior$shape
      b_tau <- noise_prior$rate
    }
    a_l_s <- 1; b_l_s <- 1

    list_priors <- list(
      "a_sigmab" = a_sigmab,
      "b_sigmab" = b_sigmab,
      "a_sigmabs" = a_sigmabs,
      "b_sigmabs" = b_sigmabs,
      "spatial_sd_prior" = spatial_sd_prior,
      "intercept_prior" = intercept_prior,
      "a_sigmah" = a_sigmah,
      "b_sigmah" = b_sigmah,
      "a_tau" = a_tau,
      "b_tau" = b_tau,
      "noise_prior" = noise_prior,
      "a_l_s" = a_l_s,
      "b_l_s" = b_l_s
    )
  }

  # chain output
  {
    if(model %in% c("occupancy","two_stage")){
      beta_theta_output <- array(NA, dim = c(ncov_theta, S, niter, nchain))
      theta0_output <- array(NA, dim = c(S, niter, nchain))
      if(!summarisedLatentPresences){
        z_output <- array(NA, dim = c(n, S, niter, nchain))
        psi_output <- array(NA, dim = c(n, S, niter, nchain))
      }

    } else {
      beta_theta_output <- NULL
      theta0_output <- NULL
      if(summarisedLatentPresences){
        z_output_mean <- NULL
        psi_output_mean <- NULL
      } else {
        z_output <- NULL
        psi_output <- NULL
      }
      theta_output_mean <- NULL

    }

    if(model == "two_stage"){
      p_output <- array(NA, dim = c(maxP, S, niter, nchain))
      q_output <- array(NA, dim = c(maxP, S, niter, nchain))
    } else {
      p_output <- NULL
      q_output <- NULL
      w_output_mean <- NULL
    }

    # jsdm params

    {
      B0_output <- array(NA, dim = c(S, niter, nchain))
      B_output <- array(NA, dim = c(ncov_psi, S, niter, nchain))
      G_output <- array(NA, dim = c(g, ncov_psi, niter, nchain))
      A_output <- array(NA, dim = c(S, gt, niter, nchain))
      C_output <- array(NA, dim = c(gt, ncov_psi, niter, nchain))
      Bs_output <- array(NA, dim = c(ps, S, niter, nchain))
      Gs_output <- array(NA, dim = c(g, ps, niter, nchain))
      As_output <- array(NA, dim = c(S, gt, niter, nchain))
      Cs_output <- array(NA, dim = c(gt, ps, niter, nchain))
      U_output  <- array(NA, dim = c(n, d, niter, nchain))
      L_output  <- array(NA, dim = c(d, S, niter, nchain))
      tau_output <- array(NA, dim = c(S, niter, nchain))
      sigmab_output <- array(NA, dim = c(niter, nchain))
      sigmabs_output <- array(NA, dim = c(niter, nchain))
      sigmah_output <- array(NA, dim = c(niter, nchain))
      idx_ls_output <- array(NA, dim = c(niter, nchain))
      varPart_output <- array(NA, dim = c(S, 4, niter, nchain))
      mu_b0_output <- if (identical(intercept_prior$type, "community"))
        array(NA, dim = c(niter, nchain)) else NULL
      mu_theta_output <- if (identical(collection_intercept_prior$type, "community"))
        array(NA, dim = c(niter, nchain)) else NULL
    }

  }

  # run MCMC
  message("Running MCMC")

  # Every chain gets its own R seed and C++ seed, drawn here from R's RNG, so
  # set.seed() before the fit reproduces it and consecutive fits differ. A
  # chain's draws then depend on its two seeds alone, not on which chains ran
  # before it or in which process, which is what makes a fit identical whatever
  # MCMCparams$ncores says. The session's RNG is reset and restored on every
  # exit. Spec: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.
  seeds <- draw_chain_seeds(nchain)
  rng_saved <- save_session_rng()
  on.exit(restore_session_rng(rng_saved), add = TRUE)

  setup <- mget(chain_setup_fields(), envir = environment(), ifnotfound = list(NULL))
  set_chain_setup(setup)
  on.exit(clear_chain_setup(), add = TRUE)

  budget <- plan_chain_budget(nchain, ncores)
  handle <- NULL
  progress_file <- NULL
  if (budget$workers > 1L) {
    progress_file <- tempfile("occJSDM-progress-", fileext = ".log")
    # Registered before the cluster exists: handle is read at exit and a NULL
    # handle is a no-op, so no moment leaves a started cluster without a stop.
    on.exit(stop_chain_cluster(handle), add = TRUE)
    handle <- start_chain_cluster(budget$workers, budget$threads_per_worker, progress_file)
    check_worker_builds(handle)
    send_chain_setup(handle, setup)
    message(sprintf("Running %d chains on %d workers, %d RcppParallel thread%s each. Follow progress with: tail -f %s",
                    nchain, budget$workers, budget$threads_per_worker,
                    if (budget$threads_per_worker == 1L) "" else "s", progress_file))
  } else if (budget$ncores > 1L) {
    restore_threads <- pin_session_threads(budget$ncores)
    on.exit(restore_threads(), add = TRUE)
  }

  latent_sums <- NULL
  waic_parts <- vector("list", nchain)
  elapsed <- rep(NA_real_, nchain)
  batches <- split(seq_len(nchain), ceiling(seq_len(nchain) / budget$workers))
  for (batch in batches) {
    results <- run_chain_batch(handle, batch, seeds)
    for (res in results) {
      res <- stop_if_chain_failed(res)
      chain <- res$chain
      dr <- res$draws

      if(model %in% c("occupancy","two_stage")){
        beta_theta_output[,,,chain] <- dr$beta_theta
        theta0_output[,,chain] <- dr$theta0
        if (!is.null(mu_theta_output)) mu_theta_output[,chain] <- dr$mu_theta
        if(!summarisedLatentPresences){
          z_output[,,,chain] <- dr$z
          psi_output[,,,chain] <- dr$psi
        }
      }

      if(model == "two_stage"){
        p_output[,,,chain] <- dr$p
        q_output[,,,chain] <- dr$q
      }

      # save jsdm params
      {
        B0_output[,,chain] <- dr$B0
        B_output[,,,chain] <- dr$B
        G_output[,,,chain] <- dr$G
        A_output[,,,chain] <- dr$A
        C_output[,,,chain] <- dr$C
        Bs_output[,,,chain] <- dr$Bs
        Gs_output[,,,chain] <- dr$Gs
        As_output[,,,chain] <- dr$As
        Cs_output[,,,chain] <- dr$Cs
        U_output[,,,chain] <- dr$U
        L_output[,,,chain] <- dr$L
        sigmab_output[,chain] <- dr$sigmab
        sigmabs_output[,chain] <- dr$sigmabs
        sigmah_output[,chain] <- dr$sigmah
        idx_ls_output[,chain] <- dr$idx_ls
        varPart_output[,,,chain] <- dr$varPart
        tau_output[,,chain] <- dr$tau
        if (!is.null(mu_b0_output)) mu_b0_output[,chain] <- dr$mu_b0
      }

      latent_sums <- add_latent_sums(latent_sums, res$latent_sums)
      waic_parts[chain] <- list(res$waic)
      elapsed[chain] <- res$elapsed
      rm(dr)
    }
    rm(results, res)
  }

  # Every chain must have returned a result before anything is combined.
  missing_chains <- which(is.na(elapsed) | vapply(waic_parts, is.null, logical(1)))
  if (length(missing_chains) > 0L) {
    stop(sprintf("No result was returned for chain(s) %s", paste(missing_chains, collapse = ", ")),
         call. = FALSE)
  }
  for (chain in seq_len(nchain)) message(sprintf("Chain %d finished in %.1f s", chain, elapsed[chain]))

  # latent-state means, from the per-chain sums
  if (model %in% c("occupancy","two_stage")) {
    if (summarisedLatentPresences) {
      z_output_mean <- latent_sums$z / (niter * nchain)
      psi_output_mean <- latent_sums$psi / (niter * nchain)
    }
    theta_output_mean <- latent_sums$theta / (niter * nchain)
  }
  if (model == "two_stage") w_output_mean <- latent_sums$w / (niter * nchain)

  # compute WAIC
  {
    numIters <- niter * nchain
    for (chain in seq_len(nchain)) {
      for (nm in names(waic_parts[[chain]])) {
        part <- waic_parts[[chain]][[nm]]
        if (!is.null(part) && !isTRUE(part$n == niter)) {
          stop(sprintf("Current wAIC iter wrong: chain %d recorded %s %s draws, expected %s",
                       chain, deparse(part$n), nm, format(niter)), call. = FALSE)
        }
      }
    }
    WAIC_jsdm <- compute_waic(merge_waic_parts(lapply(waic_parts, `[[`, "jsdm")), numIters)

    if(model %in% c("occupancy","two_stage")){
      WAIC_w <- compute_waic(merge_waic_parts(lapply(waic_parts, `[[`, "w")), numIters)
    } else {
      WAIC_w <- 0
    }
    if(model %in% "two_stage"){
      WAIC_y <- compute_waic(merge_waic_parts(lapply(waic_parts, `[[`, "y")), numIters)
    } else {
      WAIC_y <- 0
    }

    WAIC <- WAIC_jsdm + WAIC_w + WAIC_y
  }

  # reparametrise factor coefficients
  {
    if(d > 0){

      list_UL_output_reparm <- reparamFactorModel(U_output, L_output)
      U_output <- list_UL_output_reparm$U_output
      L_output <- list_UL_output_reparm$L_output

    }

    if(ncov_psi > 0 & gt > 0){

      list_AC_output_reparm <- reparamFactorModel(A_output, C_output)
      A_output <- list_AC_output_reparm$U_output
      C_output <- list_AC_output_reparm$L_output

    }

  }

  jsdm_results_output <- list(
    "B0_output" = B0_output,
    "B_output" = B_output,
    "G_output" = G_output,
    "A_output" =  A_output,
    "C_output" = C_output,
    "Bs_output" = Bs_output,
    "Gs_output" = Gs_output,
    "As_output" = As_output,
    "Cs_output" = Cs_output,
    "U_output" = U_output,
    "L_output" = L_output,
    "sigmab_output" = sigmab_output,
    "sigmabs_output" = sigmabs_output,
    "sigmah_output" = sigmah_output,
    "idx_ls_output" = idx_ls_output,
    "varPart_output" = varPart_output,
    "tau_output" = tau_output
  )
  # Present only when b0_prior = "community", so default fits keep their shape.
  jsdm_results_output$mu_b0_output <- mu_b0_output

  results_output <- list(
    "jsdm_output" = jsdm_results_output,
    "beta_theta_output" = beta_theta_output,
    "p_output" = p_output,
    "q_output" = q_output,
    "theta0_output" = theta0_output,
    "WAIC" = WAIC
  )
  # Present only when theta_intercept_prior = "community", as for mu_b0_output.
  results_output$mu_theta_output <- mu_theta_output

  if(summarisedLatentPresences){
    results_output$z_output <- z_output_mean
    results_output$psi_output <- psi_output_mean
  } else {
    results_output$z_output <- z_output
    results_output$psi_output <- psi_output
  }
  results_output$w_output <- w_output_mean
  results_output$theta_output <- theta_output_mean

  # print diagnostics
  computeDiagnostics(results_output)

  infos <- list(
    "S" = S,
    "g" = g,
    "M" = M,
    "n" = n,
    "K" = K,
    "ps" = ps,
    "n_factors" = d,
    "list_idx" = list_idx,
    "data_info" = data_info,
    "speciesNames" = speciesNames,
    "primerNames" = primerNames,
    "siteNames" = siteNames,
    "ncov_theta" = ncov_theta,
    "ncov_psi" = ncov_psi,
    "OTU" = OTU,
    "threshold" = if (model %in% c("occupancy", "two_stage")) threshold else NULL,
    "X0_psi" = X0_psi,
    "list_Xs" = list_Xs,
    "list_X_psi_mat" = list_Xpsi_mat,
    "list_Xs_mat" = list_Xs_mat,
    "list_X_theta_mat" = list_X_theta_mat,
    "list_Tr_mat" = list_Tr_mat,
    "l_s_grid" = l_s_grid,
    "model" = model,
    "jsdmModel" = jsdmModel
  )
  infos$noise_prior <- noise_prior
  infos$spatial_sd_prior <- spatial_sd_prior
  infos$intercept_prior <- intercept_prior
  infos$collection_intercept_prior <- collection_intercept_prior
  infos$MCMCparams <- list(nchain = as.integer(nchain), nburn = as.integer(nburn),
                           niter = as.integer(niter), nthin = as.integer(nthin))
  infos$chain_seeds <- seeds
  infos$execution <- list(
    ncores = budget$ncores,
    workers = budget$workers,
    threads_per_worker = if (budget$thread_setting == "ambient") ambient_rcpp_threads() else budget$threads_per_worker,
    thread_setting = budget$thread_setting,
    blas_threads_requested = if (is.null(handle)) NULL else handle$blas_threads_requested,
    worker_pids = if (is.null(handle)) NULL else handle$pids,
    blas_declared = nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")),
    blas = unname(extSoftVersion()[["BLAS"]]),
    lapack = La_library(),
    elapsed = elapsed,
    progress_file = progress_file)

  list(
    "results_output" = results_output,
    "infos" = infos,
    "Tr" = Tr,
    "X_theta" = X_theta,
    "Xs" = Xs,
    "X_psi" = X_psi)

}

# The occupancy intercept prior is Normal(0, sd^2); it is not the collection intercept.
# The experimental b0_prior = "community" makes it Normal(mu_b0, sd^2) with
# mu_b0 ~ Normal(0, mu_b0_sd^2) estimated. The default list is unchanged.
read_intercept_prior <- function(priors) {
  sd <- get_param(priors, "sigma_b0", 1)
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0)
    stop("sigma_b0 must be a finite positive number")
  type <- get_param(priors, "b0_prior", "fixed")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("fixed", "community"))
    stop("b0_prior must be 'fixed' or 'community'")
  if (type == "fixed") {
    if ("mu_b0_sd" %in% names(priors))
      stop("mu_b0_sd requires b0_prior = 'community'")
    return(list(mean = 0, sd = unname(sd)))
  }
  mu_b0_sd <- get_param(priors, "mu_b0_sd", 2)
  if (!is.numeric(mu_b0_sd) || length(mu_b0_sd) != 1L || !is.finite(mu_b0_sd) ||
      mu_b0_sd <= 0)
    stop("mu_b0_sd must be a finite positive number")
  list(type = "community", sd = unname(sd), mu_b0_mean = 0,
       mu_b0_sd = unname(mu_b0_sd))
}

# The collection intercepts' prior is Normal(0, 1) (prior_beta_theta and
# prior_beta_theta_sd in runOccJSDM()). The experimental theta_intercept_prior
# = "community" makes it Normal(mu_theta, 1) with mu_theta ~ Normal(0,
# mu_theta_sd^2) estimated. The default returns NULL, so default fits record
# nothing new.
read_collection_intercept_prior <- function(priors, model) {
  type <- get_param(priors, "theta_intercept_prior", "fixed")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("fixed", "community"))
    stop("theta_intercept_prior must be 'fixed' or 'community'")
  if (type == "fixed") {
    if ("mu_theta_sd" %in% names(priors))
      stop("mu_theta_sd requires theta_intercept_prior = 'community'")
    return(NULL)
  }
  if (!model %in% c("occupancy", "two_stage"))
    stop("theta_intercept_prior = 'community' requires an occupancy or two-stage model")
  mu_theta_sd <- get_param(priors, "mu_theta_sd", 2)
  if (!is.numeric(mu_theta_sd) || length(mu_theta_sd) != 1L ||
      !is.finite(mu_theta_sd) || mu_theta_sd <= 0)
    stop("mu_theta_sd must be a finite positive number")
  list(type = "community", sd = 1, mu_theta_mean = 0,
       mu_theta_sd = unname(mu_theta_sd))
}

# This prior controls spatial coefficients, not continuous-response noise.
read_spatial_sd_prior <- function(priors,spatial) {
  supplied <- any(c("sigma_bs_prior","sigma_bs_scale") %in% names(priors))
  if (!spatial) {
    if (supplied) stop("sigma_bs prior settings require spatial covariates")
    return(NULL)
  }
  type <- get_param(priors,"sigma_bs_prior","inverse_gamma")
  if (!is.character(type) || length(type)!=1L || is.na(type) ||
      !type %in% c("half_cauchy","inverse_gamma"))
    stop("sigma_bs_prior must be 'half_cauchy' or 'inverse_gamma'")
  type <- unname(as.character(type))
  if (type=="inverse_gamma") {
    if ("sigma_bs_scale" %in% names(priors))
      stop("sigma_bs_scale requires sigma_bs_prior = 'half_cauchy'")
    return(list(type=type,shape=10,rate=1))
  }
  scale <- get_param(priors,"sigma_bs_scale",1)
  if (!is.numeric(scale) || length(scale)!=1L || !is.finite(scale) || scale<=0)
    stop("sigma_bs_scale must be a finite positive number")
  list(type=type,scale=unname(scale))
}

# Continuous-response noise priors are separate from the detection priors.
read_noise_prior <- function(priors) {
  type <- get_param(priors,"tau_prior","inverse_gamma")
  if (!is.character(type) || length(type)!=1L || is.na(type) ||
      !type %in% c("half_cauchy","inverse_gamma"))
    stop("tau_prior must be 'half_cauchy' or 'inverse_gamma'")
  type <- unname(as.character(type))
  positive_scalar <- function(name,default) {
    value <- get_param(priors,name,default)
    if (!is.numeric(value) || length(value)!=1L || !is.finite(value) || value<=0)
      stop(name," must be a finite positive number")
    value
  }
  if (type == "half_cauchy") {
    list(type=type,scale=positive_scalar("tau_scale",1))
  } else {
    list(type=type,shape=positive_scalar("a_tau",5),rate=positive_scalar("b_tau",5))
  }
}
