

logit <- function(x) log(x / (1-x))

logistic <- function(x) 1 / (1 + exp(-x))


#' simulateOccJSDMData
#'
#' Simulate data for any supported model type
#'
#' @param list_datasettings A list of survey dimensions:
#' \describe{
#'   \item{n}{Number of sites.}
#'   \item{S}{Number of species.}
#'   \item{g}{Number of measured species traits (the columns of the returned
#'   \code{traits} table).}
#'   \item{M}{Number of samples at each site: a vector of length \code{n},
#'   one entry per site. A single number is not recycled; use
#'   \code{rep(m, n)} for the same number everywhere. Used by the
#'   \code{"occupancy"} and \code{"two_stage"} models.}
#'   \item{P}{Number of primers, a single number. Every sample is analysed
#'   with all \code{P} primers. Used by \code{"two_stage"} only.}
#'   \item{K}{Number of PCR replicates for each sample and primer: a vector
#'   of length \code{P * sum(M)}, ordered by sample (sites in order, samples
#'   in order within a site) and by primer within a sample. Used by
#'   \code{"two_stage"} only.}
#'   \item{ncov_psi}{Number of occupancy (environmental) covariates, at
#'   least 1. Each is drawn for every site from a normal distribution with
#'   mean 0 and standard deviation 10. When \code{ncov_psi} is 3 or more,
#'   \code{round(0.2 * ncov_psi)} of them are instead categorical, with
#'   levels \code{"1"}, \code{"2"} and \code{"3"}. Numeric covariates are
#'   centred and scaled before they enter the model, so each coefficient is
#'   the change in the occupancy linear predictor per standard deviation of
#'   its covariate; a categorical covariate enters as indicator columns.}
#'   \item{ncov_theta}{Number of collection covariates, drawn from a
#'   standard normal distribution independently for every sample. Used by
#'   the \code{"occupancy"} and \code{"two_stage"} models.}
#' }
#' @param list_params A list of collection and laboratory settings. It is
#'   not used by the \code{"binary"} and \code{"continuous"} models and can
#'   then be an empty list.
#' \describe{
#'   \item{p}{A \code{P} by \code{S} matrix: the probability that one PCR
#'   replicate with a given primer detects a species whose DNA is in the
#'   sample. \code{"two_stage"} only.}
#'   \item{q}{A \code{P} by \code{S} matrix: the probability that one PCR
#'   replicate with a given primer gives a false-positive detection of a
#'   species whose DNA is not in the sample. \code{"two_stage"} only.}
#'   \item{theta0}{A vector of length \code{S}: the probability that a
#'   sample collects a species' DNA although the species does not occupy
#'   the site (contamination in the field).}
#'   \item{theta_baseline}{A vector of length \code{S}, or one number used
#'   for every species: the probability that a sample collects the DNA of a
#'   species that occupies the site, when the collection covariates are at
#'   their mean of 0. Its logit is the intercept row of
#'   \code{beta_theta_true}; the slopes on the collection covariates are
#'   drawn from -1, 0 and 1.}
#'   \item{mu1, sigma1, mu0, sigma0}{Optional read-count settings for
#'   \code{"two_stage"} (see Details). If omitted they default to
#'   \code{mu1 = 5}, \code{sigma1 = 1}, \code{mu0 = 1.5} and
#'   \code{sigma0 = 1}.}
#' }
#' @param list_jsdmParams A list of settings for the occupancy part of the
#'   model, shared by all four model types:
#' \describe{
#'   \item{gt}{Number of unmeasured traits: the columns of \code{A} in the
#'   coefficient formula under Details.}
#'   \item{d}{Number of hidden factors. Their scores \code{U} (\code{n} by
#'   \code{d}) are normal with standard deviation \code{sigma_h}; their
#'   loadings \code{L} (\code{d} by \code{S}) are drawn from -1, 0 and 1,
#'   with 1 on the diagonal and 0 below it. The term \code{U \%*\% L} makes
#'   species occur together, or apart, more than the environment explains.}
#'   \item{ds}{How the spatial field is shared between species; used only
#'   when \code{useSpatField = TRUE}. \code{ds = 0} draws one field that
#'   every species shares. \code{ds > 0} draws one field per species, with
#'   the fields correlated across species through a correlation matrix of
#'   rank \code{ds}: they are built from \code{ds} underlying patterns, so
#'   some species do well in the same patches and others in opposite ones.}
#'   \item{sigma_b}{Standard deviation of \code{Bt}, the part of each
#'   species' environmental coefficients that its traits do not explain.}
#'   \item{sigma_bs, sigma_ts}{Accepted but have no effect on the simulated
#'   data. The simulator has no residual spatial coefficients to draw with
#'   \code{sigma_bs} (it only stores the value in \code{jsdmParams_true}),
#'   and it passes \code{sigma_ts} on without using it. Either can be left
#'   out.}
#'   \item{sigma_h}{Standard deviation of the hidden-factor scores
#'   \code{U}.}
#'   \item{sigma_s}{Scale of the spatial field. The field's covariance
#'   between two sites at distance \code{dist} is
#'   \code{sigma_s * exp(-dist^2 / (2 * l_s^2))}, so \code{sigma_s} is the
#'   field's marginal variance at any one site.}
#'   \item{l_s}{Length scale of the spatial field, in the units of the site
#'   coordinates. The coordinates are drawn uniformly on the unit square,
#'   one pair per site, whether or not the field is switched on. Small
#'   values give small patches; values near 1 or more make the field nearly
#'   flat across the sites.}
#'   \item{tau}{A vector of length \code{S}: the residual standard deviation
#'   of each species' response. Required by \code{"continuous"} and not
#'   used by the other models.}
#'   \item{rnb}{Accepted but has no effect: it would set the dispersion of
#'   a count model, which this simulator does not offer.}
#'   \item{useSpatField}{\code{TRUE} or \code{FALSE}; required. The spatial
#'   field is added to the occupancy linear predictor only when it is
#'   \code{TRUE}. When it is \code{FALSE}, \code{ds}, \code{sigma_s} and
#'   \code{l_s} have no effect.}
#' }
#' @param model The kind of data to simulate: \code{"binary"},
#'   \code{"continuous"}, \code{"occupancy"} or \code{"two_stage"}. See
#'   Details.
#'
#' @details
#' General-purpose simulation function that supports multiple model types.
#' Every model starts from the same occupancy part. For each site and
#' species, the linear predictor \code{eta} is an intercept
#' \code{B0} (standard normal), plus the site's covariates times the
#' species' coefficients \code{B}, plus the hidden-factor term
#' \code{U \%*\% L}, plus the spatial field when \code{useSpatField = TRUE}.
#' The four models then differ in what they return as \code{OTU}:
#' \describe{
#'   \item{\code{"binary"}}{One row per site: the true occupancy state
#'   \code{z}, drawn as Bernoulli with probability \code{plogis(eta)}.}
#'   \item{\code{"continuous"}}{One row per site: \code{eta} plus normal
#'   noise with standard deviation \code{tau}.}
#'   \item{\code{"occupancy"}}{One row per sample: whether the sample
#'   collected the species' DNA (\code{w}). Given \code{z = 1}, this is
#'   Bernoulli with probability \code{theta}, where \code{logit(theta)} is
#'   an intercept plus the sample's collection covariates times their
#'   slopes (the rows of \code{beta_theta_true}, intercept first); given
#'   \code{z = 0}, the probability is \code{theta0}.}
#'   \item{\code{"two_stage"}}{One row per PCR replicate: read counts. Each
#'   replicate detects the species with probability \code{p} when
#'   \code{w = 1}, and gives a false positive with probability \code{q}
#'   when \code{w = 0}.}
#' }
#'
#' Read counts (\code{"two_stage"} only). Given a true detection,
#' \code{log(y + 1)} is normal with mean \code{mu1} and standard deviation
#' \code{sigma1}; given a false-positive detection, it is normal with mean
#' \code{mu0} and standard deviation \code{sigma0}. A replicate without a
#' detection has \code{y = 0}. For a normal draw \code{v}, the count is
#' \code{round(exp(v) - 1)} with negative values set to 0, so a weak
#' detection can also give 0 reads.
#'
#' How the environmental coefficients are built. Species' coefficients
#' come partly from their traits:
#' \code{B = t(Tr \%*\% G + A \%*\% C + Bt)},
#' where \code{Tr} is the \code{S} by \code{g} table of measured traits
#' (standard normal, returned as \code{traits}), \code{G} is the effect of
#' each measured trait on each covariate's coefficient, \code{A} is an
#' \code{S} by \code{gt} matrix of unmeasured traits, \code{C} is their
#' effect, and \code{Bt} is the part the traits do not explain. \code{G}
#' and the free entries of \code{C} are drawn from -1, 0 and 1 (\code{C}
#' has 1 on its diagonal and 0 below it), \code{A} is standard normal, and
#' \code{Bt} is normal with standard deviation \code{sigma_b}. \code{B} has
#' one row per covariate column and one column per species.
#'
#' The spatial field. Each site has a pair of coordinates drawn uniformly
#' on the unit square, returned as \code{Xs.1} and \code{Xs.2} in
#' \code{info} whether or not the field is used. With
#' \code{useSpatField = TRUE}, a Gaussian field with the squared-exponential
#' covariance \code{sigma_s * exp(-dist^2 / (2 * l_s^2))} is drawn at the
#' sites and added to \code{eta}: one field shared by every species when
#' \code{ds = 0}, or one field per species, correlated across species,
#' when \code{ds > 0}. To make the spatial signal stronger, raise
#' \code{sigma_s}, or lower \code{sigma_h} to weaken the hidden factors.
#' \code{l_s} sets the patch size; values near 1 or more flatten the field
#' across the unit square and so weaken it.
#'
#' Variance partitioning. \code{true_params$jsdmParams_true$varPart} has
#' one row per species. Its \code{Environmental}, \code{Spatial} and
#' \code{Biotic} (hidden-factor) columns are relative contributions that
#' add up to one for each species; they say which source matters most, but
#' they are not pieces of \code{Total}. \code{Total} is the standard
#' deviation across sites of the species' true occupancy probability
#' \code{plogis(eta)} (of \code{eta} itself for \code{"continuous"}). The
#' \code{Spatial} column is 0 when \code{useSpatField = FALSE}.
#'
#' For a worked example, see
#' \code{vignette("simulateOccJSDMData", package = "occJSDM")}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{true_params}{The values that generated the data, to compare with
#'   a fit. \code{jsdmParams_true} is a list that includes \code{B0},
#'   \code{B}, \code{G}, \code{A}, \code{C}, \code{Bt}, \code{U}, \code{L},
#'   \code{spatField} (\code{n} by \code{S}; zero when the field is off),
#'   \code{eta} (\code{n} by \code{S}) and \code{varPart} (see Details).
#'   For \code{"occupancy"} and \code{"two_stage"} there are also
#'   \code{beta_theta_true} (collection coefficients, intercept row first)
#'   and \code{z_true} (the \code{n} by \code{S} true occupancy states). For
#'   \code{"two_stage"} there are also \code{w_true} (one row per sample:
#'   whether the sample collected the species' DNA), \code{p_true} and
#'   \code{q_true}.}
#'   \item{data_list}{The simulated survey, in the form
#'   \code{runOccJSDM()} takes as its \code{data} argument: \code{info}, a
#'   data frame with one row per observation (the \code{Site},
#'   \code{Sample} and \code{Primer} ids the model has, the occupancy
#'   covariates \code{X_psi.EnvCov.1}, \code{X_psi.EnvCov.2}, and so on, on
#'   their original scale, the coordinates \code{Xs.1} and \code{Xs.2},
#'   and, for \code{"occupancy"} and \code{"two_stage"}, the collection
#'   covariates \code{X_theta}, or \code{X_theta.1}, \code{X_theta.2}, and
#'   so on when there are several); \code{OTU}, a
#'   matrix with one column per species, named \code{OTU_1},
#'   \code{OTU_2}, and so on; and \code{traits}, the \code{S} by \code{g}
#'   trait table with rows named by species and columns \code{Trait_1},
#'   \code{Trait_2}, and so on.}
#' }
#'
#' @export
#'
simulateOccJSDMData <- function(list_datasettings,
                                list_params,
                                list_jsdmParams,
                                model){

  if(!(model %in% c("binary","occupancy","continuous","two_stage"))){
    stop("The model used is not supported")
  }

  # read data settings
  {
    n <- list_datasettings$n
    S <- list_datasettings$S
    g <- list_datasettings$g
    M <- list_datasettings$M
    P <- list_datasettings$P
    K <- list_datasettings$K
    ncov_psi <- list_datasettings$ncov_psi
    ncov_theta <- list_datasettings$ncov_theta
  }

  # read jsdm params
  {
    gt <- list_jsdmParams$gt
    d <- list_jsdmParams$d
    ds <- list_jsdmParams$ds
    sigma_b <- list_jsdmParams$sigma_b
    sigma_bs <- list_jsdmParams$sigma_bs
    sigma_ts <- list_jsdmParams$sigma_ts
    sigma_h <- list_jsdmParams$sigma_h
    sigma_s <- list_jsdmParams$sigma_s
    l_s <- list_jsdmParams$l_s
    tau <- list_jsdmParams$tau
    rnb <- list_jsdmParams$rnb
    useSpatField <- list_jsdmParams$useSpatField
  }

  # read param settings
  {
    p <- list_params$p
    q <- list_params$q
    theta0 <- list_params$theta0
    theta_baseline <- list_params$theta_baseline
  }

  if(model %in% c("binary","occupancy","two_stage")){
    jSDMsimModel <- "binary"
  } else if(model == "continuous"){
    jSDMsimModel <- "continuous"
  }

  list_simjSDMData <- simulateData(
    n, S, ncov_psi,
    g, gt, d, tau, rnb, ds, n,
    sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
    useSpatField = useSpatField, usingSplines = F, model = jSDMsimModel)

  # data
  {
    z <- list_simjSDMData$data$z
    X_psi <- list_simjSDMData$data$X
    Tr <- list_simjSDMData$data$Tr
    Xs <- list_simjSDMData$data$Xs
  }

  # create data structure (for occupancy or two-stage model)
  {
    N <- sum(M)

    if(model %in% c("occupancy","two_stage")){
      N2 <- P * N
      N3 <- sum(K)

      maxP <- P
      P <- rep(P, N)

      sumM <- c(0, cumsum(M)[-n])
      sumP <- c(0, cumsum(P)[-N])
      sumK <- c(0, cumsum(K)[-N2])


      list_idx <- createDataIdx(n, M, P, K, model == "two_stage",
                                primerId_p = rep(1:maxP, times = N))
      idx_z_w <- list_idx$idx_z_w
      idx_z_k <- list_idx$idx_z_k
      idx_w_p <- list_idx$idx_w_p
      idx_z_p <- list_idx$idx_z_p
      idx_w_k <- list_idx$idx_w_k
      idx_p_k <- list_idx$idx_p_k
    }

  }

  if(model %in% c("occupancy","two_stage")){
    X_theta <- cbind(1, matrix(rnorm(N * ncov_theta), N, ncov_theta))

    beta_theta_true <- matrix(sample(c(-1,1,0), (ncov_theta + 1) * S, replace = T), ncov_theta + 1, S)
    beta_theta_true[1,] <- logit(theta_baseline)
    theta_true <- logistic(X_theta %*% beta_theta_true)

    theta0_true <- theta0

    z_rep <- z[idx_z_w,]
    w <- sapply(1:S, function(s){
      sapply(1:N, function(i){
        if(z_rep[i,s] == 1){
          rbinom(1, 1, theta_true[i,s])
        } else {
          rbinom(1, 1, theta0_true[s])
        }
      })
    })

    if(model == "two_stage"){
      p_true <- p
      q_true <- q

      y <- matrix(NA, N3, S)
      w_rep <- w[idx_w_k,]

      cimk_true <- sapply(1:S, function(s){
        sapply(1:N3, function(i){
          if(w_rep[i,s] == 1){
            rbinom(1, 1, p_true[idx_p_k[i],s])
          } else {
            2 * rbinom(1, 1, q_true[idx_p_k[i],s])
          }
        })
      })

      # Simulate read counts; runOccJSDM() converts them to binary detections
      # at a positive threshold rather than fitting continuous intensities. Given
      # a true detection (cimk_true == 1), log(y + 1) ~ Normal(mu1, sigma1); given
      # a false-positive/contamination detection (cimk_true == 2),
      # log(y + 1) ~ Normal(mu0, sigma0). No reads (cimk_true == 0) gives y = 0.
      mu1_true <- get_param(list_params, "mu1", 5)
      sigma1_true <- get_param(list_params, "sigma1", 1)
      mu0_true <- get_param(list_params, "mu0", 1.5)
      sigma0_true <- get_param(list_params, "sigma0", 1)

      logy1 <- matrix(0, N3, S)
      logy1[cimk_true == 1] <- rnorm(sum(cimk_true == 1), mu1_true, sigma1_true)
      logy1[cimk_true == 2] <- rnorm(sum(cimk_true == 2), mu0_true, sigma0_true)

      y <- round(exp(logy1) - 1)
      y[y < 0] <- 0
    }

  }

  if(model %in% c("continuous","binary")){
    y <- z
  } else if(model == "occupancy"){
    y <- w
  } else if(model == "two_stage"){
    y <- y
  }

  speciesNames <- paste0("OTU_", 1:S)

  colnames(y) <- speciesNames
  rownames(Tr) <- speciesNames
  colnames(Tr) <- paste0("Trait_", seq_len(ncol(Tr)))

  if(model %in% c("binary","continuous")){
    data_info <- data.frame(
      X_psi = X_psi,
      Xs = Xs
    )
  } else if(model %in% c("occupancy")){
    data_info <- data.frame(
      Site = idx_z_w,
      X_psi = X_psi[idx_z_w,],
      Xs = Xs[idx_z_w,],
      X_theta = X_theta[,-1]
    )
  } else if (model == "two_stage"){
    data_info <- data.frame(
      Site = idx_z_k,
      Sample = idx_w_k,
      Primer = idx_p_k,
      X_psi = X_psi[idx_z_k,],
      Xs = Xs[idx_z_k,],
      X_theta = X_theta[idx_w_k,-1]
    )
  }

  data <- list(info = data_info,
               OTU = y,
               traits = Tr)

  true_params <- list(
    "jsdmParams_true" = list_simjSDMData$trueParams
  )

  if(model %in% c("occupancy","two_stage")){
    true_params$beta_theta_true <- beta_theta_true
    true_params$z_true <- z
  }

  if(model == "two_stage"){
    true_params$w_true <- w
    true_params$p_true <- p_true
    true_params$q_true <- q_true
  }

  list(true_params = true_params,
       data_list = data)
}
