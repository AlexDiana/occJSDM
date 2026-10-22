#' Return posterior draws with chains kept separate
#'
#' Extract one occupancy or detection parameter from a fitted model, ready for
#' [plotTraceplot()] or summaries calculated separately for each chain.
#'
#' @param fitmodel A model object returned by [runOccJSDM()].
#' @param parameter A single parameter name: \code{"beta0_psi"} (occupancy
#'   intercepts), \code{"beta_psi"} (occupancy covariate effects),
#'   \code{"beta_theta"} (collection intercepts and covariate effects),
#'   \code{"p"} (PCR true-positive probabilities), \code{"q"} (PCR
#'   false-positive probabilities), or \code{"theta0"} (field-contamination
#'   probabilities).
#' @param species Optional character vector of species names. \code{NULL}
#'   returns every species, in fitted order.
#' @param covariate Optional character vector of fitted covariate names, for
#'   \code{"beta_psi"} or \code{"beta_theta"}. Use \code{"(Intercept)"}
#'   for the collection intercept. \code{NULL} returns every coefficient.
#' @param primer Optional vector of primer identifiers, for \code{"p"} or
#'   \code{"q"}. Numeric identifiers and their character forms are equivalent;
#'   they identify primers, not positions in an array. \code{NULL} returns every
#'   primer.
#'
#' @details
#' Selections preserve the requested order and never drop dimensions, including
#' when selecting one species or when the fit has only one chain. Unknown names,
#' selectors that do not apply to the parameter, and parameters without saved
#' draws produce an error. A plain JSDM fit has no detection-parameter draws.
#'
#' Coefficients retain their fitted scale: for binary, occupancy and two-stage
#' models they are on the log-odds scale, and for continuous models occupancy
#' coefficients are on the response scale. Numeric covariates use the fitting
#' function's standardization; factor coefficients use its fitted contrasts.
#' The \code{p}, \code{q} and \code{theta0} draws are probabilities.
#'
#' The iteration axis indexes saved draws after burn-in, from \code{1} to the
#' number retained per chain. It does not give the original sampler iteration
#' numbers. Draws are neither pooled nor thinned by this function. Posterior
#' means of latent states are not draws and are not supported here.
#'
#' @return A labelled numeric array with four dimensions:
#'   \code{[coefficient or primer, species, saved iteration, chain]}.
#'   For \code{beta0_psi} and \code{theta0}, the first dimension has length one
#'   and is labelled with the parameter name. Dimension labels are read by
#'   [plotTraceplot()] unless explicitly overridden.
#'
#' @examples
#' data(sampleresults)
#' draws <- returnPosteriorDraws(sampleresults, "theta0")
#' # Average over saved iterations, keeping species and chains separate.
#' apply(draws, c(2, 4), mean)
#' plotTraceplot(draws, param_name = "Field-contamination probability")
#'
#' @md
#' @export
returnPosteriorDraws <- function(fitmodel, parameter, species = NULL,
                                 covariate = NULL, primer = NULL) {
  supported <- c("beta0_psi", "beta_psi", "beta_theta", "p", "q", "theta0")
  if (!is.character(parameter) || length(parameter) != 1L ||
      is.na(parameter) || !parameter %in% supported) {
    stop("'parameter' must be a single name from: ", paste(supported, collapse = ", "))
  }
  if (!is.list(fitmodel) || !is.list(fitmodel$results_output)) {
    stop("'fitmodel' must contain results_output from runOccJSDM()")
  }
  has_covariate <- parameter %in% c("beta_psi", "beta_theta")
  has_primer <- parameter %in% c("p", "q")
  if (!is.null(covariate) && !has_covariate) {
    stop("'covariate' applies only to beta_psi and beta_theta")
  }
  if (!is.null(primer) && !has_primer) {
    stop("'primer' applies only to p and q")
  }

  results <- fitmodel$results_output
  draws <- switch(parameter,
    beta0_psi = results$jsdm_output$B0_output,
    beta_psi = results$jsdm_output$B_output,
    beta_theta = results$beta_theta_output,
    p = results$p_output,
    q = results$q_output,
    theta0 = results$theta0_output
  )
  if (is.null(draws) || length(draws) == 0L) {
    stop("Parameter '", parameter, "' has no saved draws in this fit")
  }
  species_only <- parameter %in% c("beta0_psi", "theta0")
  if (!is.numeric(draws) || length(dim(draws)) != if (species_only) 3L else 4L) {
    stop("Saved draws for '", parameter, "' have an unexpected array shape")
  }
  if (species_only) draws <- array(draws, dim = c(1L, dim(draws)))

  first_names <- switch(parameter,
    beta_psi = colnames(fitmodel$X_psi),
    beta_theta = colnames(fitmodel$X_theta),
    p = as.character(fitmodel$infos$primerNames),
    q = as.character(fitmodel$infos$primerNames),
    parameter
  )
  species_names <- as.character(fitmodel$infos$speciesNames)
  if (length(first_names) != dim(draws)[1L] ||
      length(species_names) != dim(draws)[2L]) {
    stop("Parameter '", parameter, "' has missing or inconsistent fitted labels")
  }
  dimnames(draws) <- list(
    first_names, species_names, as.character(seq_len(dim(draws)[3L])),
    paste0("chain_", seq_len(dim(draws)[4L])))
  names(dimnames(draws)) <- c(if (has_covariate) "covariate" else if (has_primer)
    "primer" else "parameter", "species", "iteration", "chain")

  select_names <- function(selection, available, argument, numeric_ids = FALSE) {
    if (is.null(selection)) return(seq_along(available))
    if (!(is.character(selection) || (numeric_ids && is.numeric(selection))) ||
        length(selection) == 0L || anyNA(selection)) {
      stop("'", argument, "' must contain one or more fitted names")
    }
    indices <- match(as.character(selection), available)
    if (anyNA(indices)) {
      stop("'", argument, "' not found: ",
           paste(selection[is.na(indices)], collapse = ", "),
           ". Available ", argument, " names: ", paste(available, collapse = ", "))
    }
    indices
  }
  first <- select_names(if (has_covariate) covariate else primer,
                        first_names, if (has_covariate) "covariate" else "primer",
                        numeric_ids = has_primer)
  selected_species <- select_names(species, species_names, "species")
  draws[first, selected_species, , , drop = FALSE]
}
