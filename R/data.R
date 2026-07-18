#' Simulated example dataset for runOccJSDM()
#'
#' A simulated eDNA metabarcoding dataset, in the \code{info}/\code{OTU} list
#' format expected by \code{\link{runOccJSDM}}. Generated via
#' \code{\link{simulateOccJSDMData}} and converted with
#' \code{\link{toRunOccJSDMFormat}}.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{info}{A data.frame with one row per PCR replicate (1800 rows in
#'     the shipped example) and columns:
#'     \describe{
#'       \item{Site}{Integer site id.}
#'       \item{Sample}{Integer sample id (unique across sites).}
#'       \item{Primer}{Integer primer id (1:P within each sample).}
#'       \item{X_psi.*}{Occupancy covariates (site-level, repeated across
#'         each site's replicates).}
#'       \item{X_theta}{Collection covariate(s) (sample-level, repeated
#'         across each sample's replicates). \code{runOccJSDM()} adds its
#'         own intercept internally, so no intercept column is included
#'         here.}
#'       \item{Xs.*}{Site coordinates (site-level, repeated across each
#'         site's replicates), for use as \code{spatCovariates} in
#'         \code{runOccJSDM()}.}
#'     }
#'   }
#'   \item{OTU}{A numeric matrix of read counts, with \code{nrow(OTU) ==
#'     nrow(info)} and one column per species (named \code{OTU_1},
#'     \code{OTU_2}, ...).}
#'   \item{spatCovariates}{A character vector of the \code{Xs.*} column
#'     names in \code{info} (\code{c("Xs.1", "Xs.2")}), ready to pass as the
#'     \code{spatCovariates} argument of \code{runOccJSDM()}.}
#'   \item{traitsMatrix}{A species (rows) by trait (columns) matrix, with
#'     row names matching \code{colnames(OTU)} and column names
#'     \code{Trait_1}, \code{Trait_2}, ..., ready to pass as the
#'     \code{traitsMatrix} argument of \code{runOccJSDM()}.}
#' }
#'
#' @seealso \code{\link{simulateOccJSDMData}}, \code{\link{toRunOccJSDMFormat}},
#'   \code{\link{runOccJSDM}}
#'
#' @examples
#' str(sampledata, max.level = 1)
#' head(sampledata$info)
"sampledata"

#' Example fitted model for the occJSDM vignette
#'
#' Output of \code{\link{runOccJSDM}} fit to \code{\link{sampledata}} (10
#' species, 100 sites, up to 3 primers per sample), with
#' \code{MCMCparams = list(nchain = 2, nburn = 5000, niter = 5000, nthin =
#' 1)}. Shipped so that \code{vignettes/occJSDM.Rmd} can load a precomputed
#' fit rather than re-running \code{runOccJSDM()} (which takes several
#' minutes) at vignette build time.
#'
#' @format A list with six elements:
#' \describe{
#'   \item{results_output}{A list of posterior summaries/samples, including
#'     \code{jsdm_output} (JSDM coefficient posteriors), \code{beta_theta_output},
#'     \code{p_output}, \code{q_output}, \code{theta0_output} (detection-process
#'     posteriors, dimensioned \code{[primer, species, niter, nchain]} or
#'     similar), \code{WAIC}, and posterior means \code{z_output} (latent
#'     occupancy), \code{w_output} (latent collection), \code{psi_output}
#'     (occupancy probabilities), and \code{theta_output} (collection
#'     probabilities).}
#'   \item{infos}{A list of metadata about the fitted model (e.g. \code{S},
#'     \code{n}, \code{M}, \code{K}, \code{speciesNames}, \code{primerNames},
#'     \code{siteNames}, \code{ncov_psi}, \code{ncov_theta}, the \code{OTU}
#'     matrix used for fitting, and standardisation info for \code{X_psi}/
#'     \code{Xs}) used by the plotting/output helpers (e.g.
#'     \code{\link{plotOccupancyRates}}, \code{\link{plotDetectionRates}}).}
#'   \item{Tr}{The species trait matrix (species by trait) used to fit the
#'     model.}
#'   \item{X_theta}{The collection covariate design matrix (including
#'     intercept) used to fit the model.}
#'   \item{X_s}{The site coordinates used for the spatial field.}
#'   \item{X_psi}{The (standardised) occupancy covariate design matrix used
#'     to fit the model.}
#' }
#'
#' @seealso \code{\link{runOccJSDM}}, \code{\link{sampledata}}
#'
#' @examples
#' str(sampleresults, max.level = 1)
#' plotDetectionRates(sampleresults, idx_species = 1:5)
"sampleresults"
