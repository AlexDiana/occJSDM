#' Simulated example dataset for runOccJSDM()
#'
#' A simulated eDNA metabarcoding dataset, in the \code{info}/\code{OTU}/
#' \code{traits} list format expected by \code{\link{runOccJSDM}}. Generated
#' via \code{\link{simulateOccJSDMData}(model = "two_stage")}, which returns
#' \code{data_list} already in this ready-to-use shape.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{info}{A data.frame with one row per PCR replicate (1800 rows in
#'     the shipped example) and columns:
#'     \describe{
#'       \item{Site}{Integer site id.}
#'       \item{Sample}{Integer sample id. Unique across the whole survey
#'         in the shipped example; \code{runOccJSDM()} identifies a sample
#'         by its (Site, Sample) pair, so ids may also restart at each site.}
#'       \item{Primer}{Integer primer id (1:P within each sample).}
#'       \item{X_psi.*}{Occupancy covariates (site-level, repeated across
#'         each site's replicates); named \code{X_psi.EnvCov.1},
#'         \code{X_psi.EnvCov.2} in the shipped example.}
#'       \item{X_theta.*}{Collection covariate(s) (sample-level, repeated
#'         across each sample's replicates). \code{runOccJSDM()} adds its
#'         own intercept internally, so no intercept column is included
#'         here.}
#'       \item{Xs.*}{Site coordinates (site-level, repeated across each
#'         site's replicates), for use as \code{spatCovariates} in
#'         \code{runOccJSDM()} (e.g. \code{c("Xs.1", "Xs.2")}).}
#'     }
#'   }
#'   \item{OTU}{A numeric matrix of read counts, with \code{nrow(OTU) ==
#'     nrow(info)} and one column per species (named \code{OTU_1},
#'     \code{OTU_2}, ...).}
#'   \item{traits}{A species (rows) by trait (columns) matrix, with row
#'     names matching \code{colnames(OTU)} and column names \code{Trait_1},
#'     \code{Trait_2}, .... \code{runOccJSDM()} reads this directly as
#'     \code{data$traits} (the \code{traitsMatrix} function argument is
#'     unused).}
#' }
#'
#' @return A list with elements \code{info} (data.frame), \code{OTU}
#'   (matrix), and \code{traits} (matrix), as described under Format.
#'
#' @seealso \code{\link{simulateOccJSDMData}}, \code{\link{runOccJSDM}}
#'
#' @examples
#' str(sampledata, max.level = 1)
#' head(sampledata$info)
"sampledata"

#' Example fitted model for the occJSDM vignette
#'
#' Output of \code{\link{runOccJSDM}} fit to \code{\link{sampledata}} (10
#' species, 100 sites, up to 3 primers per sample), with
#' \code{MCMCparams = list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1,
#' ncores = 2)}, seed \code{20261007} and R 4.6.1, which is the quickstart's
#' live-fit call with \code{ncores = 2} added. Shipped so that
#' \code{vignettes/occJSDM-00-quickstart.Rmd} can load a precomputed fit
#' rather than re-running \code{runOccJSDM()} (which takes several minutes) at
#' vignette build time.
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
#'   \item{Xs}{The site coordinates used for the spatial field.}
#'   \item{X_psi}{The (standardised) occupancy covariate design matrix used
#'     to fit the model.}
#' }
#'
#' @return A list with elements \code{results_output}, \code{infos},
#'   \code{Tr}, \code{X_theta}, \code{Xs}, and \code{X_psi}, as described
#'   under Format -- i.e. the object returned by \code{\link{runOccJSDM}}.
#'
#' @seealso \code{\link{runOccJSDM}}, \code{\link{sampledata}}
#'
#' @examples
#' str(sampleresults, max.level = 1)
#' plotDetectionRates(sampleresults, idx_species = 1:5)
"sampleresults"

#' Simulated river eDNA survey for the quickstart
#'
#' A simulated two-stage eDNA metabarcoding survey of 20 vertebrate species
#' at 100 river sites, modelled on the Gaoligongshan aquatic eDNA study of
#' Ji et al. (2025, \emph{Ecology Letters}). Each site has 3 water samples,
#' and each sample was tested with 2 primers, with 3 PCR replicates per
#' primer. The species are not real: each is named by class and order with a
#' running number, so that results can be grouped by taxon. Generated by
#' \code{dev/release/2026-11-01-quickstart-data/simulate-quickstartdata.R} in the source
#' repository, which plants known effects; the truth is not shipped.
#'
#' @format A list with three elements, in the format expected by
#'   \code{\link{runOccJSDM}}:
#' \describe{
#'   \item{info}{A data.frame with one row per PCR reaction (1800 rows) and
#'     columns:
#'     \describe{
#'       \item{Site}{Integer site id, 1 to 100.}
#'       \item{Sample}{Integer water-sample id, unique across the survey.}
#'       \item{Primer}{Integer primer id, 1 or 2.}
#'       \item{elevation}{Site elevation in metres (occupancy covariate).}
#'       \item{hmi}{Human Modification Index of the site, 0 to 1 (occupancy
#'         covariate).}
#'       \item{protected}{Factor with levels \code{"not_protected"} and
#'         \code{"protected"}: whether the stream rises inside a protected
#'         area (occupancy covariate). \code{"not_protected"} is the
#'         reference level.}
#'       \item{water_volume}{Litres of water filtered for the sample
#'         (collection covariate).}
#'     }
#'   }
#'   \item{OTU}{A numeric matrix of read counts, one row per row of
#'     \code{info} and one column per species (for example
#'     \code{Mammalia_Carnivora_1}).}
#'   \item{traits}{A species by trait matrix with row names matching
#'     \code{colnames(OTU)} and columns \code{log10_body_mass} (log10 of body
#'     mass in kg) and
#'     \code{red_listed} (1 if the species is on the IUCN Red List, else 0).}
#' }
#'
#' @seealso \code{\link{quickstartresults}}, \code{\link{runOccJSDM}}
#'
#' @examples
#' str(quickstartdata, max.level = 1)
#' head(quickstartdata$info)
"quickstartdata"

#' Example fitted model for the quickstart
#'
#' Output of \code{\link{runOccJSDM}} fit to \code{\link{quickstartdata}}
#' without a spatial field and with the default priors:
#' \preformatted{
#' set.seed(20261101)
#' runOccJSDM(
#'   data = quickstartdata,
#'   listParams = list(n_factors = 2),
#'   threshold = 1,
#'   occCovariates = c("elevation", "hmi", "protected"),
#'   collCovariates = c("water_volume"),
#'   MCMCparams = list(nchain = 2, nburn = 2000, niter = 2000, nthin = 2,
#'                     ncores = 2),
#'   summarisedLatentPresences = TRUE
#' )
#' }
#' run with R 4.6.1 and Apple's Accelerate BLAS. The chains are deliberately
#' short enough that a few parameters miss the quickstart's convergence rule
#' of thumb, so that the quickstart can show what that looks like. Shipped so
#' that \code{vignettes/occJSDM-00-quickstart.Rmd} can load a precomputed fit.
#' The script that made it, with its record, is in
#' \code{dev/release/2026-11-01-quickstart-data/} in the source repository.
#'
#' @format A list as returned by \code{\link{runOccJSDM}}.
#'
#' @seealso \code{\link{quickstartdata}}, \code{\link{returnConvergenceDiagnostics}}
#'
#' @examples
#' plotDetectionRates(quickstartresults)
"quickstartresults"
