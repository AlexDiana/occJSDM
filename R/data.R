#' Simulated example dataset for runOccPlus()
#'
#' A simulated eDNA metabarcoding dataset, in the \code{info}/\code{OTU} list
#' format expected by \code{\link{runOccPlus}}. Generated via
#' \code{\link{simulateOccPlusData}} and converted with
#' \code{\link{toRunOccPlusFormat}}.
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
#'         across each sample's replicates). \code{runOccPlus()} adds its
#'         own intercept internally, so no intercept column is included
#'         here.}
#'       \item{Xs.*}{Site coordinates (site-level, repeated across each
#'         site's replicates), for use as \code{spatCovariates} in
#'         \code{runOccPlus()}.}
#'     }
#'   }
#'   \item{OTU}{A numeric matrix of read counts, with \code{nrow(OTU) ==
#'     nrow(info)} and one column per species (named \code{OTU_1},
#'     \code{OTU_2}, ...).}
#'   \item{spatCovariates}{A character vector of the \code{Xs.*} column
#'     names in \code{info} (\code{c("Xs.1", "Xs.2")}), ready to pass as the
#'     \code{spatCovariates} argument of \code{runOccPlus()}.}
#'   \item{traitsMatrix}{A species (rows) by trait (columns) matrix, with
#'     row names matching \code{colnames(OTU)} and column names
#'     \code{Trait_1}, \code{Trait_2}, ..., ready to pass as the
#'     \code{traitsMatrix} argument of \code{runOccPlus()}.}
#' }
#'
#' @seealso \code{\link{simulateOccPlusData}}, \code{\link{toRunOccPlusFormat}},
#'   \code{\link{runOccPlus}}
#'
#' @examples
#' str(sampledata, max.level = 1)
#' head(sampledata$info)
"sampledata"
