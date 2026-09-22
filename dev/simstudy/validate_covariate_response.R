# Validate the corrected output using an unchanged archived Lesson 3 fit.
# No MCMC. Use an installation of the branch under review via R_LIBS.
# Rscript dev/simstudy/validate_covariate_response.R /path/to/default-fit.rds
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
library(occJSDM)

saved <- readRDS(args[1])
fit <- saved$fit
covariate <- "X_psi.EnvCov.1"
species <- c(1L, 10L)
confidence <- .8
meta <- fit$infos$list_X_psi_mat
stopifnot(all(meta$is_numeric), all(vapply(meta$bs_info, is.null, logical(1))))

# Recover original input values independently of the output helper's X0 field.
info <- fit$infos$data_info
raw <- info[!duplicated(info$Site), meta$names_df, drop = FALSE]
x <- seq(min(raw[[covariate]]), max(raw[[covariate]]), length.out = 200L)
reference <- vapply(raw, median, numeric(1))
design <- matrix(reference, nrow = 200L, ncol = length(reference), byrow = TRUE)
colnames(design) <- names(reference)
design[, covariate] <- x
for (name in meta$names_df) {
  design[, name] <- (design[, name] - meta$mean_df[[name]]) / meta$sd_df[[name]]
}
design <- design[, colnames(fit$X_psi), drop = FALSE]

result <- returnCovariateEffect(fit, covariate, species, confidence)
stopifnot(all(result$lower >= 0), all(result$upper <= 1))
posterior <- fit$results_output$jsdm_output
probes <- c(1L, 100L, 200L)
for (sp in species) {
  # Explicit chain-by-chain extraction avoids the helper's array reshaping.
  coefficients <- do.call(rbind, lapply(seq_len(dim(posterior$B_output)[4]), function(ch) {
    t(posterior$B_output[, sp, , ch])
  }))
  intercepts <- unlist(lapply(seq_len(dim(posterior$B0_output)[3]), function(ch) {
    posterior$B0_output[sp, , ch]
  }))
  response <- plogis(sweep(design[probes, , drop = FALSE] %*% t(coefficients),
                            2L, intercepts, "+"))
  expected <- t(apply(response, 1L, quantile, probs = c(.5, .1, .9), names = FALSE))
  rows <- result[result$Species == fit$infos$speciesNames[sp], ]
  stopifnot(isTRUE(all.equal(rows$x, x, tolerance = 1e-12)),
            isTRUE(all.equal(unname(as.matrix(rows[probes, c("median", "lower", "upper")])),
                             expected, tolerance = 1e-12)))
}
plot <- plotCovariateEffect(fit, covariate, species, confidence)[[1]]
stopifnot(identical(plot$data, result), identical(plot$labels$y, "Occupancy probability"))
invisible(ggplot2::ggplot_build(plot))
print(dplyr::summarise(dplyr::group_by(result, Species),
                       minimum_median = min(median), maximum_median = max(median),
                       invalid_medians = sum(median < 0 | median > 1)))
cat("Original predictor range:", range(result$x), "\n")
cat("All checks passed; saved-fit MD5:", unname(tools::md5sum(args[1])), "\n")
