# Independent checks of the compact teaching evidence against saved scoring runs.
# Rscript dev/simstudy/vignette-lesson/verify-waic.R ARCHIVE
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
b <- readRDS('vignettes/teaching-data/site-waic-lesson.rds')
stopifnot(identical(b$source_md5, tools::md5sum(names(b$source_md5))),
          identical(b$exporter_md5, tools::md5sum(names(b$exporter_md5))),
          identical(b$validation_md5, tools::md5sum(names(b$validation_md5))))
for (n in c(1000L, 4000L)) {
  x <- readRDS(file.path(archive,paste0('waic-',n),'validation.rds'))
  stopifnot(identical(x$source_md5,b$source_md5))
  for (model in c('one_factor','two_factors')) {
    fit_path <- if (model == 'one_factor') file.path(archive,'prediction','one-factor-fit.rds') else file.path(archive,'default-fit.rds')
    stopifnot(identical(x$results[[model]]$source_fit_md5, unname(tools::md5sum(fit_path))))
    score <- x$results[[model]]$score
    ll <- score$log_lik
    row <- b$summaries[b$summaries$model == model & b$summaries$draws == n, ]
    stopifnot(nrow(row) == 1L, nrow(ll) == n, ncol(ll) == 100L,
              identical(as.integer(table(score$draw_ids$chain)), rep(n %/% 4L,4L)),
              all(score$observations$threshold == 1))
    penalty <- apply(ll,2,var)
    logmean <- apply(ll,2,function(z) {m <- max(z); m + log(mean(exp(z-m)))})
    expected <- -2*(logmean-penalty)
    stopifnot(abs(row$WAIC-sum(expected)) < 1e-8,
              abs(row$p_waic-sum(penalty)) < 1e-8,
              row$flagged_sites == sum(penalty > 0.4),
              identical(row$fit_md5,x$results[[model]]$source_fit_md5))
    if(n == 4000L) stopifnot(identical(b$pointwise[[model]],score$pointwise))
  }
  difference <- x$results$one_factor$score$pointwise$waic - x$results$two_factors$score$pointwise$waic
  row <- b$comparisons[b$comparisons$draws == n, ]
  stopifnot(abs(row$difference_one_minus_two-sum(difference)) < 1e-8,
            abs(row$paired_site_SE-sqrt(length(difference)*var(difference))) < 1e-8,
            abs(row$approximate_difference_MCSE-sqrt(sum(vapply(x$results,function(z)z$waic_mcse^2,numeric(1))))) < 1e-8)
}
cat('Current WAIC bundle: source, fits, draws, pointwise penalties, score differences and uncertainty verified.\n')
