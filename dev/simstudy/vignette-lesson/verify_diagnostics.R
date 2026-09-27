# Independent checks against the original archived arrays and simulator inputs.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
source("dev/simstudy/vignette-lesson/helpers.R")
x <- readRDS("vignettes/teaching-data/diagnostics-lesson.rds")
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
stopifnot(identical(x$source_hashes, lesson_source_hashes()),
          identical(x$lesson_md5, unname(tools::md5sum("vignettes/teaching-data/nonspatial-lesson.rds"))),
          identical(x$exporter_md5, unname(tools::md5sum("dev/simstudy/vignette-lesson/summarise_diagnostics.R"))))

for (arm in names(x$fit_manifests)) {
  m <- x$fit_manifests[[arm]]
  path <- file.path(archive, m$file)
  stopifnot(identical(unname(tools::md5sum(path)), m$md5))
  fit <- readRDS(path)$fit
  for (name in names(x$traces)) {
    a <- x$traces[[name]]
    if (a$arm != arm) next
    dims <- dim(a$draws)
    stopifnot(identical(tail(dims, 2), c(m$mcmc$niter, m$mcmc$nchain)),
              m$mcmc$nthin == 1, all(is.finite(a$draws)))
    for (i in seq_len(nrow(a$truth))) {
      species <- if (a$param == "theta0") a$truth$label1[i] else a$truth$label2[i]
      s <- match(species, fit$infos$speciesNames)
      if (a$param == "beta_theta") {
        k <- match(a$truth$label1[i], colnames(fit$X_theta))
        raw <- fit$results_output$beta_theta_output[k, s, , ]
        recorded <- a$draws[1, i, , ]
        # Separate scale check using the actual raw field-sample covariates.
        info <- lesson$input$sim$data_list$info
        samples <- unique(info[c("Site", "Sample", "X_theta")])
        truth <- lesson$input$sim$true_params$beta_theta_true[2, s] * sd(samples$X_theta)
      } else if (a$param == "p") {
        k <- match(a$truth$label1[i], as.character(fit$infos$primerNames))
        raw <- fit$results_output$p_output[k, s, , ]
        recorded <- a$draws[1, i, , ]
        p <- lesson$input$params
        truth <- p$p[k, s] * pnorm(log(1.5), p$mu1, p$sigma1, lower.tail = FALSE)
      } else {
        raw <- fit$results_output$theta0_output[s, , ]
        recorded <- a$draws[i, , ]
        truth <- lesson$input$params$theta0[s]
      }
      stopifnot(identical(recorded, raw), abs(a$truth$truth[i] - truth) < 1e-12)
      # Check the public plotting function preserves every draw and chain ID.
      plot <- occJSDM::plotTraceplot(a$draws, dimnames1 = a$label1, dimnames2 = a$label2)
      table <- plot$data
      for (ch in seq_len(ncol(raw))) {
        keep <- table$label1 == a$truth$label1[i] & table$label2 == a$truth$label2[i] &
          as.character(table$chain) == as.character(ch)
        stopifnot(identical(table$iter[keep], seq_len(nrow(raw))),
                  identical(table$value[keep], unname(raw[, ch])))
      }
    }
    cat(name, ": every exported draw, chain, iteration and truth line verified\n")
  }
}

flags <- lesson$diagnostics$alternative
field <- flags[flags$param == "theta0", ]
stopifnot(x$traces$field_contamination$label1 == field$label1[which.min(field$ess)])
stopifnot(sum(flags$rhat > 1.01 | flags$ess < 400) == 3L)
cat("Diagnostic selection, fit hashes, source hashes and plotting identities verified.\n")
