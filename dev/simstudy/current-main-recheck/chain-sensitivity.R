#!/usr/bin/env Rscript
# Numerical sensitivity for every selected fit that still fails the screen.
# These ranges are diagnostic choices among chain means, not uncertainty CIs.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
repo <- normalizePath(arg("repo", ".")); study <- normalizePath(arg("study"))
out <- file.path(repo, "dev/simstudy/current-main-recheck/results")
selection <- readRDS(file.path(study, "selection.rds"))
diagnostics <- read.csv(file.path(out, "diagnostic-summary.csv"))
verification <- read.csv(file.path(out, "verification-selected.csv"))
flagged <- with(diagnostics, warnings > 0 | max_group_rhat > 1.05 |
  max_element_rhat > 1.05 | unresolved_rhat > 0)
stopifnot(!anyNA(flagged), nrow(selection) == 110L,
  setequal(selection$key, diagnostics$key), setequal(selection$key, verification$key))
keys <- diagnostics$key[flagged]
rows <- list()
for (key in keys) {
  s <- selection[selection$key == key, ]; v <- verification[verification$key == key, ]
  fit_file <- sub("-result\\.rds$", "-fit.rds", s$result_file)
  stopifnot(unname(tools::md5sum(s$result_file)) == v$result_md5,
    unname(tools::md5sum(fit_file)) == v$fit_md5)
  r <- readRDS(s$result_file); saved <- readRDS(fit_file); fit <- saved$fit
  j <- fit$results_output$jsdm_output
  ni <- dim(j$B0_output)[2]; nc <- dim(j$B0_output)[3]; ns <- dim(j$B0_output)[1]
  stopifnot(ni == r$mcmc$niter, nc == r$mcmc$nchain)
  estimates <- array(0, c(100L, ns, nc))
  for (ch in seq_len(nc)) for (sp in seq_len(ns)) {
    beta <- matrix(j$B_output[, sp, , ch], nrow = ncol(fit$X_psi), ncol = ni)
    eta <- sweep(fit$X_psi[1:100, , drop = FALSE] %*% beta, 2, j$B0_output[sp, , ch], "+")
    for (k in seq_len(dim(j$U_output)[2]))
      eta <- eta + sweep(matrix(j$U_output[1:100, k, , ch], 100, ni), 2, j$L_output[k, sp, , ch], "*")
    estimates[, sp, ch] <- rowMeans(plogis(eta))
  }
  pooled <- apply(estimates, c(1, 2), mean)
  stopifnot(max(abs(pooled - r$estimate[1:100, , drop = FALSE])) < 1e-12)
  truth <- r$truth[1:100, , drop = FALSE]
  for (ch in 0:nc) for (group in c("all", "low", "medium", "high")) {
    e <- if (ch == 0L) pooled else estimates[, , ch]
    mask <- switch(group, all = matrix(TRUE, 100, ns), low = truth < .2,
      medium = truth >= .2 & truth <= .8, high = truth > .8)
    stopifnot(any(mask))
    rows[[paste(key, ch, group)]] <- data.frame(s[c("key", "family", "scenario", "arm", "replicate")],
      chain = ch, group = group, bias = mean(e[mask] - truth[mask]), mae = mean(abs(e[mask] - truth[mask])))
  }
  cat(key, "chain means reconstructed and pooled mean verified.\n"); flush.console()
  rm(saved, fit, j, estimates, eta); gc(FALSE)
}
chain_scores <- if (length(rows)) do.call(rbind, rows) else data.frame(key = character(),
  family = character(), scenario = character(), arm = character(), replicate = integer(),
  chain = integer(), group = character(), bias = numeric(), mae = numeric())
write.csv(chain_scores, file.path(out, "flagged-chain-scores.csv"), row.names = FALSE)

scores <- read.csv(file.path(out, "community-scores.csv"))
p <- subset(scores, version == "current" & metric == "occupancy_original_sites")
stopifnot(nrow(p) == 440L, !anyDuplicated(p[c("key", "group")]))
p$min_mae <- p$max_mae <- p$mae
p$min_bias <- p$max_bias <- p$bias
p$flagged <- p$key %in% keys
for (i in which(p$flagged)) {
  z <- chain_scores[chain_scores$key == p$key[i] & chain_scores$group == p$group[i], ]
  stopifnot(sum(z$chain == 0L) == 1L, abs(z$mae[z$chain == 0L] - p$mae[i]) < 1e-12,
    abs(z$bias[z$chain == 0L] - p$bias[i]) < 1e-12)
  p$min_mae[i] <- min(z$mae); p$max_mae[i] <- max(z$mae)
  p$min_bias[i] <- min(z$bias); p$max_bias[i] <- max(z$bias)
}
ranges <- do.call(rbind, lapply(split(p, interaction(p$family, p$scenario, p$arm, p$group, drop = TRUE)), function(d) {
  stopifnot(nrow(d) == 10L, setequal(d$replicate, 1:10))
  data.frame(d[1, c("family", "scenario", "arm", "group")], communities = 10L,
    flagged_fits = sum(d$flagged), selected_mae = mean(d$mae),
    minimum_mae = mean(d$min_mae), maximum_mae = mean(d$max_mae),
    selected_bias = mean(d$bias), minimum_bias = mean(d$min_bias), maximum_bias = mean(d$max_bias))
}))
write.csv(ranges, file.path(out, "chain-sensitivity-summary.csv"), row.names = FALSE)
contrasts <- subset(read.csv(file.path(out, "design-paired-summary.csv")), version == "current")
contrasts$minimum_reduction <- contrasts$maximum_reduction <- NA_real_
for (i in seq_len(nrow(contrasts))) {
  contrast <- contrasts[i, ]
  ref <- subset(ranges, family == contrast$family & scenario == contrast$scenario & group == contrast$group & arm == contrast$reference)
  alt <- subset(ranges, family == contrast$family & scenario == contrast$scenario & group == contrast$group & arm == contrast$alternative)
  stopifnot(nrow(ref) == 1L, nrow(alt) == 1L,
    abs((ref$selected_mae - alt$selected_mae) - contrast$reduction_mae) < 1e-12)
  contrasts$minimum_reduction[i] <- ref$minimum_mae - alt$maximum_mae
  contrasts$maximum_reduction[i] <- ref$maximum_mae - alt$minimum_mae
}
write.csv(contrasts, file.path(out, "design-chain-sensitivity.csv"), row.names = FALSE)
cat(length(keys), "flagged fits retained; all", nrow(p), "primary scores remain in the ten-community summaries.\n")
print(subset(ranges, flagged_fits > 0 & group == "all"), row.names = FALSE)
