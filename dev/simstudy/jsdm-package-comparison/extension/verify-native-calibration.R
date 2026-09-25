# Independent verification; does not source the worker or interval helpers.
# Rscript verify-native-calibration.R CALIBRATION_BUNDLE
args <- commandArgs(trailingOnly = TRUE)
bundle <- readRDS(args[1]); p <- bundle$native_provenance
stopifnot(!is.null(p), identical(tools::md5sum(names(p$base)), p$base),
  identical(tools::md5sum(names(p$jobs)), p$jobs),
  identical(tools::md5sum(names(p$worker_dependencies)), p$worker_dependencies),
  identical(tools::md5sum(names(p$sources)), p$sources))
base <- readRDS(names(p$base))
stopifnot(identical(bundle$prediction, base$prediction), identical(bundle$manifest, base$manifest),
          identical(bundle$integration, base$integration))
unchanged <- setdiff(names(base$rows), c("lower", "upper"))
stopifnot(identical(bundle$rows[unchanged], base$rows[unchanged]))
eligible <- with(bundle$manifest, scored & package %in% c("gllvm", "sjSDM") & scenario != "traits" &
                   !(scenario == "curved" & response == "linear"))
stopifnot(setequal(bundle$native_intervals$job, bundle$manifest$job[eligible]),
          !anyDuplicated(bundle$native_intervals$job))
changed <- bundle$rows$job %in% bundle$manifest$job[eligible] & bundle$rows$target == "Environmental coefficient"
stopifnot(identical(bundle$rows[!changed, c("lower", "upper")], base$rows[!changed, c("lower", "upper")]))
for (path in names(p$jobs)) {
  r <- readRDS(path); j <- r$job
  stopifnot(!r$truth_used, identical(tools::md5sum(names(r$input_hashes)), r$input_hashes),
            identical(tools::md5sum(names(r$source_hashes)), r$source_hashes))
  d <- subset(bundle$rows, job == j$job & target == "Environmental coefficient")
  status <- subset(bundle$native_intervals, job == j$job)
  stopifnot(identical(status$ok, r$ok), identical(status$numerical_pass, r$numerical_pass))
  if (!r$ok) {
    stopifnot(all(is.na(d$lower)), all(is.na(d$upper)), !status$numerical_pass)
    next
  }
  input <- readRDS(names(r$input_hashes)[2])
  fit <- readRDS(names(r$input_hashes)[1])
  k <- nrow(fit$parameters$beta); species <- ncol(fit$parameters$beta)
  A <- diag(k); A[-1,-1] <- diag(1/input$spread, k-1L)
  A[1,-1] <- -input$centre/input$spread
  for (s in seq_len(species)) {
    V <- r$audit$covariance[[s]]; idx <- (s-1L)*k+seq_len(k)
    valid <- all(is.finite(V)) && max(abs(V-t(V))) <= 1e-8*max(1,max(abs(V))) &&
      min(eigen(V,symmetric=TRUE,only.values=TRUE)$values) > 0
    if (!valid) {
      stopifnot(all(is.na(d$lower[idx])), all(is.na(d$upper[idx])))
      next
    }
    C <- A %*% V %*% t(A)
    m <- A %*% fit$parameters$beta[,s]
    se <- sqrt(diag(C))
    stopifnot(max(abs(d$estimate[idx]-m)) < 1e-10,
      max(abs(d$lower[idx]-(m-qnorm(.975)*se))) < 1e-10,
      max(abs(d$upper[idx]-(m+qnorm(.975)*se))) < 1e-10,
      max(abs(diag(V)-r$audit$native_se[,s]^2)) < 1e-8)
    if (j$package == "sjSDM") {
      H <- r$audit$regularized_hessian[s,,]
      stopifnot(max(abs(H %*% V-diag(k))) < 1e-8)
    }
  }
  stopifnot(r$audit$parameter_change <= 1e-8)
  if (j$package == "gllvm") stopifnot(r$audit$loglik_change <= 1e-8)
  if (j$package == "sjSDM") {
    stopifnot(identical(tools::md5sum(names(r$audit$standard_error_source)), r$audit$standard_error_source))
    audited_se <- function(blocks) unlist(lapply(blocks, function(v) {
      if (any(!is.finite(v)) || max(abs(v-t(v))) > 1e-8*max(1,max(abs(v))) ||
          min(eigen(v,symmetric=TRUE,only.values=TRUE)$values) <= 0) return(rep(NA_real_, k))
      raw_variance <- diag(A %*% v %*% t(A))
      if (any(raw_variance <= 0)) return(rep(NA_real_, k))
      sqrt(raw_variance)
    }))
    primary <- audited_se(r$audit$covariance)
    secondary <- audited_se(r$audit$secondary_covariance)
    change <- if (anyNA(c(primary, secondary))) Inf else max(abs(primary-secondary)/pmax(primary,secondary))
    stopifnot(if (is.infinite(change)) is.infinite(status$se_change) else abs(change-status$se_change) < 1e-10,
      !status$numerical_pass || change <= .05)
  }
}
summary_keys <- c("package", "scenario", "n_sites", "n_species", "response", "use_traits", "term")
key <- function(d) do.call(paste, c(d[summary_keys], sep="|"))
good <- with(bundle$native_intervals, ok & numerical_pass & original_diagnostic_pass)
good_rows <- subset(bundle$rows, job %in% bundle$native_intervals$job[good] & target == "Environmental coefficient")
stopifnot(setequal(unique(key(good_rows)), key(bundle$native_passed_summary)),
          !anyDuplicated(key(bundle$native_passed_summary)))
for (i in seq_len(nrow(bundle$native_passed_summary))) {
  s <- bundle$native_passed_summary[i, ]
  d <- good_rows
  for (column in summary_keys) d <- d[d[[column]] == s[[column]], ]
  stopifnot(length(unique(d$replicate)) == s$communities,
    abs(mean(tapply(d$truth >= d$lower & d$truth <= d$upper, d$replicate, mean))-s$coverage) < 1e-12,
    abs(mean(tapply(d$upper-d$lower, d$replicate, mean))-s$width) < 1e-12)
}
native_communities <- subset(bundle$per_community, package %in% c("gllvm", "sjSDM") &
  target == "Environmental coefficient" & scenario != "traits")
native_rows <- subset(bundle$rows, changed)
for (i in seq_len(nrow(native_communities))) {
  s <- native_communities[i, ]; d <- native_rows[native_rows$replicate == s$replicate, ]
  for (column in summary_keys) d <- d[d[[column]] == s[[column]], ]
  valid <- is.finite(d$lower) & is.finite(d$upper) & d$lower <= d$upper
  stopifnot(s$elements == nrow(d), s$intervals == sum(valid))
  if (any(valid)) stopifnot(abs(s$coverage - mean(with(d[valid, ], truth >= lower & truth <= upper))) < 1e-12,
                           abs(s$width - mean(with(d[valid, ], upper-lower))) < 1e-12)
  else stopifnot(is.na(s$coverage), is.na(s$width))
}
cat("Verified native intervals, native covariance identities, archived fit agreement and unchanged point/other interval results for", length(p$jobs), "jobs.\n")
