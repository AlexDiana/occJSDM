#!/usr/bin/env Rscript
# Baseline against prototype for one fix attempt: the mixing table and the check that the answer
# has not moved. Reads what repro.R wrote under results/<attempt>/base and results/<attempt>/proto.
#
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/compare.R [--attempt=attempt1] [--nburn=1500]
#
# Mixing columns are those of repro.R. Tau is kept draws per effective draw (coda::effectiveSize,
# which fits an autoregression within each chain, so it cannot see more persistence than a chain is
# long: read a tau above about a tenth of the chain length as "at least").
#
# Agreement of the answer. Each arm has two chains. For a quantity with chain means c1 and c2, the
# arm mean is (c1 + c2) / 2 and (c1 - c2)^2 / 4 estimates its variance (one degree of freedom, so
# it is only used summed over many quantities). If the two arms sample the same posterior, the mean
# square difference between the arm means equals the sum of the two arms' variances, so
#   ratio = sqrt(mean((proto - base)^2) / mean(var_base + var_proto))
# is about 1. A ratio well above 1 means the arms disagree by more than their chains disagree with
# each other. It is reported per block of quantities, with the largest standardised difference.
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
attempt <- getopt("attempt", "attempt1")
nburn <- as.integer(getopt("nburn", "1500"))
root <- file.path("dev/simstudy/grouped-mixing-20261118/results", attempt)
read_arm <- function(arm) {
  f <- list.files(file.path(root, arm), pattern = sprintf("_b%d_s[0-9]+\\.csv$", nburn), full.names = TRUE)
  d <- do.call(rbind, lapply(f, utils::read.csv))
  d$arm <- arm
  d
}
base <- read_arm("base"); proto <- read_arm("proto")
tags <- intersect(base$tag, proto$tag)
stopifnot(length(tags) > 0)
base <- base[match(tags, base$tag), ]; proto <- proto[match(tags, proto$tag), ]
cols <- c(loglam_tau_med_inactive = "tau log axis scale, inactive (median)",
          loglam_tau_med_active = "tau log axis scale, active (median)",
          loglam_rhat_gt1.1 = "log axis scales over R-hat 1.1 (of 64)",
          loglam_rhat_max = "largest log axis scale R-hat",
          logsig_rhat = "log sigma_emb R-hat", logsig_tau = "log sigma_emb tau",
          B_emb_share_rhat_gt1.1 = "share of embedding B over R-hat 1.1",
          B_emb_tau_med = "tau embedding B (median)", seconds = "seconds per fit")
cat(sprintf("Attempt %s, %d + %d, baseline %s, prototype %s\n\n", attempt, nburn, base$niter[1],
            base$package_commit[1], proto$package_commit[1]))
mixing <- do.call(rbind, lapply(seq_along(tags), function(i) data.frame(
  run = sprintf("%s seed %d", base$truth[i], base$seed[i]), quantity = unname(cols),
  baseline = signif(unlist(base[i, names(cols)]), 3), prototype = signif(unlist(proto[i, names(cols)]), 3),
  ratio = signif(unlist(base[i, names(cols)]) / unlist(proto[i, names(cols)]), 3), row.names = NULL)))
print(mixing, row.names = FALSE)

agree <- function(b1, b2, p1, p2) {
  diff <- (p1 + p2) / 2 - (b1 + b2) / 2
  v <- (b1 - b2)^2 / 4 + (p1 - p2)^2 / 4
  c(n = length(diff), rms_diff = sqrt(mean(diff^2)), rms_expected = sqrt(mean(v)),
    ratio = sqrt(mean(diff^2) / mean(v)), largest_abs_diff = max(abs(diff)))
}
agreement <- do.call(rbind, lapply(tags, function(tag) {
  b <- readRDS(file.path(root, "base", paste0(tag, "_summary.rds")))
  p <- readRDS(file.path(root, "proto", paste0(tag, "_summary.rds")))
  emb <- b$emb_rows; meas <- setdiff(seq_len(nrow(b$B_mean)), emb)
  blocks <- list(
    "B, measured rows" = agree(b$B_chain_mean[meas, , 1], b$B_chain_mean[meas, , 2],
                               p$B_chain_mean[meas, , 1], p$B_chain_mean[meas, , 2]),
    "B, embedding rows" = agree(b$B_chain_mean[emb, , 1], b$B_chain_mean[emb, , 2],
                                p$B_chain_mean[emb, , 1], p$B_chain_mean[emb, , 2]),
    "embedding spread per species" = agree(b$spread_chain_mean[, 1], b$spread_chain_mean[, 2],
                                           p$spread_chain_mean[, 1], p$spread_chain_mean[, 2]),
    "log axis scales" = agree(b$log_lambda_chain_mean[, 1], b$log_lambda_chain_mean[, 2],
                              p$log_lambda_chain_mean[, 1], p$log_lambda_chain_mean[, 2]),
    "log sigma_emb" = agree(b$log_sigma_emb_chain_mean[1], b$log_sigma_emb_chain_mean[2],
                            p$log_sigma_emb_chain_mean[1], p$log_sigma_emb_chain_mean[2]),
    "log sigma_b" = agree(b$log_sigma_b_chain_mean[1], b$log_sigma_b_chain_mean[2],
                          p$log_sigma_b_chain_mean[1], p$log_sigma_b_chain_mean[2]))
  data.frame(run = sub("_n.*_s", " seed ", tag), block = names(blocks),
             signif(do.call(rbind, blocks), 3), row.names = NULL)
}))
cat("\nAgreement of the answer between the arms (ratio about 1 = same posterior):\n")
print(agreement, row.names = FALSE)

cat("\nChain means of log sigma_emb (baseline chains 1 and 2, prototype chains 1 and 2), and posterior means of a few quantities:\n")
for (tag in tags) {
  b <- readRDS(file.path(root, "base", paste0(tag, "_summary.rds")))
  p <- readRDS(file.path(root, "proto", paste0(tag, "_summary.rds")))
  cat(sprintf("%s\n  log sigma_emb chain means: base %s | proto %s\n", tag,
              paste(round(b$log_sigma_emb_chain_mean, 2), collapse = ", "),
              paste(round(p$log_sigma_emb_chain_mean, 2), collapse = ", ")))
  cat(sprintf("  mean log axis scale over columns: base %.2f, proto %.2f; mean embedding spread: base %.3f, proto %.3f\n",
              mean(b$log_lambda_mean), mean(p$log_lambda_mean), mean(b$spread_mean), mean(p$spread_mean)))
  if (!is.null(p$axis_move))
    cat(sprintf("  axis move: acceptance after burn-in median %.2f (range %.2f to %.2f); step median %.2f (range %.2f to %.2f)\n",
                stats::median(p$axis_move$acceptance), min(p$axis_move$acceptance), max(p$axis_move$acceptance),
                stats::median(p$axis_move$step), min(p$axis_move$step), max(p$axis_move$step)))
}
utils::write.csv(mixing, file.path(root, sprintf("mixing_b%d.csv", nburn)), row.names = FALSE)
utils::write.csv(agreement, file.path(root, sprintf("agreement_b%d.csv", nburn)), row.names = FALSE)
