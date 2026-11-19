#!/usr/bin/env Rscript
# Draws figures/four-chains-on-the-gradients.svg from one fit of the rotation test: the low-rank
# truth with signal in 10 columns, seed 3, columns as simulated, 4 chains of 5,000 + 5,000 with
# attempt 3's library. Each chain's three gradients are the leading directions of its posterior
# mean of the embedding's contribution; a species' position is its response to each. Every chain's
# three gradients are turned, as a rigid set, to lie as close as possible to the truth's, which
# changes nothing about the space they span. Reads a per-run _summary.rds (kept locally, not
# tracked), so it runs only where the fit was made. Style follows the two figures in
# dev/superpowers/specs/2026-11-11-grouped-covariate-prior-figures.
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/figures/make_four_chains_figure.R
root <- "dev/simstudy/grouped-mixing-20261118"
x <- readRDS(file.path(root, "results/rotation/lowrank_s3_rot0_b5000_summary.rds"))
r <- 3L
positions <- function(M) { s <- svd(scale(M, scale = FALSE), nu = 0, nv = r); s$v %*% diag(s$d[1:r]) }
truth <- positions(x$true_contribution)
turn <- function(Z) { s <- svd(t(Z) %*% truth); Z %*% (s$u %*% t(s$v)) }      # rigid turn onto the truth
chains <- lapply(seq_len(dim(x$contribution_mean)[3]), function(k) turn(positions(x$contribution_mean[, , k])))

panel <- function(x0, y0, w, h, a, b) {
  all <- do.call(rbind, c(list(truth), chains))
  lim <- max(abs(all[, c(a, b)])) * 1.08
  px <- function(v) x0 + w / 2 + v / lim * (min(w, h) / 2 - 8)
  py <- function(v) y0 + h / 2 - v / lim * (min(w, h) / 2 - 8)
  out <- c(sprintf('<rect x="%g" y="%g" width="%g" height="%g" rx="4" fill="none" stroke="#c3c2b7" stroke-width="0.75"/>', x0, y0, w, h),
           sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#c3c2b7" stroke-width="0.75"/>', x0 + 8, py(0), x0 + w - 8, py(0)),
           sprintf('<line x1="%.1f" y1="%g" x2="%.1f" y2="%g" stroke="#c3c2b7" stroke-width="0.75"/>', px(0), y0 + 8, px(0), y0 + h - 8),
           sprintf('<text class="ts" x="%g" y="%g">↑ gradient %d</text>', x0 + 8, y0 + 18, b),
           sprintf('<text class="ts" x="%g" y="%g">gradient %d →</text>', x0 + 8, y0 + h - 8, a),
           sprintf('<circle cx="%.1f" cy="%.1f" r="5.5" fill="none" stroke="#888780" stroke-width="1.2"/>', px(truth[, a]), py(truth[, b])))
  # a faint tie from each ring to each of its four dots, so that a species' dots can be told apart from its neighbours'
  for (Z in chains) out <- c(out, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888780" stroke-opacity="0.45" stroke-width="0.6"/>', px(truth[, a]), py(truth[, b]), px(Z[, a]), py(Z[, b])))
  for (Z in chains) out <- c(out, sprintf('<circle cx="%.1f" cy="%.1f" r="2.4" fill="#1D9E75" fill-opacity="0.8"/>', px(Z[, a]), py(Z[, b])))
  out
}
svg <- c('<?xml version="1.0" encoding="UTF-8"?>',
  '<svg xmlns="http://www.w3.org/2000/svg" width="680" height="452" viewBox="0 0 680 452" role="img" font-family="Helvetica, Arial, sans-serif">',
  '<title>Where four chains put each species on the gradients</title>',
  '<desc>Two panels from one fit of simulated data with three gradients. Each of 48 species is drawn as a ring at its true position on the gradients, with four dots, one per chain, at the positions the four chains give it. The left panel shows gradients 1 and 2 and the right panel gradients 1 and 3. The four dots of a species fall close together and close to its ring.</desc>',
  '<defs><style>.th{font-size:14px;font-weight:bold;fill:#1a1a19}.ts{font-size:12px;fill:#52514e}</style></defs>',
  '<rect x="0" y="0" width="680" height="452" fill="#ffffff"/>',
  '<text class="th" x="40" y="36">One fit, four chains: where each species sits on the gradients</text>',
  '<text class="ts" x="40" y="54">The tips of the species\' arrows, seen within the sheet. Simulated data with three gradients.</text>',
  panel(40, 68, 290, 272, 1, 2), panel(350, 68, 290, 272, 1, 3),
  '<circle cx="48" cy="362" r="5.5" fill="none" stroke="#888780" stroke-width="1.2"/>',
  '<text class="ts" x="62" y="366">A species\' true position</text>',
  '<circle cx="232" cy="362" r="2.4" fill="#1D9E75" fill-opacity="0.8"/><circle cx="240" cy="360" r="2.4" fill="#1D9E75" fill-opacity="0.8"/><circle cx="236" cy="366" r="2.4" fill="#1D9E75" fill-opacity="0.8"/><circle cx="244" cy="365" r="2.4" fill="#1D9E75" fill-opacity="0.8"/>',
  '<text class="ts" x="256" y="366">Where each of the four chains puts it, tied to its ring</text>',
  '<text class="ts" x="40" y="398">The four chains agree on the set of gradients: 0.96 to 0.99, where 1 is the same set.</text>',
  '<text class="ts" x="40" y="416">Their gradients hold 95 to 98 percent of the true signal.</text>',
  '<text class="ts" x="40" y="434">In this same fit, 79 percent of single site-by-species values fail the R-hat check between chains.</text>',
  '</svg>')
writeLines(svg, file.path(root, "figures/four-chains-on-the-gradients.svg"))
spread <- sapply(1:r, function(j) mean(apply(sapply(chains, function(Z) Z[, j]), 1, stats::sd)) / stats::sd(truth[, j]))
cat("mean spread of a species across the four chains, as a share of the spread between species, by gradient:", round(spread, 3), "\n")
