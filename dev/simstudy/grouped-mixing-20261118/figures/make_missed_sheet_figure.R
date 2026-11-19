#!/usr/bin/env Rscript
# Draws figures/chains-agree-but-miss-the-sheet.svg: the other failure of the grouped prior, in
# which the leftovers carry the species' responses and no sheet is found. The top row is a drawing
# in the style of dev/superpowers/specs/2026-11-11-grouped-covariate-prior-figures. The bottom row
# is from one fit: the low-rank truth with signal in 10 columns, seed 2, attempt 1's library
# (eba6a20), 2 chains of 1,500 + 1,500, whose first chain ended with the blends carrying the
# signal and whose second with the leftovers carrying it. It plots each chain's posterior mean of
# the coefficients of the 10 active columns against the truth. Reads a per-run _summary.rds (kept
# locally, not tracked), so it runs only where the fit was made.
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/figures/make_missed_sheet_figure.R
root <- "dev/simstudy/grouped-mixing-20261118"
x <- readRDS(file.path(root, "results/attempt1/proto/lowrank_n256_S48_K64_gt4_b1500_s2_summary.rds"))
emb <- x$emb_rows; active <- emb[rowSums(x$true_B[emb, ]^2) > 0]
truth <- as.vector(x$true_B[active, ])
est <- lapply(1:2, function(k) as.vector(x$B_chain_mean[active, , k]))
slope <- sapply(est, function(e) sum(e * truth) / sum(truth^2))
stopifnot(x$log_sigma_emb_chain_mean[1] < 0.5, x$log_sigma_emb_chain_mean[2] > 0.5)   # chain 1 blends, chain 2 leftovers

# ---- top right: the species' arrows of the earlier figure, shrunk, with their true lengths behind
tips <- rbind(c(296, 235.3), c(275, 204.4), c(178, 200.3), c(128, 230.7), c(111, 275), c(151, 299.2),
              c(222, 299.7), c(279, 279.6), c(234.5, 231.5), c(330.5, 182.8))
o <- c(200, 250); shrink <- 0.42
set.seed(4); wobble <- matrix(rnorm(20, sd = 5), 10, 2)          # each species' estimate is its own
short <- sweep(sweep(tips, 2, o) * shrink + wobble, 2, o, "+")
leftover_panel <- c(
  '<line x1="200" y1="250" x2="400" y2="250" stroke="#c3c2b7" stroke-width="1.2"/>',
  '<line x1="200" y1="250" x2="200" y2="90" stroke="#c3c2b7" stroke-width="1.2"/>',
  '<line x1="200" y1="250" x2="120" y2="330" stroke="#c3c2b7" stroke-width="1.2"/>',
  '<polygon points="352,196.4 256,332.4 48,303.6 144,167.6" fill="none" stroke="#7F77DD" stroke-opacity="0.35" stroke-width="0.9" stroke-dasharray="4 4"/>',
  sprintf('<line x1="200" y1="250" x2="%.1f" y2="%.1f" stroke="#c3c2b7" stroke-width="1.6"/>', tips[, 1], tips[, 2]),
  sprintf('<line x1="200" y1="250" x2="%.1f" y2="%.1f" stroke="#1D9E75" stroke-width="1.9" marker-end="url(#arrow-teal)"/>', short[, 1], short[, 2]))

# ---- bottom: fitted against true, one panel per chain
lim <- max(abs(c(truth, unlist(est)))) * 1.04
scatter <- function(x0, y0, w, h, e, sl) {
  half <- min(w, h) / 2 - 6
  px <- function(v) x0 + w / 2 + v / lim * half; py <- function(v) y0 + h / 2 - v / lim * half
  c(sprintf('<rect x="%g" y="%g" width="%g" height="%g" rx="4" fill="none" stroke="#c3c2b7" stroke-width="0.75"/>', x0, y0, w, h),
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#c3c2b7" stroke-width="0.75"/>', px(-lim), py(0), px(lim), py(0)),
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#c3c2b7" stroke-width="0.75"/>', px(0), py(-lim), px(0), py(lim)),
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888780" stroke-width="1" stroke-dasharray="3 3"/>', px(-lim), py(-lim), px(lim), py(lim)),
    sprintf('<circle cx="%.1f" cy="%.1f" r="1.9" fill="#1D9E75" fill-opacity="0.55"/>', px(truth), py(e)),
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#D85A30" stroke-width="2.5" stroke-linecap="round"/>', px(-lim), py(-lim * sl), px(lim), py(lim * sl)),
    sprintf('<text class="ts" x="%g" y="%g">↑ fitted</text>', x0 + 8, y0 + 18),
    sprintf('<text class="ts" x="%g" y="%g">true →</text>', x0 + w - 52, y0 + h - 8))
}
svg <- c('<?xml version="1.0" encoding="UTF-8"?>',
  '<svg xmlns="http://www.w3.org/2000/svg" width="680" height="640" viewBox="0 0 680 640" role="img" font-family="Helvetica, Arial, sans-serif">',
  '<title>Chains that agree with each other and miss the sheet</title>',
  '<desc>Top left: species arrows lying in a sheet spanned by two blends, as the model is meant to find. Top right: the same species with no sheet found; each arrow is fitted on its own and comes out much shorter than its true length, drawn in grey behind it. Bottom: from one fit, the fitted coefficients of the active embedding columns against their true values, for a chain where the blends carry the signal and a chain where the leftovers carry it. The second set lies along a much shallower line.</desc>',
  '<defs>',
  '<marker id="arrow-coral" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M2 1L8 5L2 9" fill="none" stroke="#D85A30" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/></marker>',
  '<marker id="arrow-teal" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M2 1L8 5L2 9" fill="none" stroke="#1D9E75" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/></marker>',
  '<style>.th{font-size:14px;font-weight:bold;fill:#1a1a19}.ts{font-size:12px;fill:#52514e}</style>',
  '</defs>',
  '<rect x="0" y="0" width="680" height="640" fill="#ffffff"/>',
  '<text class="th" x="40" y="40">The blends carry the signal</text>',
  '<text class="th" x="380" y="40">The leftovers carry the signal</text>',
  '<text class="ts" x="40" y="58">A sheet is found, and species share it</text>',
  '<text class="ts" x="380" y="58">No sheet: every species is fitted on its own</text>',
  '<g transform="translate(2.6,-2) scale(0.78)">',
  '<line x1="200" y1="250" x2="400" y2="250" stroke="#c3c2b7" stroke-width="1.2"/>',
  '<line x1="200" y1="250" x2="200" y2="90" stroke="#c3c2b7" stroke-width="1.2"/>',
  '<line x1="200" y1="250" x2="120" y2="330" stroke="#c3c2b7" stroke-width="1.2"/>',
  '<polygon points="352,196.4 256,332.4 48,303.6 144,167.6" fill="#7F77DD" fill-opacity="0.10" stroke="#7F77DD" stroke-opacity="0.5" stroke-width="0.9"/>',
  sprintf('<line x1="200" y1="250" x2="%.1f" y2="%.1f" stroke="#1D9E75" stroke-width="1.9" marker-end="url(#arrow-teal)"/>', tips[, 1], tips[, 2]),
  '<line x1="330.5" y1="216.8" x2="330.5" y2="186" stroke="#898781" stroke-width="1.2" stroke-dasharray="3 3"/>',
  '<line x1="200" y1="250" x2="297.5" y2="263.5" stroke="#D85A30" stroke-width="3.2" stroke-linecap="round" marker-end="url(#arrow-coral)"/>',
  '<line x1="200" y1="250" x2="245" y2="186.3" stroke="#D85A30" stroke-width="3.2" stroke-linecap="round" marker-end="url(#arrow-coral)"/>',
  '</g>',
  '<text class="th" x="242" y="212">Blend 1</text>', '<text class="th" x="170" y="132">Blend 2</text>',
  '<g transform="translate(342.6,-2) scale(0.78)">', leftover_panel, '</g>',
  '<text class="ts" x="540" y="120">True length, in grey</text>',
  '<text class="ts" x="402" y="270">The fitted arrows are far too short</text>',
  '<line x1="40" y1="282" x2="640" y2="282" stroke="#c3c2b7" stroke-width="0.75"/>',
  '<text class="th" x="40" y="306">One fit: a chain where the blends carry it</text>',
  '<text class="th" x="380" y="306">The other chain of the same fit</text>',
  scatter(40, 318, 260, 230, est[[1]], slope[1]), scatter(380, 318, 260, 230, est[[2]], slope[2]),
  sprintf('<text class="th" x="%g" y="%g">Slope %.2f</text>', 216, 350, slope[1]),
  sprintf('<text class="th" x="%g" y="%g">Slope %.2f</text>', 556, 350, slope[2]),
  '<text class="ts" x="40" y="576">Each dot is one species\' coefficient on one of the 10 embedding columns that carry signal.</text>',
  '<text class="ts" x="40" y="594">Where the leftovers carry the signal, many are shrunk to almost nothing: overall, a third of their true size.</text>',
  '<text class="ts" x="40" y="612">Chains that all end there agree with each other, so a check between chains passes.</text>',
  '<line x1="40" y1="628" x2="68" y2="628" stroke="#D85A30" stroke-width="2.5" stroke-linecap="round"/>',
  '<text class="ts" x="78" y="632">Line through the dots</text>',
  '<line x1="232" y1="628" x2="260" y2="628" stroke="#888780" stroke-width="1" stroke-dasharray="3 3"/>',
  '<text class="ts" x="270" y="632">Fitted equals true</text>',
  '<circle cx="412" cy="628" r="1.9" fill="#1D9E75" fill-opacity="0.55"/><circle cx="420" cy="626" r="1.9" fill="#1D9E75" fill-opacity="0.55"/><circle cx="428" cy="629" r="1.9" fill="#1D9E75" fill-opacity="0.55"/>',
  '<text class="ts" x="438" y="632">One coefficient</text>',
  '</svg>')
writeLines(svg, file.path(root, "figures/chains-agree-but-miss-the-sheet.svg"))
cat(sprintf("slopes: blends chain %.3f, leftovers chain %.3f; log sigma_emb %.2f and %.2f; %d coefficients\n",
            slope[1], slope[2], x$log_sigma_emb_chain_mean[1], x$log_sigma_emb_chain_mean[2], length(truth)))
