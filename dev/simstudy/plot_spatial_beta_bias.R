#!/usr/bin/env Rscript
# Recreate the static figure from the compact, saved analysis results.
script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")])
stopifnot(length(script) == 1L)
root <- dirname(normalizePath(script, mustWork = TRUE))
results <- file.path(root, "spatial-beta-binary-results")
main <- read.csv(file.path(results, "main-strata.csv"))
summary <- read.csv(file.path(results, "three-dataset-summary.csv"))
groups <- c("low", "medium", "high")
summary <- summary[match(groups, summary$stratum), ]
stopifnot(all(summary$datasets == 3L), nrow(main) == 9L)
png(file.path(root, "spatial-beta-probability-bias.png"),
    width = 1800, height = 1150, res = 180)
par(mar = c(7.3, 5.1, 4.7, 1.5), family = "sans", las = 1,
    col.axis = "#334155", col.lab = "#0f172a")
plot(NA, xlim = c(.55, 3.45), ylim = c(-6.2, 6.2), xaxt = "n",
     xlab = "", ylab = "Average error (percentage points)", bty = "n")
rect(.55, -5, 3.45, 5, col = "#eef6f1", border = NA)
abline(h = c(-5, 5), col = "#718a77", lty = 2)
abline(h = 0, col = "#94a3b8")
axis(1, at = 1:3, labels = c("Low\ntrue probability <20%",
     "Medium\n20% to 80%", "High\n>80%"), tick = FALSE, padj = .6)
for (g in seq_along(groups)) {
  values <- main$mean[main$stratum == groups[g]] * 100
  points(g + c(-.16, -.09, -.02), values, pch = 1, cex = 1.35,
         col = "#607d8b", lwd = 1.6)
  center <- summary$mean_bias[g] * 100
  se <- summary$between_dataset_se[g] * 100
  arrows(g + .13, center - se, g + .13, center + se,
         angle = 90, code = 3, length = .055, col = "#164e63", lwd = 2)
  points(g + .13, center, pch = 19, cex = 1.35, col = "#164e63")
  text(g + .13, center + if (center > 0) .9 else -.9,
       labels = sprintf("%+.2f", center), col = "#164e63", font = 2)
}
title("Spatial probability errors meet the provisional five-point target",
      sub = "", adj = 0, cex.main = 1.12, col.main = "#0f172a", line = 2.5)
mtext("Three independent datasets, one per spatial range; common species with abundant observations",
      side = 3, line = 1.1, adj = 0, cex = .85, col = "#475569")
legend("topleft", c("Individual dataset", "Mean across datasets, +/- one standard error"),
       pch = c(1, 19), col = c("#607d8b", "#164e63"), bty = "n", cex = .8)
mtext("Shading: provisional +/-5-point target. Positive = overestimate; negative = underestimate.",
      side = 1, line = 4.2, adj = 0, cex = .85, col = "#475569")
mtext("Bars reflect variation across only three datasets, not prediction intervals. Rare species were not tested.",
      side = 1, line = 5.5, adj = 0, cex = .8, col = "#475569")
invisible(dev.off())
