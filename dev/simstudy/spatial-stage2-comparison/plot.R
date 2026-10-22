# Figures and compact results for stage 2. Run through run-r:
#   run-r -f plot.R --args --study=DIR
# Copies the scores, timing records and a runtime record into results/, and draws headroom-share.png and pairs.png.
D2 <- "dev/simstudy/spatial-stage2-comparison"
arm_labels <- c("occJSDM/env" = "env", "occJSDM/gp" = "gp", "Hmsc/spatial2" = "Hmsc 2", "Hmsc/spatial6" = "Hmsc 6", "gllvm/spatial2" = "gllvm 2")
cell_labels <- function(cell) {
  truth <- if (startsWith(cell, "t1")) "Truth 1, range 0.08" else "Truth 2"
  sprintf("%s, n = %s", truth, sub(".*-n", "", cell))
}

plot_headroom_share <- function(agg, cfg, file) {
  keys <- paste(cfg$arms$package, cfg$arms$arm, sep = "/"); clip <- 1.5
  png(file, width = 1500, height = 1150, res = 150)
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.6, 2.6, 0.8), oma = c(0, 0, 2.4, 0)); on.exit({ par(op); dev.off() })
  for (cell in cfg$cells) {
    d <- agg[agg$cell == cell, ]
    plot(NA, xlim = c(0.5, length(keys) + 0.5), ylim = c(-clip, clip), xaxt = "n", xlab = "", ylab = "Headroom share (oracle = 1, env arm = 0)", main = cell_labels(cell))
    axis(1, at = seq_along(keys), labels = arm_labels[keys], cex.axis = 0.9)
    abline(h = 0, col = "grey40", lty = 2); abline(h = 1, col = "grey40", lty = 2)
    for (i in seq_along(keys)) {
      r <- d[d$key == keys[i], ]; y <- r$headroom_share; off <- (seq_along(y) - (length(y) + 1) / 2) * 0.05
      keep <- is.finite(y); beyond <- keep & abs(y) > clip; yy <- pmax(pmin(y, clip), -clip)
      pass <- r$status == "pass"
      points(i + off[keep & !beyond], yy[keep & !beyond], pch = ifelse(pass[keep & !beyond], 16, 1), col = "#2c5d8f", cex = 1.1)
      if (any(beyond)) points(i + off[beyond], yy[beyond], pch = ifelse(y[beyond] > 0, 2, 6), col = "#b24a2a", cex = 1.1)
      if (any(keep)) segments(i - 0.3, median(y[keep]), i + 0.3, median(y[keep]), lwd = 3, col = "black")
      if (any(!keep)) mtext(sprintf("%d NA", sum(!keep)), side = 1, at = i, line = 2.2, cex = 0.7)
    }
  }
  mtext("Held-out headroom share by arm: one point per community, bar is the median", outer = TRUE, cex = 1.05, font = 2)
  mtext("Filled = status pass, open = any other status; triangles are beyond the plotted range", outer = TRUE, side = 1, line = -0.9, cex = 0.7, adj = 0.5)
  invisible(file)
}

plot_pairs <- function(pairs, cfg, file) {
  png(file, width = 1500, height = 1250, res = 150)
  op <- par(mar = c(4.2, 17.5, 6.6, 1)); on.exit({ par(op); dev.off() })
  pairs <- pairs[order(match(pairs$cell, cfg$cells), match(paste(pairs$a, pairs$b), vapply(cfg$comparisons, paste, "", collapse = " "))), ]
  n <- nrow(pairs); y <- rev(seq_len(n)); lim <- max(0.45, max(abs(c(pairs$lower, pairs$upper, pairs$estimate)), na.rm = TRUE) * 1.05)
  plot(NA, xlim = c(-lim, lim), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "Paired difference in headroom share (second arm minus first arm)", ylab = "")
  title(main = "Paired comparisons within each cell", line = 4.8)
  rect(-cfg$delta, 0.5, cfg$delta, n + 0.5, col = "#e8eef5", border = NA); abline(v = 0, col = "grey40")
  labels <- sprintf("%s: %s vs %s", pairs$cell, arm_labels[pairs$a], arm_labels[pairs$b]); axis(2, at = y, labels = labels, las = 1, cex.axis = 0.75)
  col <- ifelse(pairs$verdict == "inconclusive", "grey55", ifelse(pairs$verdict == "no clear difference", "#2c5d8f", "#b24a2a"))
  segments(pairs$lower, y, pairs$upper, y, col = col, lwd = 2); points(pairs$estimate, y, pch = 16, col = col, cex = 1.1)
  text(rep(lim, n), y, sprintf("n=%d", pairs$paired), adj = c(1, 0.5), cex = 0.65, col = "grey30")
  legend("topleft", inset = c(0, -0.115), xpd = NA, legend = c("clear (estimate beyond 0.10, interval excludes 0)", "no clear difference", "inconclusive (fewer than six pairs)"), col = c("#b24a2a", "#2c5d8f", "grey55"), pch = 16, lwd = 2, cex = 0.7, bg = "white")
  invisible(file)
}

copy_results <- function(study, results_dir) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  for (f in c("species.csv", "aggregate.csv", "pairs.csv", "jobs.csv", "reading.json")) file.copy(file.path(study, "scores", f), file.path(results_dir, f), overwrite = TRUE)
  for (f in c(list.files(study, pattern = "^timing-step[0-9]+\\.json$"), "microtiming.json")) if (file.exists(file.path(study, f))) file.copy(file.path(study, f), file.path(results_dir, f), overwrite = TRUE)
  cfg <- readRDS(file.path(study, "config.rds"))
  runtime <- list(r_version = cfg$r_version, packages = as.list(cfg$packages), frozen = cfg$frozen, frozen_at = format(cfg$frozen_at), frozen_from = cfg$frozen_from,
                  git_at_prepare = cfg$git, pilot_study = cfg$pilot_study, created = format(cfg$created),
                  provisioning = "Route B: Hmsc 3.3.7 and gllvm 2.0.15 installed as CRAN binaries into the study library under the pinned R 4.5.0, TMB 1.9.22 from the user library; the versions the first comparison used",
                  source_hashes = as.list(cfg$sources))
  writeLines(jsonlite::toJSON(runtime, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null"), file.path(results_dir, "runtime.json"))
  invisible(results_dir)
}

if (sys.nframe() == 0L && !interactive()) {
  source(file.path(D2, "config.R"))
  args <- commandArgs(trailingOnly = TRUE); study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
  results <- file.path(D2, "results"); cfg <- readRDS(file.path(study, "config.rds"))
  copy_results(study, results)
  agg <- read.csv(file.path(study, "scores", "aggregate.csv"), stringsAsFactors = FALSE); pairs <- read.csv(file.path(study, "scores", "pairs.csv"), stringsAsFactors = FALSE)
  plot_headroom_share(agg, cfg, file.path(results, "headroom-share.png")); plot_pairs(pairs, cfg, file.path(results, "pairs.png"))
  cat("Plots and results written to", results, "\n")
}
