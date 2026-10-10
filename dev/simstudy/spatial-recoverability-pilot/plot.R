# Copies compact results into the study folder and draws the recoverability map.
# run-r -f plot.R --args --study=DIR
args <- commandArgs(trailingOnly = TRUE)
study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
out <- "dev/simstudy/spatial-recoverability-pilot/results"; dir.create(out, showWarnings = FALSE)
for (f in c("species.csv", "aggregate.csv", "gate.csv", "map.json", "factors.csv", "conditionals.csv", "geometry.csv", "screen.csv", "jobs.csv"))
  file.copy(file.path(study, "scores", f), file.path(out, f), overwrite = TRUE)
file.copy(file.path(study, "input-manifest.csv"), file.path(out, "input-manifest.csv"), overwrite = TRUE)
file.copy(file.path(study, "timing.json"), file.path(out, "timing.json"), overwrite = TRUE)
cfg <- readRDS(file.path(study, "config.rds"))
writeLines(jsonlite::toJSON(cfg[c("created", "communities", "thresholds", "budget", "r_version", "r_home", "packages", "sources", "git")],
                            auto_unbox = TRUE, pretty = TRUE), file.path(out, "config.json"))

agg <- read.csv(file.path(study, "scores", "aggregate.csv"), stringsAsFactors = FALSE)
gate <- read.csv(file.path(study, "scores", "gate.csv"), stringsAsFactors = FALSE)
cols <- c(pass = "#8fd3c7", fail = "#f4a582", unresolved = "#cccccc")
png(file.path(out, "recoverability-map.png"), width = 2200, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 6, 3.5, 1))
for (tr in c("t1", "t2")) {
  g <- gate[gate$truth == tr, ]; a <- agg[agg$truth == tr, ]
  rows <- if (tr == "t1") sort(unique(g$range)) else c("short", "medium", "long")
  key <- if (tr == "t1") g$range else g$group
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.5, length(rows) + 0.5), xaxt = "n", yaxt = "n", xlab = "training sites (random uniform)",
       ylab = if (tr == "t1") "shared-field range" else "dispersal class",
       main = if (tr == "t1") "Truth 1: common gradients (joint oracle)" else "Truth 2: trait-linked dispersal fields (per-species oracle)")
  axis(1, 1:3, c(100, 300, 900)); axis(2, seq_along(rows), if (tr == "t1") sprintf("l = %.2f", rows) else rows, las = 1)
  for (i in seq_len(nrow(g))) {
    x <- match(g$n[i], c(100, 300, 900)); y <- match(key[i], rows)
    rect(x - 0.47, y - 0.47, x + 0.47, y + 0.47, col = cols[g$decision[i]], border = "white")
    text(x, y + 0.3, toupper(g$decision[i]), font = 2, cex = 0.8)
    if (g$decision[i] != "unresolved") {
      text(x, y + 0.1, sprintf("pattern r = %.2f (%d/10 >= 0.4)", g$pattern_median[i], g$pattern_count[i]), cex = 0.68)
      text(x, y - 0.07, sprintf("slope %.2f", g$slope_median[i]), cex = 0.68)
      text(x, y - 0.22, sprintf("headroom %.1f pp [%.1f, %.1f]", g$headroom_mean[i], g$headroom_low[i], g$headroom_high[i]), cex = 0.68)
    }
    sel <- a$arrangement == "clustered" & a$n == g$n[i] & (if (tr == "t1") a$range %in% g$range[i] else a$group == g$group[i])
    text(x, y - 0.38, sprintf("clustered: r = %.2f, slope %.2f", median(a$contrast_correlation[sel]), median(a$field_slope[sel])), cex = 0.62, col = "grey25")
  }
}
dev.off()
cat("Wrote", file.path(out, "recoverability-map.png"), "\n")
