#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) { hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }; substring(hit, nchar(name) + 4L) }
repo <- normalizePath(option("repo", ".")); study <- normalizePath(option("study"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
out <- file.path(repo, "dev/simstudy/spatial-design-sweep/results"); dir.create(out, showWarnings = FALSE)
summary <- file.path(study, "summary-final")
for (f in c("aggregate.csv", "reading.csv", "range-reading.csv", "paired.csv", "selected-fits.csv", "groups.csv",
            "species.csv", "field.csv", "range.csv", "amplitude.csv", "lattice.csv", "oracle.csv", "long-run-sensitivity.csv"))
  stopifnot(file.copy(file.path(summary, f), out, overwrite = TRUE))
# The committed copies of the two file-listing tables carry paths relative to the study directory;
# verify.R reads the study copies, which keep the absolute paths.
relative_to_study <- function(path) { prefix <- paste0(study, "/")
  stopifnot(startsWith(path, prefix)); substring(path, nchar(prefix) + 1L) }
for (spec in list(c("selected-fits.csv", "selected_file"), c("oracle.csv", "file"))) {
  d <- read.csv(file.path(out, spec[1])); d[[spec[2]]] <- relative_to_study(d[[spec[2]]])
  write.csv(d, file.path(out, spec[1]), row.names = FALSE)
}
stopifnot(file.copy(file.path(study, "design-statistics.csv"), out, overwrite = TRUE))
stopifnot(file.copy(file.path(study, "audit/audit.csv"), file.path(out, "audit.csv"), overwrite = TRUE))
order_arr <- c("spread", "pairs", "clustered", "grid")
agg <- read.csv(file.path(out, "aggregate.csv")) |> mutate(arrangement = factor(arrangement, order_arr),
  arm = factor(arm, c("oracle", "binary", "two_stage"), c("Oracle: true states, known parameters", "occJSDM: true states", "occJSDM: eDNA survey")),
  group = factor(group, c("prevalence_5pct", "prevalence_25pct", "prevalence_75pct")),
  reduction = 1 - centred_rmse / zero_field_rmse)
p <- ggplot(agg, aes(arrangement, reduction, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(position = position_dodge(.5), size = 2.5) +
  facet_wrap(~ group) + theme_bw(base_size = 12) +
  labs(x = "Site arrangement, 100 sites each", y = "Field RMSE reduction relative to a flat field",
       colour = NULL, title = "How much of the spatial field each design recovers",
       caption = "Mean of three communities per point. Zero means no better than assuming no field.") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "field-recovery.png"), p, width = 10, height = 6, dpi = 150)
lat <- read.csv(file.path(out, "lattice.csv")) |> filter(bin != "all") |>
  mutate(arrangement = factor(arrangement, order_arr), bin = factor(bin, c("up to 0.02", "0.02 to 0.05", "0.05 to 0.1", "above 0.1"))) |>
  group_by(arrangement, arm, spatial_term, bin) |> summarise(mae = mean(mae), .groups = "drop")
p <- ggplot(lat, aes(bin, 100 * mae, colour = spatial_term, group = spatial_term)) + geom_line() + geom_point() +
  facet_grid(arm ~ arrangement) + theme_bw(base_size = 11) +
  labs(x = "Distance from the nearest surveyed site", y = "Mean absolute prediction error (points)", colour = "Spatial term",
       title = "Prediction at unsurveyed locations by distance from the survey") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(out, "lattice-prediction.png"), p, width = 11, height = 6, dpi = 150)
rng <- read.csv(file.path(out, "range.csv")) |> mutate(arrangement = factor(arrangement, order_arr))
p <- ggplot(rng, aes(arrangement, mass_within_one_step, colour = arm)) + geom_hline(yintercept = .5, linetype = "dashed") +
  geom_point(position = position_dodge(.4), size = 2.5) + theme_bw(base_size = 12) + ylim(0, 1) +
  labs(x = NULL, y = "Posterior mass within one grid step of the true range", colour = "Arm", title = "Range recovery by arrangement")
ggsave(file.path(out, "range-recovery.png"), p, width = 8, height = 5, dpi = 150)
cat("Figures and compact results written to", out, "\n")
