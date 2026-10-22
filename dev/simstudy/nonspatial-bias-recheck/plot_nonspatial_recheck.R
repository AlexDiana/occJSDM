#!/usr/bin/env Rscript
# Full render: Rscript plot_nonspatial_recheck.R STUDY_ROOT
# Complete collection figures only: append --collection-only.
args <- commandArgs(trailingOnly = TRUE)
collection_only <- "--collection-only" %in% args
unknown_flags <- setdiff(args[startsWith(args, "--")], "--collection-only")
if (length(unknown_flags)) stop("Unknown option: ", paste(unknown_flags, collapse = ", "))
paths <- args[!startsWith(args, "--")]
root <- normalizePath(if (length(paths)) paths[1] else ".")
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(root, "main-results/summary")
s <- read.csv(file.path(out, "summary.csv"))

check_complete <- function(data, expected_rows, label) {
  if (nrow(data) != expected_rows || any(is.na(data$datasets)) ||
      any(data$datasets != 10L) || any(!is.finite(data$between_dataset_se))) {
    stop(label, " needs all expected arms, each with 10 completed datasets. ",
         "Use --collection-only while PCR comparisons are still running.")
  }
}
with_error_units <- function(data) {
  data$bias_pp <- 100 * data$bias
  data$half_width <- 100 * qt(.975, data$datasets - 1) * data$between_dataset_se
  data
}
base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0), axis.text.y = element_text(size = 10))
model_names <- c(occupancy = "Classical occupancy", traits_isolated = "Two-stage eDNA")
bands <- c(all = "All probabilities", low = "True probability below 20%",
           medium = "True probability 20-80%", high = "True probability above 80%")

a <- s[s$scenario %in% names(model_names) & s$prior == "var2" &
         s$metric == "occupancy" & s$group %in% names(bands), ]
b <- read.csv(file.path(out, "known-occupancy-summary.csv"))
b <- b[b$scenario %in% names(model_names) & b$metric == "occupancy" &
         b$group %in% names(bands), ]
check_complete(a, 8L, "Imperfect-detection occupancy figure")
check_complete(b, 8L, "Exact-state occupancy figure")
a$observations <- "Imperfect detection"
b$observations <- "Exact occupancy states"
cols <- c("scenario", "group", "datasets", "bias", "between_dataset_se", "observations")
x <- rbind(a[cols], b[cols])
x$model <- unname(model_names[x$scenario])
x$band <- factor(bands[x$group], levels = rev(bands))
x <- with_error_units(x)
p <- ggplot(x, aes(bias_pp, band, colour = observations)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = bias_pp - half_width, xmax = bias_pp + half_width),
                orientation = "y", width = .12, position = position_dodge(width = .4)) +
  geom_point(size = 2.6, position = position_dodge(width = .4)) +
  facet_wrap(~model, ncol = 2) +
  scale_colour_manual(values = c("Exact occupancy states" = "#3b7b73", "Imperfect detection" = "#c2612e")) +
  labs(title = "Small overall occupancy errors hide errors at low and high probabilities",
       subtitle = "100 sites, 10 species, two field samples per site; current priors; no spatial effects",
       x = "Average error in probability (percentage points)\nNegative: too low                 0: correct on average                 Positive: too high",
       y = NULL, colour = NULL,
       caption = "The x-axis shows error, not the probability itself. Each point averages errors across 10 simulated datasets.\nExample: true probability 70%, estimated 73% = +3 percentage points; estimated 68% = -2 points.\nBars: 95% t intervals for average error across 10 simulated datasets. These are not model credible intervals.\nThe exact-state control removes detection uncertainty; each site still supplies only one occupancy state per species.") +
  base_theme
ggsave(file.path(out, "occupancy-bias.png"), p, width = 11, height = 6.5, dpi = 170, bg = "white")

c <- read.csv(file.path(out, "collection-contrast-summary.csv"))
c <- c[c$scenario %in% names(model_names) & c$metric == "theta_change_mean_to_plus1SD" &
         c$slope_sign %in% c("negative", "zero", "positive") &
         c$prior %in% c("var2", "var05", "var01"), ]
check_complete(c, 18L, "Collection-prior figure")
c$model <- unname(model_names[c$scenario])
c$sign <- factor(c$slope_sign, levels = c("negative", "zero", "positive"),
                 labels = c("Negative true effects", "Zero true effects", "Positive true effects"))
c$prior_label <- factor(c$prior, levels = c("var2", "var05", "var01"),
                        labels = c("Variance 2 (current)", "Variance 0.5", "Variance 0.1"))
c <- with_error_units(c)
p <- ggplot(c, aes(bias_pp, sign, colour = prior_label)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = bias_pp - half_width, xmax = bias_pp + half_width),
                orientation = "y", width = .1, position = position_dodge(width = .5)) +
  geom_point(size = 2.6, position = position_dodge(width = .5)) +
  facet_wrap(~model, ncol = 2) +
  scale_colour_manual(values = c("#336c91", "#ad8130", "#a44556")) +
  labs(title = "Tighter collection priors weaken real collection effects",
       subtitle = "Effect on collection probability when the covariate increases from its mean by one standard deviation",
       x = "Estimated change minus true change (percentage points)", y = NULL, colour = NULL,
       caption = "The x-axis shows error in the collection effect, averaged across 10 datasets, not the collection probability.\nExample: a true increase of 10 percentage points estimated as 13 gives an error of +3 points.\nZero means the change is correct on average; positive means the estimated change is higher, negative means lower.\nBars: 95% t intervals for average error across 10 paired datasets.\nFor positive effects, a negative error means an effect that is too small; for negative effects, a positive error means too small a decrease.") +
  base_theme
ggsave(file.path(out, "collection-prior-bias.png"), p, width = 11, height = 6.1, dpi = 170, bg = "white")

if (collection_only) {
  cat("Saved the two complete collection figures; PCR comparisons were not rendered.\n")
  quit(status = 0)
}

# The historical comparison includes K30 as an information diagnostic. It
# must never be presented as an operational recommendation under a cap of 6.
q_scenarios <- c("qnear_K3", "qnear_K6", "qnear_K30", "qfar_K3", "qfar_K6", "qfar_K30")
q_labels <- c("Low contamination: 3 per primer", "Low contamination: 6 per primer",
              "Low contamination: 30 (diagnostic)", "High contamination: 3 per primer",
              "High contamination: 6 per primer", "High contamination: 30 (diagnostic)")
q <- s[s$scenario %in% q_scenarios & s$prior == "q20" &
         s$metric %in% c("q", "q_nominal") & s$group == "all", ]
check_complete(q, 12L, "False-positive truth figure (K3, K6, K30)")
q$comparison <- ifelse(q$metric == "q", "Correct: positive read result", "Old: contamination event")
q$scenario_label <- factor(q$scenario, levels = rev(q_scenarios), labels = rev(q_labels))
q <- with_error_units(q)
p <- ggplot(q, aes(bias_pp, scenario_label, colour = comparison)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = bias_pp - half_width, xmax = bias_pp + half_width),
                orientation = "y", width = .1, position = position_dodge(width = .4)) +
  geom_point(size = 2.7, position = position_dodge(width = .4)) +
  scale_colour_manual(values = c("Correct: positive read result" = "#3b7b73", "Old: contamination event" = "#a44556")) +
  labs(title = "Compare false-positive estimates with the probability of positive reads",
       subtitle = "PCR counts are per primer; two primers. Operational maximum: 6 per primer; K30 is diagnostic only.",
       x = "Average false-positive estimate minus comparison truth (percentage points)",
       y = NULL, colour = NULL,
       caption = "The x-axis shows the estimated probability minus the comparison value, averaged across 10 datasets.\nFor the correct (green) comparison: true 20%, estimated 23% = +3 percentage points; estimated 18% = -2 points.\nFor green points, zero means correct on average; positive means too high and negative means too low.\nBars: 95% t intervals for average error across 10 datasets, using current priors.\nSome contamination events generate zero reads. The fitted binary model estimates the chance of a positive read result.") +
  base_theme
ggsave(file.path(out, "false-positive-truth.png"), p, width = 12, height = 7.6, dpi = 170, bg = "white")

# Operational comparison: each panel shows both feasible PCR allocations
# and all three prior choices. All truth probabilities are post-threshold.
practical_scenarios <- c("qnear_K3", "qnear_K6", "qfar_K3", "qfar_K6")
practical <- s[s$scenario %in% practical_scenarios & s$prior %in% c("q20", "q9", "p32") &
                ((s$metric %in% c("p", "q") & s$group == "all") |
                 (s$metric == "occupancy" & s$group %in% c("low", "high"))), ]
check_complete(practical, 48L, "Operational PCR figure (K3 versus K6)")
practical$contamination <- ifelse(startsWith(practical$scenario, "qnear"), "Low contamination", "High contamination")
practical$contamination <- factor(practical$contamination, levels = c("Low contamination", "High contamination"))
practical$replicates <- factor(sub(".*_K", "", practical$scenario), levels = c("6", "3"),
                               labels = c("6 per primer", "3 per primer"))
metric_labels <- c(q = "False positives (q)", p = "True detections (p)",
                   occupancy_low = "Occupancy: true probability below 20%",
                   occupancy_high = "Occupancy: true probability above 80%")
metric_key <- ifelse(practical$metric == "occupancy", paste0("occupancy_", practical$group), practical$metric)
practical$quantity <- factor(metric_labels[metric_key], levels = metric_labels)
practical$prior_label <- factor(practical$prior, levels = c("q20", "q9", "p32"),
                                labels = c("Default priors: p Beta(5,1); q Beta(1,20)",
                                           "Change q only (false positives): p Beta(5,1); q Beta(1,9)",
                                           "Change p only (true detections): p Beta(3,2); q Beta(1,20)"))
practical <- with_error_units(practical)
p <- ggplot(practical, aes(bias_pp, replicates, colour = prior_label, shape = prior_label)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = bias_pp - half_width, xmax = bias_pp + half_width),
                orientation = "y", width = .14, position = position_dodge(width = .6)) +
  geom_point(size = 2.5, position = position_dodge(width = .6)) +
  facet_wrap(vars(quantity, contamination), ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c("#336c91", "#3b7b73", "#a44556")) +
  scale_shape_manual(values = c(16, 15, 17)) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(title = "Three or six PCR replicates per primer: detection and occupancy errors",
       subtitle = "Two primers; no spatial effects. Six PCR replicates per primer is the operational maximum.",
       x = "Average error in probability (percentage points)\nNegative: too low                 0: correct on average                 Positive: too high",
       y = NULL, colour = NULL, shape = NULL,
       caption = "The x-axis shows error, not the probability itself. Each point averages errors across 10 simulated datasets.\nExample: true detection probability 70%, estimated 73% = +3 percentage points; estimated 68% = -2 points.\nAn x-value of 3 means three percentage points too high, not a detection probability of 3%.\nBars: 95% t intervals for average error across 10 datasets. Each panel uses an axis range suited to its errors.\np: positive PCR result when species DNA is present in the field sample; q: positive result when it is absent.\nAll other priors are unchanged. Beta(a,b) describes a prior distribution, not a fixed detection rate.\nLow/high occupancy panels group by the true site probability, not by species prevalence. K30 is excluded from this operational comparison.") +
  base_theme + theme(panel.spacing = grid::unit(1.2, "lines"), strip.text = element_text(size = 10),
                     legend.text = element_text(size = 10), plot.caption = element_text(size = 10))
ggsave(file.path(out, "operational-pcr-bias.png"), p, width = 12, height = 13, dpi = 170, bg = "white")
cat("Saved four figures, including the K3/K6 operational comparison and the separately labeled K30 diagnostic.\n")
