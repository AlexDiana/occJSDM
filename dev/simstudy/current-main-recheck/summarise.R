#!/usr/bin/env Rscript
# Compare the already selected fits. Selection never uses these errors.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
repo <- normalizePath(arg("repo", "."))
study <- normalizePath(arg("study"))
out <- file.path(repo, "dev/simstudy/current-main-recheck/results")
dir.create(out, showWarnings = FALSE)
selection <- readRDS(file.path(study, "selection.rds"))
verification <- read.csv(file.path(study, "verification-selected.csv"))
stopifnot(nrow(selection) == 110L, !anyDuplicated(selection$key),
          nrow(verification) == 110L, setequal(selection$key, verification$key),
          all(verification$max_probability_difference < 1e-12),
          all(verification$groups_checked == 4L))

normalise_groups <- function(r, job) {
  g <- r$groups
  if (job$family == "jsdm") {
    g$metric <- ifelse(g$scope == "original100", "occupancy_original_sites", "occupancy_all_sites")
    g$n_elements <- g$cells
  } else {
    # The baseline contains exactly the original 100 sites.
    if (job$arm == "baseline" && !"occupancy_original_sites" %in% g$metric) {
      extra <- g[g$metric == "occupancy" & g$group %in% c("all", "low", "medium", "high"), ]
      stopifnot(nrow(extra) == 4L)
      extra$metric <- "occupancy_original_sites"
      g <- rbind(g, extra)
    }
  }
  g$group[g$group == "middle"] <- "medium"
  keep <- c("metric", "group", "n_elements", "truth", "estimate", "bias", "mae", "rmse", "rhat", "ess_mean", "mcse", "chain_gap")
  stopifnot(all(keep %in% names(g)), !anyDuplicated(g[c("metric", "group")]))
  g[, keep]
}

rows <- list(); provenance <- list()
for (i in seq_len(nrow(selection))) {
  s <- selection[i, ]
  paths <- c(historical = s$old_result_file, current = s$result_file,
             historical_initial = s$old_initial_file,
             current_initial = file.path(study, "initial", paste0(s$key, "-result.rds")))
  v <- verification[match(s$key, verification$key), ]
  stopifnot(identical(unname(tools::md5sum(s$result_file)), v$result_md5))
  for (version in names(paths)) {
    r <- readRDS(paths[[version]])
    g <- normalise_groups(r, s)
    rows[[paste(s$key, version)]] <- cbind(s[c("key", "family", "arm", "scenario", "replicate")],
                                         version = version, g, row.names = NULL)
    provenance[[paste(s$key, version)]] <- data.frame(key = s$key, version = version,
      file = paths[[version]], md5 = unname(tools::md5sum(paths[[version]])),
      fit_seconds = as.numeric(difftime(r$finished, r$started, units = "secs")),
      warnings = length(r$warnings))
  }
}
scores <- do.call(rbind, rows); rownames(scores) <- NULL
provenance <- do.call(rbind, provenance); rownames(provenance) <- NULL
primary <- subset(scores, metric == "occupancy_original_sites")
stopifnot(nrow(primary) == 110L * 4L * 4L)

paired <- function(old, current, label) {
  keys <- c("key", "family", "arm", "scenario", "replicate", "metric", "group")
  p <- merge(subset(scores, version == old), subset(scores, version == current),
             by = keys, suffixes = c("_historical", "_current"), sort = FALSE)
  stopifnot(all(p$n_elements_historical == p$n_elements_current),
            max(abs(p$truth_historical - p$truth_current)) < 1e-12)
  for (m in c("bias", "mae", "rmse")) p[[paste0("delta_", m)]] <- p[[paste0(m, "_current")]] - p[[paste0(m, "_historical")]]
  p$comparison <- label
  p
}
pairs <- rbind(paired("historical", "current", "selected"),
               paired("historical_initial", "current_initial", "initial"))
stopifnot(sum(pairs$metric == "occupancy_original_sites") == 110L * 4L * 2L)

mean_ci <- function(x) {
  stopifnot(length(x) == 10L, all(is.finite(x)))
  half <- qt(.975, df = 9L) * sd(x) / sqrt(10)
  c(mean = mean(x), lower = mean(x) - half, upper = mean(x) + half)
}
# Rare-species groups occur in only one community; preserve their raw rows but
# do not treat their repeated fits as ten independent rare-species communities.
aggregate_pairs <- subset(pairs, !group %in% c("rare_below_20pct", "common_20pct_or_more") & metric != "q_nominal")
summary_rows <- lapply(split(aggregate_pairs, interaction(aggregate_pairs$comparison, aggregate_pairs$family, aggregate_pairs$scenario,
                                               aggregate_pairs$arm, aggregate_pairs$metric, aggregate_pairs$group, drop = TRUE)), function(p) {
  stopifnot(nrow(p) == 10L, setequal(p$replicate, 1:10))
  ans <- p[1, c("comparison", "family", "scenario", "arm", "metric", "group")]
  ans$communities <- 10L
  for (m in c("bias", "mae", "rmse")) {
    label <- if (m == "rmse") "mean_community_rmse" else m
    for (version in c("historical", "current")) ans[[paste0(version, "_", label)]] <- mean(p[[paste0(m, "_", version)]])
    ci <- mean_ci(p[[paste0("delta_", m)]])
    for (n in names(ci)) ans[[paste0("delta_", label, "_", n)]] <- unname(ci[n])
  }
  ans
})
summary <- do.call(rbind, summary_rows); rownames(summary) <- NULL
write.csv(scores, file.path(out, "community-scores.csv"), row.names = FALSE)
write.csv(pairs, file.path(out, "paired-community-changes.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "paired-summary.csv"), row.names = FALSE)

# Preserve within-version design comparisons as well as between-version changes.
# Positive reduction means the alternative design has lower absolute error.
design_pairs <- list()
for (version in c("historical", "current")) for (family in c("jsdm", "design")) {
  comparisons <- if (family == "jsdm") {
    list(c("sites100", "sites300"), c("sites100", "sites1000"), c("sites300", "sites1000"))
  } else {
    list(c("baseline", "field4"), c("baseline", "sites300"), c("baseline", "knownU"), c("field4", "sites300"))
  }
  for (comparison in comparisons) {
    reference <- primary[primary$version == version & primary$family == family & primary$arm == comparison[1], ]
    alternative <- primary[primary$version == version & primary$family == family & primary$arm == comparison[2], ]
    d <- merge(reference, alternative, by = c("family", "scenario", "replicate", "group"),
               suffixes = c("_reference", "_alternative"))
    stopifnot(max(abs(d$truth_reference - d$truth_alternative)) < 1e-12,
              all(d$n_elements_reference == d$n_elements_alternative))
    design_pairs[[paste(version, family, paste(comparison, collapse = ":"))]] <- data.frame(
      d[c("family", "scenario", "replicate", "group")], version = version,
      reference = comparison[1], alternative = comparison[2],
      reduction_mae = d$mae_reference - d$mae_alternative)
  }
}
design_pairs <- do.call(rbind, design_pairs); rownames(design_pairs) <- NULL
design_summary <- do.call(rbind, lapply(split(design_pairs, interaction(design_pairs$version,
  design_pairs$family, design_pairs$scenario, design_pairs$reference, design_pairs$alternative,
  design_pairs$group, drop = TRUE)), function(d) {
    stopifnot(nrow(d) == 10L, setequal(d$replicate, 1:10))
    ci <- mean_ci(d$reduction_mae)
    data.frame(d[1, c("version", "family", "scenario", "reference", "alternative", "group")],
               communities = 10L, improved = sum(d$reduction_mae > 0),
               reduction_mae = ci["mean"], lower = ci["lower"], upper = ci["upper"])
  }))
write.csv(design_pairs, file.path(out, "design-paired-community-changes.csv"), row.names = FALSE)
write.csv(design_summary, file.path(out, "design-paired-summary.csv"), row.names = FALSE)
write.csv(provenance, file.path(out, "result-provenance.csv"), row.names = FALSE)
write.csv(selection, file.path(out, "selected-manifest.csv"), row.names = FALSE)
file.copy(file.path(study, c("verification-selected.csv", "long-selection.csv")), out, overwrite = TRUE)
settings <- readRDS(file.path(study, "initial/settings.rds"))
postprocessing <- file.path(repo, "dev/simstudy/current-main-recheck",
                           c("select.R", "verify.R", "summarise.R", "diagnostics.R", "chain-sensitivity.R", "protocol.md"))
hashes <- c(settings$source_hashes, tools::md5sum(postprocessing))
write.csv(data.frame(path = names(hashes), md5 = unname(hashes)),
          file.path(out, "source-hashes.csv"), row.names = FALSE)
writeLines(trimws(c(paste("Production revision:", settings$source_revision), capture.output(settings$session)), which = "right"),
           file.path(out, "environment.txt"))

scenario_labels <- c(perfect = "Binary JSDM: perfect detection",
                     qnear_K6 = "Two-stage: low contamination",
                     qfar_K6 = "Two-stage: high contamination")
arm_labels <- c(sites100 = "100 sites", sites300 = "300 sites", sites1000 = "1,000 sites",
                baseline = "Baseline", field4 = "4 field samples", knownU = "Known site factors")
display <- function(x) {
  probability_columns <- intersect(c("bias", "mae", "rmse", "historical_mae", "current_mae",
                                    "delta_mae_mean", "delta_mae_lower", "delta_mae_upper"), names(x))
  x[probability_columns] <- lapply(x[probability_columns], function(y) y * 100)
  x$Scenario <- factor(scenario_labels[x$scenario], levels = unname(scenario_labels))
  x$Arm <- factor(arm_labels[x$arm], levels = unname(arm_labels[c("baseline", "sites100", "field4", "sites300", "sites1000", "knownU")]))
  x
}
g <- display(subset(primary, version %in% c("historical", "current") & group == "all"))
plot_rows <- do.call(rbind, lapply(split(g, interaction(g$family, g$scenario, g$arm, g$version, drop = TRUE)), function(d) {
  z <- d[1, c("Scenario", "Arm", "version")]; ci <- mean_ci(d$mae)
  z$mae <- ci["mean"]; z$lower <- ci["lower"]; z$upper <- ci["upper"]; z
}))
suppressPackageStartupMessages(library(ggplot2))
base_theme <- theme_bw(base_size = 11) + theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1))
p <- ggplot(plot_rows, aes(Arm, mae, colour = version, group = version)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(.35), width = .15) +
  geom_point(position = position_dodge(.35), size = 2.5) +
  facet_wrap(~ Scenario, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(historical = "#666666", current = "#0072B2"),
                      breaks = c("historical", "current"), labels = c("Archived code", "Current main")) +
  labs(x = NULL, y = "Mean absolute error (percentage points)", colour = NULL,
       title = "How much did fitted-site probability recovery change?",
       subtitle = "Same original 100 sites and 10 communities in each arm",
       caption = "Bars: 95% t intervals across community errors. These are not posterior coverage estimates.") + base_theme
ggsave(file.path(out, "probability-error.png"), p, width = 12, height = 4.8, dpi = 180)

d <- display(subset(summary, comparison == "selected" & metric == "occupancy_original_sites" & group == "all"))
d <- d[order(d$Scenario, d$Arm), ]
p <- ggplot(d, aes(Arm, delta_mae_mean)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = delta_mae_lower, ymax = delta_mae_upper), width = .15, colour = "#0072B2") +
  geom_point(size = 2.5, colour = "#0072B2") + facet_wrap(~ Scenario, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Change in absolute error (percentage points)",
       title = "Paired change on identical communities",
       subtitle = "Below zero: lower error with current code; above zero: higher error",
       caption = "Bars: 95% t intervals for the mean paired change across 10 communities.") + base_theme
ggsave(file.path(out, "paired-change.png"), p, width = 12, height = 4.8, dpi = 180)

b <- display(subset(primary, version %in% c("historical", "current") & group %in% c("low", "high")))
b$Band <- factor(b$group, levels = c("low", "high"), labels = c("True probability below 0.2", "True probability above 0.8"))
p <- ggplot(b, aes(Arm, bias, colour = version, group = version)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  stat_summary(fun = mean, geom = "point", position = position_dodge(.35), size = 2.5) +
  facet_grid(Band ~ Scenario, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c(historical = "#666666", current = "#0072B2"),
                      breaks = c("historical", "current"), labels = c("Archived code", "Current main")) +
  labs(x = NULL, y = "Mean signed error (percentage points)", colour = NULL,
       title = "Do extreme probabilities still move toward the middle?",
       caption = "Positive error in the low band and negative error in the high band indicate inward bias.") + base_theme
ggsave(file.path(out, "probability-band-bias.png"), p, width = 12, height = 6, dpi = 180)

# A generated table for the report, never a manually maintained copy of results.
tab <- d[c("Scenario", "Arm", "historical_mae", "current_mae", "delta_mae_mean", "delta_mae_lower", "delta_mae_upper")]
names(tab) <- c("Scenario", "Design", "Archived MAE", "Current MAE", "Paired change", "95% lower", "95% upper")
writeLines(c("All errors and changes are in percentage points.", "",
  as.character(knitr::kable(tab, digits = 3, row.names = FALSE, format = "pipe"))), file.path(out, "comparison-table.md"))
cat("Exported", nrow(scores), "community scores and", nrow(summary), "paired summaries.\n")
print(tab, row.names = FALSE)
