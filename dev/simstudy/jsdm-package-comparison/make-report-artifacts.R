args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
out <- file.path(code, "pilot-results")
dir.create(out, showWarnings = FALSE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
files <- list.files(file.path(root, "checks"), "-summary[.]rds$", full.names = TRUE)
get_number <- function(s, name) if (is.null(s[[name]])) NA_real_ else as.numeric(s[[name]])
rows <- lapply(files, function(path) {
  label <- sub("-summary[.]rds$", "", basename(path))
  s <- readRDS(path)
  r <- readRDS(file.path(root, "fits", paste0(label, "-record.rds")))
  data.frame(fit = label, package = s$package, seconds = r$elapsed_seconds,
    process_note = if (startsWith(label, "sjSDM-attempt-2"))
      "Nonzero exit after numeric fit and success record saved; fresh-process checks passed" else "Completed",
    rhat = get_number(s, "max_rhat"), bulk_ess = get_number(s, "min_bulk_ess"),
    tail_ess = get_number(s, "min_tail_ess"), objective = if (s$package == "sjSDM")
      get_number(s, "checked_loglik") else get_number(s, "native_loglik"),
    maximum_gradient = get_number(s, "max_gradient"),
    maximum_coefficient = get_number(s, "max_abs_beta"),
    maximum_residual_sd = get_number(s, "max_residual_sd"),
    warnings = paste(r$warnings, collapse = "; "))
})
write.csv(bind_rows(rows), file.path(out, "all-fit-attempts.csv"), row.names = FALSE)
selected <- readRDS(file.path(root, "results/selection.rds"))$selected
for (package in c("occJSDM", "Hmsc")) {
  f <- file.path(root, "checks", paste0(selected[[package]], "-diagnostics.csv"))
  file.copy(f, file.path(out, paste0(package, "-diagnostics.csv")), overwrite = TRUE)
}
for (f in c("overall-errors.csv", "errors-by-band.csv", "errors-by-species.csv",
            "integration-checks.csv", "selected-fits.tsv", "predictions.rds")) {
  stopifnot(file.copy(file.path(root, "results", f), out, overwrite = TRUE))
}
file.copy(file.path(root, "truth/generating-parameters.csv"), out, overwrite = TRUE)
file.copy(file.path(root, "checks/diagnostic-environment.tsv"), out, overwrite = TRUE)
cells <- readRDS(file.path(out, "predictions.rds")) |>
  mutate(package = factor(package, levels = c("occJSDM", "gllvm", "sjSDM", "Hmsc")),
    target = recode(target, sampled_site_recovery = "Sampled sites\nHidden conditions\npartly learned",
                   new_site_prediction = "New sites\nHidden conditions\naveraged over"))
plot <- ggplot(cells, aes(truth, estimate)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = .4) +
  geom_point(alpha = .16, size = .7, colour = "#007A87") +
  facet_grid(target ~ package) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = c(0,.5,1)) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = c(0,.5,1)) +
  coord_equal() + theme_bw(base_size = 11) +
  labs(x = "True probability for the stated question", y = "Estimated probability",
       title = "Four models, the same simulated community",
       subtitle = "100 training sites, 300 independent test sites, 10 species; no observation error",
       caption = "A point is one species at one site. The diagonal means exact recovery.\nOne community does not rank packages. sjSDM retains some optimisation uncertainty; gllvm uses the checked VA fit.") +
  theme(strip.text.y = element_text(angle = 0), plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "probability-recovery.png"), plot, width = 13, height = 7, dpi = 160)

errors <- read.csv(file.path(out, "overall-errors.csv")) |>
  mutate(package = factor(package, levels = c("occJSDM", "gllvm", "sjSDM", "Hmsc")),
         target = recode(target, sampled_site_recovery = "Reconstruct sampled sites", new_site_prediction = "Predict new sites"))
plot <- ggplot(errors, aes(package, absolute_error_pp, fill = package)) +
  geom_col(width = .65) +
  geom_text(aes(label = sprintf("%.1f", absolute_error_pp)), vjust = -.45, size = 4) +
  facet_wrap(~target) + scale_y_continuous(expand = expansion(mult = c(0,.13))) +
  scale_fill_manual(values = c("#007A87", "#4F5D95", "#CB6E17", "#648747"), guide = "none") +
  theme_bw(base_size = 12) + labs(x = NULL, y = "Average absolute error (percentage points)",
    title = "How far from the true probability, ignoring the direction of error?",
    subtitle = "Lower bars mean smaller errors in this one simulated community",
    caption = "These are different prediction questions, not a train-versus-test overfitting comparison.\nBars show actual saved-fit results. sjSDM remains provisional pending optimisation checks.\nThere is no uncertainty interval across independently simulated communities.") +
  theme(plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "absolute-errors.png"), plot, width = 10, height = 5.7, dpi = 160)

# Preserve the repeated-start differences independently of generating truth.
for (prefix in c("sjSDM-attempt-2", "gllvm-VA")) {
  fs <- files[startsWith(basename(files), prefix)]
  s <- lapply(fs, readRDS)
  predictions <- simplify2array(lapply(s, `[[`, "grid_prediction"))
  cat(prefix, "maximum fixed-grid spread (pp):", 100*max(apply(predictions,1:2,function(x) diff(range(x)))), "\n")
}
cat("Compact results and figures written to", out, "\n")
