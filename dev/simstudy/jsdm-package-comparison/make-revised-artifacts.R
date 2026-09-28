# Compact records and figures for the revised comparison (sjSDM replaced by its
# selected weak-penalty fit). Written beside the stability records; the
# original pilot-results/ directory is left as the historical record.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
revised <- file.path(root, "revised/results")
out <- file.path(code, "stability-resolution-results")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
for (f in c("overall-errors.csv", "errors-by-band.csv", "errors-by-species.csv",
            "integration-checks.csv", "selected-fits.tsv", "sjsdm-revision-comparison.csv",
            "sjsdm-revision-prediction-differences.csv", "sjsdm-revision-species-comparison.csv")) {
  stopifnot(file.copy(file.path(revised, f), file.path(out, paste0("revised-", f)), overwrite = TRUE))
}
cells <- readRDS(file.path(revised, "predictions.rds")) |>
  mutate(package = factor(package, levels = c("occJSDM", "gllvm", "sjSDM", "Hmsc")),
    target = recode(target, sampled_site_recovery = "Sampled sites\nHidden conditions\npartly learned",
                   new_site_prediction = "New sites\nHidden conditions\naveraged over"))
plot <- ggplot(cells, aes(truth, estimate)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = .4) +
  geom_point(alpha = .16, size = .7, colour = "#007A87") +
  facet_grid(target ~ package) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = c(0, .5, 1)) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = c(0, .5, 1)) +
  coord_equal() + theme_bw(base_size = 11) +
  labs(x = "True probability for the stated question", y = "Estimated probability",
       title = "Four models, the same simulated community",
       subtitle = "100 training sites, 300 independent test sites, 10 species; no observation error",
       caption = "A point is one species at one site. The diagonal means exact recovery.\nOne community does not rank packages. sjSDM uses its selected weak-penalty fit; gllvm uses the checked VA fit.") +
  theme(strip.text.y = element_text(angle = 0), plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "revised-probability-recovery.png"), plot, width = 13, height = 7, dpi = 160)
errors <- read.csv(file.path(revised, "overall-errors.csv")) |>
  mutate(package = factor(package, levels = c("occJSDM", "gllvm", "sjSDM", "Hmsc")),
         target = recode(target, sampled_site_recovery = "Reconstruct sampled sites", new_site_prediction = "Predict new sites"))
plot <- ggplot(errors, aes(package, absolute_error_pp, fill = package)) +
  geom_col(width = .65) +
  geom_text(aes(label = sprintf("%.1f", absolute_error_pp)), vjust = -.45, size = 4) +
  facet_wrap(~target) + scale_y_continuous(expand = expansion(mult = c(0, .13))) +
  scale_fill_manual(values = c("#007A87", "#4F5D95", "#CB6E17", "#648747"), guide = "none") +
  theme_bw(base_size = 12) + labs(x = NULL, y = "Average absolute error (percentage points)",
    title = "How far from the true probability, ignoring the direction of error?",
    subtitle = "Lower bars mean smaller errors in this one simulated community",
    caption = "These are different prediction questions, not a train-versus-test overfitting comparison.\nBars show actual saved-fit results with the selected weak-penalty sjSDM fit.\nThere is no uncertainty interval across independently simulated communities.") +
  theme(plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "revised-absolute-errors.png"), plot, width = 10, height = 5.7, dpi = 160)
cat("Revised compact results and figures written to", out, "\n")
