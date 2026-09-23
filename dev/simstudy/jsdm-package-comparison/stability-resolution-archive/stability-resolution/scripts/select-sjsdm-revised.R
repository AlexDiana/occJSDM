# Record the revised sjSDM selection from the multistart assessment, BEFORE any
# truth or test outcome is read. Builds a separately versioned results root
# ("revised/") that shares the original inputs, truth and parameter files by
# symbolic link, so the original selection and exports stay untouched.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1])
check <- readRDS(file.path(root, "stability-resolution/multistart/multistart-check.rds"))
stopifnot(isTRUE(check$selection$rule_satisfied))
start <- check$selection$selected_start
label <- sprintf("sjSDM-multistart-start-%d", start)
revised <- file.path(root, "revised")
stopifnot(!dir.exists(revised))
dir.create(file.path(revised, "results"), recursive = TRUE)
dir.create(file.path(revised, "checks"))
for (d in c("inputs", "truth")) file.symlink(file.path(root, d), file.path(revised, d))
original <- readRDS(file.path(root, "results/selection.rds"))
for (package in c("occJSDM", "Hmsc", "gllvm")) {
  f <- paste0(original$selected[[package]], "-parameters.rds")
  file.symlink(file.path(root, "checks", f), file.path(revised, "checks", f))
}
# The original sjSDM parameter files are plain lists with beta and loading;
# save the revised native endpoint in the same shape.
template <- readRDS(file.path(root, "checks", paste0(original$selected[["sjSDM"]], "-parameters.rds")))
stopifnot(identical(names(template), c("beta", "loading")))
fit <- readRDS(file.path(root, "stability-resolution/multistart", paste0(label, ".rds")))
record <- readRDS(file.path(root, "stability-resolution/multistart", paste0(label, "-record.rds")))
stopifnot(isTRUE(record$ok), identical(dim(fit$parameters$beta), dim(template$beta)),
          identical(dim(fit$parameters$loading), dim(template$loading)))
saveRDS(list(beta = fit$parameters$beta, loading = fit$parameters$loading),
        file.path(revised, "checks", paste0(label, "-parameters.rds")))
selected <- original$selected
selected[["sjSDM"]] <- label
selection <- list(selected = selected, original_selected = original$selected,
  recorded_at = Sys.time(), truth_used = FALSE,
  sjSDM_configuration = c(record$settings[c("epochs", "learning_rate", "sampling", "batch",
    "optimiser", "weight_decay", "polish_epochs", "polish_learning_rate", "fit_seed", "polish_seed")],
    backend = record$settings$backend),
  criteria = c(original$criteria[1:2],
    paste("sjSDM: default weight decay 0.0001 arm; twelve independent native starts each with the",
          "1,000-epoch low-step continuation; endpoints classified by deterministic refinement into",
          "the two verified local maxima; the selected fit has the highest precise penalised training",
          "objective, its basin was reached by four independent starts, and within that basin the",
          "native objectives agree within 0.1 and fixed-grid predictions within 1 percentage point")),
  multistart_assessment = check$across, multistart_selection = check$selection,
  basins = check$basins)
saveRDS(selection, file.path(revised, "results/selection.rds"))
write.table(data.frame(package = names(selected), fit = unname(selected)),
            file.path(revised, "results/selected-fits.tsv"), sep = "\t", row.names = FALSE)
print(selection[c("selected", "original_selected", "truth_used", "sjSDM_configuration")])
cat("Revised selection recorded without reading truth.\n")
