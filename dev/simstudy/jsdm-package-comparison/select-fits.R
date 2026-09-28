# Selection uses training fit checks only; this script never reads truth.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
files <- list.files(file.path(root, "checks"), "summary[.]rds$", full.names = TRUE)
summaries <- setNames(lapply(files, readRDS), sub("-summary[.]rds$", "", basename(files)))
stopifnot(all(paste0("sjSDM-attempt-2-start-", 1:3) %in% names(summaries)))
stopifnot(!file.exists(file.path(root, "results", "selection.rds")))
selected <- c(occJSDM = "occJSDM-attempt-1-start-1", Hmsc = "Hmsc-attempt-1-start-1")
stopifnot(all(vapply(summaries[selected], function(x) isTRUE(x$passed), logical(1))))
g <- summaries[startsWith(names(summaries), "gllvm-VA-")]
g <- g[vapply(g, function(s) isTRUE(s$convergence) && s$max_gradient < 0.01, logical(1))]
selected["gllvm"] <- names(g)[which.max(vapply(g, `[[`, numeric(1), "native_loglik"))]
s <- summaries[startsWith(names(summaries), "sjSDM-")]
stopifnot(all(vapply(s, function(a) a$joint_difference < .001 &&
                      a$conditional_difference < 1e-4, logical(1))))
selected["sjSDM"] <- names(s)[which.max(vapply(s, `[[`, numeric(1), "checked_loglik"))]
selection <- list(selected = selected, recorded_at = Sys.time(),
  criteria = c("Bayesian: all recorded Rhat and bulk/tail ESS checks pass",
               "gllvm: best VA objective among native-converged fits with max absolute gradient < .01; keep EVA failures separate",
               "sjSDM: best numerically checked integrated training log likelihood among all recorded starts"),
  truth_used = FALSE)
saveRDS(selection, file.path(root, "results", "selection.rds"))
write.table(data.frame(package = names(selected), fit = unname(selected)),
            file.path(root, "results", "selected-fits.tsv"), sep = "\t", row.names = FALSE)
print(selection)
