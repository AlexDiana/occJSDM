# Compares the summaries agree_small.R wrote: baseline, a second baseline run on other chain seeds
# (the yardstick: how much two runs of the same sampler differ), and the prototype.
#   Rscript agree_small_compare.R <base.rds> <base2.rds> <proto.rds>
a <- commandArgs(TRUE); b <- readRDS(a[1]); b2 <- readRDS(a[2]); p <- readRDS(a[3])
cmp <- function(x, y) {
  z <- (y$mean - x$mean) / sqrt(x$sd^2 / x$ess + y$sd^2 / y$ess)
  d <- (y$mean - x$mean) / sqrt((x$sd^2 + y$sd^2) / 2)
  c(rms_z = sqrt(mean(z^2)), max_z = max(abs(z)), rms_diff_sd = sqrt(mean(d^2)), max_diff_sd = max(abs(d)))
}
out <- do.call(rbind, lapply(c("B_meas", "B_emb", "spread", "log_lambda", "log_sigma_emb", "log_sigma_b", "B0"), function(nm)
  data.frame(block = nm, pair = c("base vs base2", "base vs proto", "base2 vs proto"),
             signif(rbind(cmp(b[[nm]], b2[[nm]]), cmp(b[[nm]], p[[nm]]), cmp(b2[[nm]], p[[nm]])), 3))))
print(out, row.names = FALSE)
blocks <- c("B_meas", "B_emb", "spread", "log_lambda", "log_sigma_emb", "log_sigma_b", "B0")
worst <- function(fit) max(unlist(lapply(fit[blocks], function(x) x$rhat)))
cat(sprintf("\nlargest R-hat: base %.3f, base2 %.3f, proto %.3f\n", worst(b), worst(b2), worst(p)))
if (length(a) > 3) utils::write.csv(out, a[4], row.names = FALSE)
