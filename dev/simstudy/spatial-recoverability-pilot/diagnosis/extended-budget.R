# The extension phase of the truth 1 re-run for one job, run on its own: the blockwise study
# sampler through run_joint_oracle at the frozen extended budget with the job's own extension
# seed, screened with the study's screen. Compares with the stored joint-sampler job if present.
# Run from the worktree root:
#   run-r -f dev/simstudy/spatial-recoverability-pilot/diagnosis/extended-budget.R --args JOB
job <- commandArgs(TRUE)[1]
PILOT <- "dev/simstudy/spatial-recoverability-pilot"
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R")) source(file.path(PILOT, f))
compile_pilot_samplers(".")
S0 <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"   # the study directory, as in compare.R
budget <- readRDS(file.path(S0, "config.rds"))$budget
cm <- readRDS(file.path(S0, "inputs", paste0(sub("-(uniform|clustered)$", "", job), ".rds"))); arr <- sub(".*-(uniform|clustered)$", "\\1", job)
idx <- cm$index[[arr]]; n <- length(idx); S <- nrow(cm$species)
t0 <- Sys.time()
o <- run_joint_oracle(cm, arr, nburn = budget$ext_nburn, niter = budget$ext_niter, thin = budget$thin * 4L,
                      seed = pilot_seed(paste0(pilot_label("oracle", job), "-extension")))
secs <- as.numeric(Sys.time() - t0, units = "secs")
sc <- screen_traces(o$traces)
gate_sp <- which(cm$species$occupancy > 0.05)
cat(sprintf("RESULT %s n=%d range=%.2f secs=%.0f pass=%d/36 all_pass=%s gate_species_pass=%d/%d min_bulk=%.0f min_tail=%.0f max_rhat=%.4f median_bulk=%.0f proposals_per_block_step=%.2f\n",
            job, n, cm$range, secs, sum(sc$pass), all(sc$pass), sum(sc$pass[sc$species %in% gate_sp]), sum(sc$species %in% gate_sp),
            min(sc$ess_bulk), min(sc$ess_tail), max(sc$rhat), median(sc$ess_bulk), mean(o$mean_proposals)))
fl <- sc[!sc$pass, ]
for (i in seq_len(nrow(fl))) cat(sprintf("  FAIL %s species %d occupancy=%.2f bulk=%.0f tail=%.0f rhat=%.4f\n", fl$quantity[i], fl$species[i], cm$species$occupancy[fl$species[i]], fl$ess_bulk[i], fl$ess_tail[i], fl$rhat[i]))
stored <- file.path(S0, "jobs", paste0(job, ".rds"))
if (file.exists(stored)) {
  r <- readRDS(stored); w <- joint_weights(cm); ci <- sqrt(1 - cm$shared_fraction); blocks <- 2L + S
  truth <- cm$spatial[idx, , drop = FALSE]; tc <- sweep(truth, 2, colMeans(truth), "-")
  cat(sprintf("  stored job: phases=%s min_bulk=%.0f max_rhat=%.4f median_bulk=%.0f elapsed=%.0f\n", paste(r$phases, collapse = ","), min(r$screen_detail$ess_bulk), max(r$screen_detail$rhat), median(r$screen_detail$ess_bulk), r$elapsed_seconds))
  jm <- sapply(seq_len(blocks), function(k) rowMeans(r$site_draws[[k]])); bm <- sapply(seq_len(blocks), function(k) rowMeans(o$site_draws[[k]]))
  hj <- sapply(seq_len(S), function(s) w[s, 1] * jm[, 1] + w[s, 2] * jm[, 2] + ci * jm[, 2 + s]); hb <- sapply(seq_len(S), function(s) w[s, 1] * bm[, 1] + w[s, 2] * bm[, 2] + ci * bm[, 2 + s])
  sl <- function(h) sapply(gate_sp, function(s) sum((h[, s] - mean(h[, s])) * tc[, s]) / sum(tc[, s]^2)); cr <- function(h) sapply(gate_sp, function(s) cor(h[, s], tc[, s]))
  cat(sprintf("  species effect at sites, stored vs this run: posterior mean correlation min=%.4f median=%.4f; RMS difference median=%.3f against RMS mean median=%.3f\n",
              min(diag(cor(hj, hb))), median(diag(cor(hj, hb))), median(sqrt(colMeans((hj - hb)^2))), median(sqrt(colMeans(hb^2)))))
  cat(sprintf("  training slope (gate species) median stored=%.3f this=%.3f max|diff|=%.3f; training correlation median stored=%.3f this=%.3f max|diff|=%.3f\n",
              median(sl(hj)), median(sl(hb)), max(abs(sl(hj) - sl(hb))), median(cr(hj)), median(cr(hb)), max(abs(cr(hj) - cr(hb)))))
}
cat("DONE", job, "\n")
