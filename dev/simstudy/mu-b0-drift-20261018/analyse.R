# Which term of the occupancy linear predictor moves against the species'
# intercepts in the three fits where mu_b0 failed to converge (README.md).
# Run from the repository root:
#   Rscript dev/simstudy/mu-b0-drift-20261018/analyse.R
# Reads draws/main/*.rds (refit.R) and, for the reproduction check, the main
# run's draws in MAIN_DRAWS (default: draws/main of the 100-species check in
# the worktree .claude/worktrees/nostalgic-mendel-8b7449).
# Writes:
#   reproduction.csv     the stored draws against the main run's at every
#                        10th iteration
#   chain-difference.csv per fit: chain 2 minus chain 1 of the mean over
#                        draws of each quantity, its ratio to the B0
#                        difference, and the posterior sd
#   species.csv          per fit and species: the chain difference of B0 and
#                        of the site-mean of each term, of occupied sites and
#                        of psi, with the truth and the detections
#   partners.csv         per fit: across species, the correlation and slope
#                        of each term's difference on the B0 difference;
#                        within chains, the correlation of mean B0 and mu_b0
#                        with each term; R-hat and bulk ESS of each mean
#   basis.csv            per fit and range: how well the spatial basis fits a
#                        constant over sites, and the coefficient norm of a
#                        unit shift; the share of draws at that range
#   traces.rds           per fit, chain and stored draw: every quantity
#   traces.png           batch means per chain of each quantity, every fit
#   species.png          per fit: each species' B0 shift against its spatial,
#                        factor and eta shifts
here <- Sys.getenv("RUN_DIR", "dev/simstudy/mu-b0-drift-20261018")
out <- Sys.getenv("OUT_DIR", here)
main_draws <- Sys.getenv("MAIN_DRAWS", paste0("/Users/douglasyu/src/occJSDM/.claude/worktrees/",
  "nostalgic-mendel-8b7449/dev/simstudy/b0-community-s100-20261018/draws/main"))
source(file.path(here, "terms.R"))
rhat <- posterior::rhat; ess_bulk <- posterior::ess_bulk

files <- sort(list.files(file.path(here, "draws", "main"), "[.]rds$", full.names = TRUE))
stopifnot(length(files) > 0)
fits <- setNames(lapply(files, readRDS), sub("[.]rds$", "", basename(files)))

# --- Reproduction: the stored draws against the main run's ----------------
thin_draws <- function(a, keep) {
  d <- length(dim(a))
  idx <- seq(keep, dim(a)[d - 1], by = keep)
  if (d == 2) a[idx, , drop = FALSE] else if (d == 3) a[, idx, , drop = FALSE] else
    a[, , idx, , drop = FALSE]
}
reproduction <- do.call(rbind, lapply(names(fits), function(key) {
  x <- fits[[key]]
  qs <- c("mu_b0", "B0", "beta_theta0", "theta0", "p")
  mf <- file.path(main_draws, paste0(key, ".rds"))
  if (!file.exists(mf)) {
    return(data.frame(fit = key, quantity = qs, n = NA, max_abs_difference = NA, identical = NA))
  }
  m <- readRDS(mf)$draws
  do.call(rbind, lapply(qs, function(q) {
    a <- thin_draws(m[[q]], x$mcmc$nthin); b <- x$draws[[q]]
    data.frame(fit = key, quantity = q, n = length(b), max_abs_difference = max(abs(a - b)),
               identical = identical(dim(a), dim(b)) && all(a == b))
  }))
}))
write.csv(reproduction, file.path(out, "reproduction.csv"), row.names = FALSE)
cat(sprintf("reproduction: %d of %d blocks identical to the main run's draws (%d not compared)\n",
            sum(reproduction$identical, na.rm = TRUE), nrow(reproduction), sum(is.na(reproduction$identical))))

# --- Per draw: the mean of every term; per species: its site mean ---------
summarise_fit <- function(x) {
  d <- x$draws; des <- x$design
  S <- dim(d$B0)[1]; Tn <- dim(d$B0)[2]; nch <- dim(d$B0)[3]; n <- nrow(des$X_psi)
  nf <- dim(d$U)[2]
  qn <- c("mu_b0", "mean_B0", "mean_XB", "mean_SE", "mean_UL", "mean_eta", "mean_psi",
          "occupied", "mean_bt0", "mean_logit_theta0", "mean_logit_p",
          paste0("mean_U", seq_len(nf)), "sigmah", "sigmabs", "l_s", "eta_check")
  draws <- array(NA_real_, c(Tn, nch, length(qn)), dimnames = list(NULL, NULL, qn))
  spn <- c("B0", "XB", "SE", "UL", "eta", "occ", "psi")
  species <- array(NA_real_, c(S, Tn, nch, length(spn)), dimnames = list(NULL, NULL, NULL, spn))
  for (ch in seq_len(nch)) for (t in seq_len(Tn)) {
    e <- eta_terms(d, des, t, ch)
    psi <- d$psi[, , t, ch]; z <- d$z[, , t, ch]
    species[, t, ch, ] <- cbind(e$B0, colMeans(e$XB), colMeans(e$SE), colMeans(e$UL),
                                colMeans(e$eta), colSums(z), colMeans(psi))
    draws[t, ch, ] <- c(d$mu_b0[t, ch], mean(e$B0), mean(e$XB), mean(e$SE), mean(e$UL),
                        mean(e$eta), mean(psi), sum(z), mean(d$beta_theta0[, t, ch]),
                        mean(qlogis(d$theta0[, t, ch])), mean(qlogis(d$p[, , t, ch])),
                        colMeans(matrix(d$U[, , t, ch], nrow = n)),
                        d$sigmah[t, ch], d$sigmabs[t, ch], des$l_s_grid[d$idx_ls[t, ch]],
                        max(abs(plogis(e$eta) - psi)))
  }
  list(draws = draws, species = species)
}
summaries <- lapply(fits, summarise_fit)
eta_check <- vapply(summaries, function(sm) max(sm$draws[, , "eta_check"]), 0)
cat("stored psi against logistic(eta) rebuilt from the stored terms, max |difference| per fit:\n")
print(signif(eta_check, 3))
stopifnot(all(eta_check < 1e-10))

# --- Chain 2 minus chain 1 ---------------------------------------------------
chain_difference <- do.call(rbind, lapply(names(fits), function(key) {
  dr <- summaries[[key]]$draws
  m1 <- colMeans(dr[, 1, ]); m2 <- colMeans(dr[, 2, ])
  dB0 <- unname(m2["mean_B0"] - m1["mean_B0"])
  data.frame(fit = key, quantity = names(m1), chain1 = unname(m1), chain2 = unname(m2),
             difference = unname(m2 - m1), ratio_to_B0 = unname(m2 - m1) / dB0,
             posterior_sd = apply(dr, 3, sd), row.names = NULL)
}))
chain_difference <- chain_difference[chain_difference$quantity != "eta_check", ]
write.csv(chain_difference, file.path(out, "chain-difference.csv"), row.names = FALSE)

# --- Species by species --------------------------------------------------------
species <- do.call(rbind, lapply(names(fits), function(key) {
  x <- fits[[key]]; s <- summaries[[key]]$species
  m <- function(ch) apply(s[, , ch, ], c(1, 3), mean)
  dd <- m(2) - m(1)
  sds <- apply(s, c(1, 4), sd)
  detected <- vapply(seq_len(ncol(x$OTU)), function(j)
    sum(tapply(x$OTU[, j] > 0, x$info$Site, any)), 0)
  data.frame(fit = key, species = seq_len(nrow(dd)),
             true_B0 = x$truth$jsdmParams_true$B0, occ_true = colMeans(x$truth$z_true),
             detected_sites = detected,
             d_B0 = dd[, "B0"], d_XB = dd[, "XB"], d_SE = dd[, "SE"], d_UL = dd[, "UL"],
             d_eta = dd[, "eta"], d_occupied = dd[, "occ"], d_psi = dd[, "psi"],
             sd_B0 = sds[, "B0"], sd_SE = sds[, "SE"], sd_UL = sds[, "UL"], sd_eta = sds[, "eta"],
             row.names = NULL)
}))
write.csv(species, file.path(out, "species.csv"), row.names = FALSE)

# --- Partners -------------------------------------------------------------------
partners <- do.call(rbind, lapply(names(fits), function(key) {
  dr <- summaries[[key]]$draws
  sp <- species[species$fit == key, ]
  row <- function(statistic, term, value) data.frame(fit = key, statistic = statistic, term = term,
                                                     value = unname(value))
  rows <- list()
  for (tm in c("XB", "SE", "UL", "eta", "occupied", "psi")) {
    y <- sp[[paste0("d_", tm)]]
    rows[[length(rows) + 1]] <- row("across_species_cor_with_d_B0", tm, cor(sp$d_B0, y))
    rows[[length(rows) + 1]] <- row("across_species_slope_on_d_B0", tm, coef(lm(y ~ sp$d_B0))[2])
  }
  for (tm in c("B0", "eta")) {
    rows[[length(rows) + 1]] <- row("across_species_mean_of_d", tm, mean(sp[[paste0("d_", tm)]]))
    rows[[length(rows) + 1]] <- row("across_species_sd_of_d", tm, sd(sp[[paste0("d_", tm)]]))
  }
  for (tm in c("true_B0", "occ_true", "detected_sites")) {
    rows[[length(rows) + 1]] <- row("across_species_cor_of_d_B0_with", tm, cor(sp$d_B0, sp[[tm]]))
  }
  within <- function(a, b) mean(vapply(seq_len(dim(dr)[2]), function(ch) cor(dr[, ch, a], dr[, ch, b]), 0))
  terms <- c("mean_XB", "mean_SE", "mean_UL", "mean_eta", "mean_psi", "occupied", "mean_bt0",
             "mean_logit_theta0", "mean_U1", "mean_U2")
  for (tm in intersect(terms, dimnames(dr)[[3]])) {
    rows[[length(rows) + 1]] <- row("within_chain_cor_with_mean_B0", tm, within("mean_B0", tm))
    rows[[length(rows) + 1]] <- row("within_chain_cor_with_mu_b0", tm, within("mu_b0", tm))
  }
  rows[[length(rows) + 1]] <- row("within_chain_cor_with_mu_b0", "mean_B0", within("mu_b0", "mean_B0"))
  # Species by species within each chain: does a species' own B0 trade off
  # against the site mean of its own spatial or factor term, or against its
  # own occupied sites? The mean over species of the within-chain correlation.
  s4 <- summaries[[key]]$species
  for (tm in c("XB", "SE", "UL", "eta", "occ", "psi")) {
    r <- vapply(seq_len(dim(s4)[1]), function(j) mean(vapply(seq_len(dim(s4)[3]), function(ch)
      cor(s4[j, , ch, "B0"], s4[j, , ch, tm]), 0)), 0)
    rows[[length(rows) + 1]] <- row("within_chain_species_cor_B0_with_own", tm, mean(r, na.rm = TRUE))
  }
  for (q in intersect(c("mu_b0", "mean_B0", terms), dimnames(dr)[[3]])) {
    rows[[length(rows) + 1]] <- row("rhat", q, rhat(dr[, , q]))
    rows[[length(rows) + 1]] <- row("ess_bulk", q, ess_bulk(dr[, , q]))
  }
  do.call(rbind, rows)
}))
write.csv(partners, file.path(out, "partners.csv"), row.names = FALSE)

saveRDS(lapply(summaries, function(sm) signif(sm$draws, 7)), file.path(out, "traces.rds"))

# --- Can the spatial term absorb a constant shift? ----------------------------
# For each range of the grid, the least-squares fit of a constant 1 over the
# sites by the spatial basis Ks (the SoR basis the sampler uses at that
# range): the R-squared, and the squared norm of the coefficients that give a
# unit shift, which is the prior penalty of a unit shift in units of
# 1 / (2 sigma_bs^2). The share of draws at each range says which ranges the
# chains used.
basis <- do.call(rbind, lapply(names(fits), function(key) {
  x <- fits[[key]]; des <- x$design; n <- nrow(des$X_psi)
  used <- table(factor(x$draws$idx_ls, levels = seq_along(des$l_s_grid))) / length(x$draws$idx_ls)
  do.call(rbind, lapply(seq_along(des$l_s_grid), function(j) {
    Ks <- matrix(des$Ks_all[, , j], nrow = n)
    # the basis as it multiplies Bs: site i, support point Xs_centers[i, k]
    m <- ncol(des$Xs_centers); ps <- dim(x$draws$Bs)[1]
    A <- matrix(0, n, ps)
    for (k in seq_len(m)) A[cbind(seq_len(n), des$Xs_centers[, k])] <- A[cbind(seq_len(n), des$Xs_centers[, k])] + Ks[, k]
    fit1 <- lm.fit(A, rep(1, n))
    data.frame(fit = key, range = des$l_s_grid[j], share_of_draws = as.numeric(used[j]),
               r_squared_constant = 1 - sum(fit1$residuals^2) / n,
               coef_norm2_unit_shift = sum(fit1$coefficients^2, na.rm = TRUE))
  }))
}))
write.csv(basis, file.path(out, "basis.csv"), row.names = FALSE)

# --- Plots ----------------------------------------------------------------------
show <- c("mu_b0", "mean_B0", "mean_XB", "mean_SE", "mean_UL", "mean_eta", "mean_psi", "occupied")
png(file.path(out, "traces.png"), width = 400 * length(show), height = 320 * length(fits), res = 110)
par(mfrow = c(length(fits), length(show)), mar = c(3.5, 4, 2.5, 0.5), mgp = c(2.2, 0.7, 0))
for (key in names(fits)) {
  dr <- summaries[[key]]$draws; Tn <- dim(dr)[1]; nthin <- fits[[key]]$mcmc$nthin
  b <- 25; nb <- Tn %/% b
  for (q in show) {
    y <- dr[, , q]
    plot(NA, xlim = c(0, Tn * nthin), ylim = range(y), xlab = "iteration after burn-in", ylab = q,
         main = if (q == show[1]) key else "")
    for (ch in seq_len(ncol(y))) {
      lines(seq_len(Tn) * nthin, y[, ch], col = adjustcolor(c("black", "red")[ch], 0.25))
      lines((seq_len(nb) - 0.5) * b * nthin, colMeans(matrix(y[seq_len(nb * b), ch], nrow = b)),
            col = c("black", "red")[ch], lwd = 2.5)
    }
  }
}
dev.off()
png(file.path(out, "species.png"), width = 1500, height = 420 * length(fits), res = 110)
par(mfrow = c(length(fits), 4), mar = c(4, 4, 2.5, 0.5), mgp = c(2.2, 0.7, 0))
for (key in names(fits)) {
  sp <- species[species$fit == key, ]
  for (tm in c("d_SE", "d_UL", "d_eta", "d_occupied")) {
    plot(sp$d_B0, sp[[tm]], pch = 19, col = adjustcolor("navy", 0.6),
         xlab = "chain 2 minus chain 1 in B0", ylab = tm, main = if (tm == "d_SE") key else "")
    abline(h = 0, v = 0, col = "grey70")
    if (tm != "d_occupied") abline(0, -1, lty = 2, col = "grey40")
  }
}
dev.off()

# --- Print the headline numbers ----------------------------------------------
options(width = 200)
cat("\nchain 2 minus chain 1, mean over sites, species and draws:\n")
wide <- reshape(chain_difference[, c("fit", "quantity", "difference")], idvar = "quantity",
                timevar = "fit", direction = "wide")
print(wide, digits = 3, row.names = FALSE)
cat("\npartners:\n")
pw <- reshape(partners, idvar = c("statistic", "term"), timevar = "fit", direction = "wide")
print(pw, digits = 3, row.names = FALSE)
