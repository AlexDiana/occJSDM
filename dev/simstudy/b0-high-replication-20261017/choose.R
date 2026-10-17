# Applies README.md's spread criterion and remedy rule to the diagnosis refits, after the diagnosis step and before
# the intervention step. Writes choose.csv and draws/interventions.txt (the variants the intervention step runs).
# Run from the repository root or from this directory.
here <- if (file.exists("refit_spread.R")) "." else "dev/simstudy/b0-high-replication-20261017"
source(file.path(here, "common.R"))
arms <- c("M20", "M20_cm", "base_rare_cm_m20")
d <- species_table(read_fits(file.path(here, "draws", "diag")))
d <- d[d$variant == "base" & d$arm %in% arms, ]

res <- do.call(rbind, lapply(arms, function(k) {
  x <- d[d$arm == k, ]
  sl <- cluster_slope(x$e, x$P, x$replicate)
  a_holds <- unname(sl["slope"] >= 0.5 && sl["slope"] <= 1.5 && abs(sl["slope"]) > 1.96 * sl["se"])
  bands <- lapply(extreme_bands(k), function(b) {
    y <- x[x$band_index == b, ]
    list(band = y$band[1], n = nrow(y), e = mean(y$e), P = mean(y$P))
  })
  b_ok <- vapply(bands, function(v) sign(v$P) == sign(v$e) && abs(v$P) >= 0.5 * abs(v$e), TRUE)
  data.frame(arm = k, replicates = length(unique(x$replicate)), species_fits = nrow(x),
             slope = unname(sl["slope"]), slope_se = unname(sl["se"]), r2 = unname(sl["r2"]),
             band_1 = bands[[1]]$band, n_1 = bands[[1]]$n, e_1 = bands[[1]]$e, P_1 = bands[[1]]$P,
             band_2 = bands[[2]]$band, n_2 = bands[[2]]$n, e_2 = bands[[2]]$e, P_2 = bands[[2]]$P,
             a_holds = a_holds, b_holds = all(b_ok), holds = a_holds && all(b_ok), stringsAsFactors = FALSE)
}))
k_hold <- sum(res$holds)
answer <- if (k_hold >= 2) "yes" else if (k_hold == 1) "partly" else "no"
ext <- d[d$extreme, ]
shortfall_cov <- mean(ext$dvar_XB)
shortfall_site <- mean(ext$dvar_site)
remedy <- if (answer == "no") "" else if (shortfall_site >= shortfall_cov) "d4" else "gt4"
res$answer <- answer
res$shortfall_cov <- shortfall_cov
res$shortfall_site <- shortfall_site
res$extreme_species_fits <- nrow(ext)
res$remedy <- remedy
print(res, row.names = FALSE)
write.csv(res, file.path(here, "choose.csv"), row.names = FALSE)
variants <- c("sb0", if (nzchar(remedy)) remedy)
writeLines(paste(variants, collapse = ","), file.path(here, "draws", "interventions.txt"))
message(sprintf("spread answer: %s; covariate shortfall %.3f, site-level shortfall %.3f; interventions: %s",
                answer, shortfall_cov, shortfall_site, paste(variants, collapse = ",")))
