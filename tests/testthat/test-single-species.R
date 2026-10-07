# A one-column state matrix must stay a matrix when expanded from sites to
# field samples. Dropping that dimension breaks sample_theta0() and prevents
# otherwise valid single-species occupancy and two-stage fits.
test_that("field false-positive sampling retains the single species", {
  z <- matrix(c(0, 1, 0, 1), ncol = 1L)
  idx_z <- c(1L, 1L, 2L, 3L, 3L, 4L)
  w <- matrix(c(1, 0, 1, 0, 0, 1), ncol = 1L)

  # Among the four samples from absent sites, one contains DNA and three
  # do not: the Beta(1, 20) prior therefore becomes Beta(2, 23).
  set.seed(7)
  expected <- stats::rbeta(1L, 2, 23)
  set.seed(7)
  expect_equal(occJSDM:::sample_theta0(z, w, idx_z, 1, 20), expected)
})

single_species_survey <- function(model, species) {
  site <- rep(seq_len(12L), each = 3L)
  sample <- seq_along(site)
  if (model == "two_stage") {
    site <- rep(site, each = 2L)
    sample <- rep(sample, each = 2L)
  }
  info <- data.frame(Site = site, Sample = sample, Primer = 1L,
                     X_psi.env = seq(-1, 1, length.out = 12L)[site])
  otu <- vapply(seq_len(species), function(s) {
    as.numeric((seq_along(site) + s) %% 3L == 0L) * 3
  }, numeric(length(site)))
  colnames(otu) <- paste0("species_", seq_len(species))
  list(info = info, OTU = otu)
}

for (design in c("occupancy", "two_stage")) {
  for (species_count in c(1L, 2L)) {
    test_that(paste(design, "fitting supports", species_count, "species without latent factors"), {
      set.seed(23)
      fit <- suppressMessages(suppressWarnings(runOccJSDM(
        single_species_survey(design, species_count),
        occCovariates = "X_psi.env",
        listParams = list(n_factors = 0L, n_lattrait = 0L),
        MCMCparams = list(nchain = 2L, nburn = 10L, niter = 20L, nthin = 1L)
      )))
      expect_identical(fit$infos$model, design)
      expect_equal(fit$infos$S, species_count)
      psi <- fit$results_output$psi_output
      expect_equal(dim(psi), c(12L, species_count))
      expect_true(all(is.finite(psi) & psi >= 0 & psi <= 1))
      theta0 <- fit$results_output$theta0_output
      expect_equal(dim(theta0), c(species_count, 20L, 2L))
      expect_true(all(is.finite(theta0) & theta0 >= 0 & theta0 <= 1))
    })
  }
}
