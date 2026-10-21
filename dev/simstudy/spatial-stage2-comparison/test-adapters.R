# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-adapters.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("upstream/adapter-contract.R", "upstream/basis.R", "upstream/adapter-occjsdm.R", "upstream/adapter-hmsc.R",
            "upstream/adapter-gllvm.R", "upstream/diagnostics.R", "upstream/select.R", "config.R", "data.R")) source(file.path(D, f))
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
# One community is enough for adapter tests; reading it directly avoids hashing every oracle file.
man <- read.csv(file.path(PILOT, "input-manifest.csv"), stringsAsFactors = FALSE)
row <- list(dataset_id = "t2-n300-rep1", input_file = file.path(PILOT, "inputs", "t2-n300-rep1.rds"), input_md5 = man$md5[man$file == "t2-n300-rep1.rds"])
cm <- readRDS(row$input_file); tr <- load_training(cm, row$dataset_id, row$input_md5); pr <- load_prediction(cm)
# A 40-site, 4-species slice keeps every fit under a minute.
slice <- function(tr, sites = 1:40, species = c(4L, 6L, 9L, 12L)) {
  tr$site_id <- tr$site_id[sites]; tr$species_id <- tr$species_id[species]; tr$y <- tr$y[sites, species, drop = FALSE]
  tr$environment <- tr$environment[sites, , drop = FALSE]; tr$xy <- tr$xy[sites, , drop = FALSE]; tr }
small <- slice(tr); pr_small <- pr; pr_small$species_id <- small$species_id; pr_small <- subset_prediction(pr_small, 1:30)
tiny <- c(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L)
settings <- function(...) list(seed = 11L, job_id = 1L, attempt_id = 1L, source_hashes = c(x = "0"), mcmc = tiny,
                               gllvm = stage2_config()$gllvm, prediction_seed = 12L, ...)

test_that("spatial6 is a valid basis and Hmsc fits six factors", {
  # Hmsc caps the factor count at the species count, so six factors need at least six species.
  six <- slice(tr, species = 4:10); pr_six <- pr; pr_six$species_id <- six$species_id; pr_six <- subset_prediction(pr_six, 1:30)
  b <- fit_basis(six, "spatial6"); expect_identical(b$arm, "spatial6"); expect_equal(dim(predict_basis(b, pr_six)), c(30L, 2L))
  f <- fit_hmsc(six, b, settings(hmsc_method = "Full"))
  expect_equal(nrow(f$native$postList[[1]][[1]]$Lambda[[1]]), 6L)
  p <- predict_hmsc(f, pr_six, b); expect_true(validate_prediction(p, pr_six))
})

test_that("Hmsc GPP route fits and predicts with recorded knots", {
  b <- fit_basis(small, "spatial2")
  f <- fit_hmsc(small, b, settings(hmsc_method = "GPP", hmsc_knots = 10L))
  expect_equal(nrow(f$native$ranLevels$site$sKnot), 10L)
  p <- predict_hmsc(f, pr_small, b); expect_true(validate_prediction(p, pr_small)); expect_true(grepl("GPP", p$conditioning))
})

test_that("occJSDM gp honours the support-point setting and records it", {
  b <- fit_basis(small, "gp")
  f <- fit_occjsdm(small, b, settings(support_points = 20L))
  expect_equal(nrow(f$native$infos$list_Xs$X_tilde), 20L); expect_equal(f$support_points, 20L)
  p <- predict_occjsdm(f, pr_small, b); expect_true(validate_prediction(p, pr_small))
  expect_error(fit_occjsdm(small, b, settings(support_points = 41L)))
})

test_that("gllvm fits a dataset with an all-zero training site", {
  zero <- small; zero$y[1, ] <- 0L   # a site where no species is present
  bg <- fit_basis(zero, "spatial2"); fz <- fit_gllvm(zero, bg, settings())
  expect_identical(fz$starting_val, "random"); expect_true(is.finite(fz$native$objective))
})

test_that("env and gllvm arms still fit, and diagnostics run for each package", {
  fe <- fit_occjsdm(small, fit_basis(small, "env"), settings(support_points = NA)); de <- diagnose_fit(fe, small, fit_basis(small, "env"))
  expect_true(all(c("failed", "reason") %in% names(de)))
  bg <- fit_basis(small, "spatial2"); fg <- fit_gllvm(small, bg, settings()); dg <- diagnose_fit(fg, small, bg)
  expect_equal(nrow(dg), 1L); pg <- predict_gllvm(fg, pr_small, bg); expect_true(validate_prediction(pg, pr_small))
  # Random starts fail without latent variables in gllvm 2.0.15, so the arms outside the study keep residual starts.
  expect_identical(fg$starting_val, "random"); expect_identical(fit_gllvm(small, fit_basis(small, "env"), settings())$starting_val, "res")
})
cat("Adapter tests passed.\n")
