# Deterministic posterior draws isolate output arithmetic from sampler behaviour.
# Expected responses below are direct evaluations of intercept + full design.
covariate_response_fixture <- function(raw = data.frame(
    temperature = c(10, 20, 30, 40, 100, 130),
    moisture = c(2, 4, 8, 16, 32, 64),
    habitat = factor(c("dry", "wet", "wood", "dry", "wet", "wood"),
                     levels = c("dry", "wet", "wood"))),
    species = c("species_1", "species_2", "species_3"),
    niter = 4L, nchain = 2L, spline_vars = FALSE, model = "binary") {
  design <- occJSDM:::create_covariates_matrix(raw, spline_vars = spline_vars)
  P <- ncol(design$X); S <- length(species); D <- niter * nchain
  beta <- array(0, c(P, S, niter, nchain))
  intercept <- array(0, c(S, niter, nchain))
  for (ch in seq_len(nchain)) for (it in seq_len(niter)) {
    draw <- (ch - 1L) * niter + it
    for (sp in seq_len(S)) {
      intercept[sp, it, ch] <- c(-2, .5, 2)[sp] + draw / 10
      beta[, sp, it, ch] <- seq_len(P) / 4 * (-1)^seq_len(P) +
        sp / 10 + draw / 30
    }
  }
  list(
    infos = list(S = S, speciesNames = species, jsdmModel = model,
                 X0_psi = design$X0, list_X_psi_mat = design$list_matrix),
    X_psi = design$X,
    results_output = list(jsdm_output = list(B_output = beta, B0_output = intercept))
  )
}

# Draw order is explicitly iteration-within-chain, matching stored MCMC arrays.
response_draws_direct <- function(fit, design, species) {
  j <- fit$results_output$jsdm_output
  ans <- matrix(NA_real_, nrow(design), dim(j$B0_output)[2] * dim(j$B0_output)[3])
  for (ch in seq_len(dim(j$B0_output)[3])) for (it in seq_len(dim(j$B0_output)[2])) {
    draw <- (ch - 1L) * dim(j$B0_output)[2] + it
    eta <- j$B0_output[species, it, ch] +
      as.vector(design %*% j$B_output[, species, it, ch])
    ans[, draw] <- if (fit$infos$jsdmModel == "continuous") eta else plogis(eta)
  }
  ans
}

response_interval_direct <- function(draws, confidence) {
  t(apply(draws, 1, quantile, probs = c(.5, (1-confidence)/2, (1+confidence)/2),
          names = FALSE))
}

test_that("numeric responses use original units, all predictors and the full link", {
  fit <- covariate_response_fixture()
  raw <- data.frame(temperature = c(10,20,30,40,100,130),
                    moisture = c(2,4,8,16,32,64))
  x <- seq(10, 130, length.out = 200)
  # Hold other numeric variables at their median, and habitat at its first level.
  design <- cbind((x - mean(raw$temperature)) / sd(raw$temperature),
                  (median(raw$moisture) - mean(raw$moisture)) / sd(raw$moisture),
                  0, 0)
  result <- returnCovariateEffect(fit, "temperature", idx_species = c(3, 1), confidence = .8)
  expect_identical(unique(result$Species), c("species_3", "species_1"))
  for (s in c(3, 1)) {
    rows <- result[result$Species == fit$infos$speciesNames[s], ]
    expected <- response_interval_direct(response_draws_direct(fit, design, s), .8)
    expect_equal(rows$x, x)
    expect_equal(as.matrix(rows[, c("mean", "lower", "upper")]),
                 expected, ignore_attr = TRUE)
  }
  expect_true(all(unlist(result[, c("mean", "lower", "upper")]) >= 0))
  expect_true(all(unlist(result[, c("mean", "lower", "upper")]) <= 1))
  expect_equal(returnCovariateEffect(fit, "temperature")$Species,
               rep(fit$infos$speciesNames, each = 200))
})

test_that("categorical responses include every level and the intercept", {
  fit <- covariate_response_fixture()
  x <- c(10,20,30,40,100,130); x2 <- c(2,4,8,16,32,64)
  design <- cbind(rep((median(x) - mean(x)) / sd(x), 3),
                  rep((median(x2) - mean(x2)) / sd(x2), 3),
                  c(0, 1, 0), c(0, 0, 1))
  result <- returnCovariateEffect(fit, "habitat", idx_species = c(3, 1))
  expect_identical(levels(result$x), c("dry", "wet", "wood"))
  expect_setequal(as.character(result$x), c("dry", "wet", "wood"))
  expect_equal(nrow(result), 3L * 8L * 2L)
  for (s in c(3, 1)) {
    expected <- response_draws_direct(fit, design, s)
    rows <- result[result$Species == fit$infos$speciesNames[s], ]
    for (k in seq_len(3)) {
      actual <- rows[as.character(rows$x) == c("dry", "wet", "wood")[k], ]
      expect_equal(actual$value, unname(expected[k, ]))
    }
  }
})

test_that("requested confidence controls numeric and categorical plot intervals", {
  fit <- covariate_response_fixture()
  narrow <- returnCovariateEffect(fit, "temperature", 1, confidence = .5)
  wide <- returnCovariateEffect(fit, "temperature", 1, confidence = .95)
  expect_true(all(narrow$lower > wide$lower))
  expect_true(all(narrow$upper < wide$upper))
  plots <- plotCovariateEffect(fit, c("temperature", "habitat"), c(3, 1), confidence = .8)
  expect_named(plots, c("temperature", "habitat"))
  expect_equal(plots$temperature$data,
               returnCovariateEffect(fit, "temperature", c(3,1), confidence = .8))
  expect_equal(plots$temperature$labels$y, "Occupancy probability")
  expect_equal(plots$habitat$labels$y, "Occupancy probability")
  draws <- returnCovariateEffect(fit, "habitat", c(3,1))
  # Check the actual layer: categorical bars must be the requested credible interval,
  # not the default quartiles/whiskers of a boxplot.
  layer <- ggplot2::ggplot_build(plots$habitat)$data[[1]]
  expected <- do.call(rbind, lapply(sort(unique(draws$Species)), function(sp) {
    do.call(rbind, lapply(levels(draws$x), function(level) {
      q <- quantile(draws$value[draws$Species == sp & draws$x == level], c(.1,.9))
      data.frame(ymin = unname(q[1]), ymax = unname(q[2]))
    }))
  }))
  layer <- layer[order(layer$PANEL, layer$x), ]
  expect_equal(layer$ymin, expected$ymin)
  expect_equal(layer$ymax, expected$ymax)
})

test_that("one species, one predictor and one retained draw retain dimensions", {
  raw <- data.frame(x = c(10,20,30,80))
  fit <- covariate_response_fixture(raw, species = "only", niter = 1L, nchain = 1L)
  result <- returnCovariateEffect(fit, "x")
  design <- matrix((seq(10,80,length.out=200) - mean(raw$x))/sd(raw$x), ncol=1)
  truth <- response_draws_direct(fit, design, 1)[,1]
  expect_equal(result$mean, truth)
  expect_equal(result$lower, truth)
  expect_equal(result$upper, truth)
  expect_no_error(plotCovariateEffect(fit, "x"))
})

test_that("continuous response curves stay on the identity scale", {
  fit <- covariate_response_fixture(model = "continuous")
  result <- returnCovariateEffect(fit, "temperature", 1, confidence = .8)
  x <- c(10,20,30,40,100,130); x2 <- c(2,4,8,16,32,64)
  design <- cbind((seq(10,130,length.out=200)-mean(x))/sd(x),
                  (median(x2)-mean(x2))/sd(x2), 0, 0)
  expect_equal(as.matrix(result[, c("mean","lower","upper")]),
               response_interval_direct(response_draws_direct(fit,design,1),.8),
               ignore_attr=TRUE)
  expect_true(any(result$mean < 0))
  expect_equal(plotCovariateEffect(fit,"temperature",1)[[1]]$labels$y,"Expected response")
})

test_that("response splines use training knots and coherent reference covariates", {
  raw <- data.frame(temperature = seq(10,80,length.out=60),
                    moisture = exp(seq(-2,2,length.out=60)))
  fit <- covariate_response_fixture(raw, spline_vars=TRUE)
  x <- seq(10,80,length.out=200)
  standardized <- scale(raw)
  # Independent training spline objects, evaluated at the new raw-scale grid.
  temp_basis <- splines::bs(standardized[,1], df=5, intercept=FALSE)
  moist_basis <- splines::bs(standardized[,2], df=5, intercept=FALSE)
  design <- cbind(predict(temp_basis,(x-mean(raw$temperature))/sd(raw$temperature)),
                  predict(moist_basis,rep((median(raw$moisture)-mean(raw$moisture))/sd(raw$moisture),200)))
  result <- returnCovariateEffect(fit,"temperature",3,confidence=.8)
  expect_equal(result$x,x)
  expect_equal(as.matrix(result[,c("mean","lower","upper")]),
               response_interval_direct(response_draws_direct(fit,design,3),.8),
               ignore_attr=TRUE)
  expect_no_warning(ggplot2::ggplot_build(plotCovariateEffect(fit,"temperature",3)[[1]]))
})

test_that("invalid response requests give useful errors", {
  fit <- covariate_response_fixture()
  for (conf in list(NA_real_,0,1,-.1,1.1,c(.8,.9),"wide")) {
    expect_error(returnCovariateEffect(fit,"temperature",confidence=conf),"confidence")
    expect_error(plotCovariateEffect(fit,"temperature",confidence=conf),"confidence")
  }
  expect_error(returnCovariateEffect(fit,"missing"),"covariate|Covariate")
  expect_error(returnCovariateEffect(fit,c("temperature","moisture")),"one|single")
  for (idx in list(0,-1,4,NA_integer_,1.5)) {
    expect_error(returnCovariateEffect(fit,"temperature",idx),"idx_species")
  }
})

# Prefix/regexp matching must never select another predictor's coefficient.
test_that("similarly named covariates are matched exactly", {
  fit <- covariate_response_fixture(raw = data.frame(x = c(1,2,4,8), x2 = c(3,9,15,21)))
  x <- seq(1,8,length.out=200)
  design <- cbind((x-mean(c(1,2,4,8)))/sd(c(1,2,4,8)), 0)
  result <- returnCovariateEffect(fit,"x",3)
  expect_equal(as.matrix(result[,c("mean","lower","upper")]),
               response_interval_direct(response_draws_direct(fit,design,3),.95),
               ignore_attr=TRUE)
})

test_that("categorical curves retain fitted contrasts, even if options change", {
  raw <- data.frame(habitat = ordered(rep(c("dry", "wet", "wood"), 2),
                                     levels = c("dry", "wet", "wood")),
                    temperature = c(10, 20, 30, 40, 100, 130))
  fits <- list(covariate_response_fixture(raw))
  raw$habitat <- factor(raw$habitat, ordered = FALSE)
  contrasts(raw$habitat) <- contr.sum(3)
  fits[[2]] <- covariate_response_fixture(raw)
  old <- options(contrasts = c("contr.helmert", "contr.helmert"))
  on.exit(options(old))
  for (fit in fits) {
    design <- fit$X_psi[1:3, , drop = FALSE]
    design[, 3] <- (median(raw$temperature) - mean(raw$temperature))/sd(raw$temperature)
    expected <- response_draws_direct(fit, design, 2)
    result <- returnCovariateEffect(fit, "habitat", 2)
    expect_equal(result$value, as.vector(t(expected)))
    # The held-at-first-level contribution must use the same fitted contrasts.
    design <- fit$X_psi[rep(1, 200), , drop = FALSE]
    design[, 3] <- (seq(10,130,length.out=200) - mean(raw$temperature))/sd(raw$temperature)
    result <- returnCovariateEffect(fit, "temperature", 2)
    expect_equal(result$mean, response_interval_direct(response_draws_direct(fit,design,2),.95)[,1])
  }
})

test_that("unused categorical levels are included when their encoding is recoverable", {
  raw <- data.frame(habitat = factor(c("wet", "wood", "wet", "wood"),
                                     levels = c("dry", "wet", "wood")))
  fit <- covariate_response_fixture(raw)
  result <- returnCovariateEffect(fit, "habitat", 1)
  expect_equal(result$value,
               as.vector(t(response_draws_direct(fit, rbind(c(0,0),c(1,0),c(0,1)), 1))))
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old))
  expect_error(returnCovariateEffect(fit, "habitat", 1), "contrast settings")
})

test_that("custom factor encoding can retain one column per level", {
  raw <- data.frame(habitat = factor(rep(c("dry", "wet", "wood"), 2)),
                    temperature = c(10,20,30,40,100,130))
  contrasts(raw$habitat, how.many = 3) <- diag(3)
  fit <- covariate_response_fixture(raw)
  design <- fit$X_psi[1:3, , drop = FALSE]
  design[, 4] <- (median(raw$temperature) - mean(raw$temperature))/sd(raw$temperature)
  result <- returnCovariateEffect(fit, "habitat", 2)
  expect_equal(result$value, as.vector(t(response_draws_direct(fit, design, 2))))
  design <- fit$X_psi[rep(1,200), , drop = FALSE]
  design[, 4] <- (seq(10,130,length.out=200) - mean(raw$temperature))/sd(raw$temperature)
  result <- returnCovariateEffect(fit, "temperature", 2)
  expect_equal(result$mean, response_interval_direct(response_draws_direct(fit,design,2),.95)[,1])
})
