# Tier 1: default colours of the categorical plotting functions.
#
# Every plot that tells groups apart by colour (primers, TP/FP, chains) uses
# the Okabe-Ito palette from okabe_ito(), so that the default figures are
# colour-blind safe. The lessons apply the same mapping by hand, so the exact
# hex values and their order are a contract, not a styling detail.
#
# Colours are read back with ggplot_build(), which is what the device draws,
# rather than from the scale call in the source.

OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
               "#56B4E9", "#D55E00", "#F0E442", "#000000")

# Map each break of the built colour scale to the colour drawn for it.
built_colour_map <- function(p, breaks) {
  b <- ggplot2::ggplot_build(p)
  sc <- b$plot$scales$get_scales("colour")
  stats::setNames(sc$map(breaks), breaks)
}

# Every colour actually drawn, across all layers.
built_colours <- function(p) {
  b <- ggplot2::ggplot_build(p)
  unique(unlist(lapply(b$data, function(d) d$colour)))
}

# The two-stage fixture with primer names deliberately out of alphabetical
# order, so that a mapping by sorted factor level would be caught.
fixture_twostage_named_primers <- function() {
  fit <- fixture_twostage()
  fit$infos$primerNames <- c("primer_Z", "primer_A")
  fit
}

test_that("okabe_ito() returns the palette in order and recycles with a warning", {
  expect_identical(okabe_ito(3), OKABE_ITO[1:3])
  expect_identical(okabe_ito(8), OKABE_ITO)
  expect_warning(pal <- okabe_ito(10), "recycled")
  expect_identical(pal, OKABE_ITO[c(1:8, 1:2)])
})

test_that("plotFPTPStage2Rates() draws TP in blue and FP in vermillion", {
  p <- plotFPTPStage2Rates(fixture_twostage())

  cols <- built_colour_map(p, c("TP rate", "FP rate"))
  expect_identical(unname(cols), c("#0072B2", "#D55E00"))
  expect_setequal(built_colours(p), c("#0072B2", "#D55E00"))
  # The legend title is unchanged.
  expect_identical(p$scales$get_scales("colour")$name, "Colour")
})

test_that("plotDetectionRates() colours primers in their stored order", {
  fit <- fixture_twostage_named_primers()
  p <- plotDetectionRates(fit)

  cols <- built_colour_map(p, fit$infos$primerNames)
  expect_identical(unname(cols), OKABE_ITO[1:2])
  expect_setequal(built_colours(p), OKABE_ITO[1:2])
  expect_identical(p$labels$colour, "Primer")
})

test_that("plotStage2FPRates() colours primers in their stored order", {
  fit <- fixture_twostage_named_primers()
  p <- plotStage2FPRates(fit)

  cols <- built_colour_map(p, fit$infos$primerNames)
  expect_identical(unname(cols), OKABE_ITO[1:2])
  expect_setequal(built_colours(p), OKABE_ITO[1:2])
  expect_identical(p$labels$colour, "Primer")
})

test_that("plotTraceplot() colours chains with the Okabe-Ito palette", {
  set.seed(1)
  arr <- array(stats::rnorm(2 * 1 * 10 * 3), dim = c(2, 1, 10, 3))
  p <- plotTraceplot(arr)

  cols <- built_colour_map(p, c("1", "2", "3"))
  expect_identical(unname(cols), OKABE_ITO[1:3])
  expect_setequal(built_colours(p), OKABE_ITO[1:3])
  expect_identical(p$labels$colour, "Chain")
})
