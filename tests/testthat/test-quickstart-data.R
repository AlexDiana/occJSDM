# The Quickstart's shipped example: a simulated river survey and its fit.

test_that("quickstartdata has the agreed survey design", {
  d <- quickstartdata
  expect_named(d, c("info", "OTU", "traits"))
  expect_equal(nrow(d$info), 1800L)
  expect_equal(length(unique(d$info$Site)), 100L)
  expect_true(all(table(unique(d$info[c("Site", "Sample")])$Site) == 3L))
  expect_equal(sort(unique(d$info$Primer)), 1:2)
  expect_true(all(table(interaction(d$info$Sample, d$info$Primer, drop = TRUE)) == 3L))
  # "not_protected" first, so the coefficient reads "protected compared with not protected".
  expect_identical(levels(d$info$protected), c("not_protected", "protected"))
  expect_true(all(c("elevation", "hmi", "water_volume") %in% names(d$info)))
  expect_equal(nrow(d$OTU), nrow(d$info))
  expect_true(ncol(d$OTU) >= 18L && ncol(d$OTU) <= 22L)
  expect_true(all(grepl("^[A-Z][a-z]+_[A-Z][a-z]+_[0-9]+$", colnames(d$OTU))))
  expect_identical(rownames(d$traits), colnames(d$OTU))
  expect_identical(colnames(d$traits), c("body_mass", "red_listed"))
  expect_true(all(d$traits[, "red_listed"] %in% c(0, 1)))
  expect_true(any(d$traits[, "red_listed"] == 1) && any(d$traits[, "red_listed"] == 0))
})
