args <- commandArgs(trailingOnly = TRUE)
code <- if (length(args)) args[1] else "dev/simstudy/jsdm-package-comparison/extension"
source(file.path(dirname(code), "pilot-math.R"))
source(file.path(code, "math.R"))
if (file.exists(file.path(code, "calibration-math.R"))) source(file.path(code, "calibration-math.R"))
stopifnot(exists("calibration_summary", mode = "function"))

# Unequal numbers of elements must not turn species/sites into independent replicates.
d <- data.frame(replicate = c(1, 2, 2, 2), estimate = c(2, 0, 0, 0),
                truth = 0, lower = c(1, -1, -1, -1), upper = c(3, 1, 1, 1))
s <- calibration_summary(d)
stopifnot(s$communities == 2, s$elements == 4, s$intervals == 4,
          s$bias == 1, abs(s$rmse - sqrt(2)) < 1e-12,
          s$coverage == .5, s$coverage_mcse == .5, s$width == 2)
# Missing or invalid uncertainty is never silently counted as a miss or zero width.
d$lower[2:4] <- NA_real_
s <- calibration_summary(d)
stopifnot(s$intervals == 1, s$interval_communities == 1, s$coverage == 0,
          is.na(s$coverage_mcse))
d$lower[1] <- 4
s <- calibration_summary(d)
stopifnot(s$intervals == 0, is.na(s$coverage), is.na(s$width), s$bias == 1)

# The passed-only sensitivity analysis must actually exclude flagged fits.
d <- data.frame(replicate = c(1, 2), estimate = c(2, 0), truth = 0,
                lower = c(1, -1), upper = c(3, 1), diagnostic_pass = c(TRUE, FALSE),
                package = "example", scenario = "baseline", n_sites = 100,
                n_species = 1, response = "linear", use_traits = FALSE,
                target = "Marginal probability", term = "grid")
g <- calibration_group_summary(d)
stopifnot(g$bias[g$population == "All scored fits"] == 1,
          g$bias[g$population == "Passed diagnostics only"] == 2,
          g$coverage[g$population == "Passed diagnostics only"] == 0,
          g$communities[g$population == "Passed diagnostics only"] == 1)

# Undo centring jointly with the slopes, rather than rescaling the intercept alone.
b <- array(c(10, 4, 6, 12, 8, 9), c(3, 1, 2, 1))
r <- raw_coefficients(b, c(1, 2), c(2, 3))
stopifnot(identical(dim(r), dim(b)),
          max(abs(as.vector(r) - c(4, 2, 2, 2, 4, 3))) < 1e-12)
stopifnot(max(abs(raw_coefficients(matrix(c(10, 4, 6), 3, 1), c(1,2), c(2,3)) - c(4,2,2))) < 1e-12)

# With no hidden variation, posterior probability bounds are known independently.
p <- list(beta = array(c(-2, -1, 0, 1, 2), c(1,1,5,1)),
          loading = array(0, c(2,1,5,1)), link = "logit")
x <- data.frame(row.names = "origin")
v <- probability_intervals(p, x)
stopifnot(abs(v$estimate - .5) < 1e-12,
          abs(v$lower - (.9 * plogis(-2) + .1 * plogis(-1))) < 1e-12,
          abs(v$upper - (.9 * plogis(2) + .1 * plogis(1))) < 1e-12)
# A probit-normal integral is analytic, including the residual variance.
p$link <- "probit"; p$beta[] <- 1
p$loading[1,1,,] <- 2
v <- probability_intervals(p, x)
stopifnot(abs(v$estimate - pnorm(1/sqrt(5))) < 1e-12,
          abs(v$lower - v$upper) < 1e-12)
cat("Community weighting, missing uncertainty, coefficient scaling and posterior interval identities passed.\n")
