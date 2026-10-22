library(testthat)
helper <- "dev/simstudy/vignette-lesson/helpers.R"
if (file.exists(helper)) source(helper)
score <- "dev/simstudy/vignette-lesson/score_lesson.R"
if (file.exists(score)) source(score)

test_that("the lesson has explicit observation and threshold helpers", {
  expect_true(exists("lesson_observations", mode = "function"))
  expect_true(exists("effective_detection_rate", mode = "function"))
})

if (exists("lesson_observations", mode = "function")) {
  test_that("source labels follow site and sample truth despite row order", {
    # A laboratory FP can occur in a sample from an occupied site.
    info <- data.frame(Site=c("b", "a", "a", "b"),
                       Sample=c("s2", "s1", "s1", "s2"),
                       Primer=c(1, 1, 1, 1))
    reads <- matrix(c(8, 1, NA, 0, 2, 9, 0, 1), 4, 2,
                    dimnames=list(NULL, c("sp2", "sp1")))
    z <- matrix(c(1, 0, 1, 1), 2, 2,
                dimnames=list(c("a", "b"), c("sp1", "sp2")))
    w <- matrix(c(1, 1, 0, 1), 2, 2,
                dimnames=list(c("s1", "s2"), c("sp1", "sp2")))
    b <- list(sim=list(data_list=list(info=info, OTU=reads),
                       true_params=list(z_true=z, w_true=w)), threshold=1)
    x <- lesson_observations(b)
    pick <- function(sp, row) x[x$species == sp & x$row == row, ]
    expect_equal(pick("sp2", 2)$source, "Laboratory false positive")
    expect_equal(pick("sp2", 2)$z, 1)
    expect_equal(pick("sp1", 1)$source, "Field-stage false positive")
    expect_equal(pick("sp1", 2)$source, "True detection")
    expect_equal(pick("sp2", 3)$source, "Missing")
    expect_true(is.na(pick("sp2", 3)$positive))
    expect_equal(pick("sp1", 3)$source, "No detection")
    b$sim$data_list$info$Sample[1] <- "unmatched"
    expect_error(lesson_observations(b), "identity")
  })

  test_that("read rounding changes the fitted-scale rate", {
    expect_equal(effective_detection_rate(.2, 1.5, 1, 1),
                 .2 * pnorm(log(1.5), 1.5, 1, lower.tail=FALSE))
    expect_equal(effective_detection_rate(.2, 1.5, 1, 3),
                 .2 * pnorm(log(3.5), 1.5, 1, lower.tail=FALSE))
    expect_lt(effective_detection_rate(.2, 1.5, 1, 3),
              effective_detection_rate(.2, 1.5, 1, 1))
  })
}

test_that("probability reconstruction has a truth-independent interface", {
  expect_true(exists("lesson_probability_draws", mode="function"))
})
if(exists("lesson_probability_draws", mode="function")) {
  test_that("occupancy is reconstructed on each draw before averaging", {
    f <- list(X_psi=matrix(c(-1,1),2,1),
              results_output=list(jsdm_output=list(
                B0_output=array(c(-1,2),c(1,2,1)),
                B_output=array(c(.5,1.5),c(1,1,2,1)),
                U_output=array(0,c(2,1,2,1)),
                L_output=array(1,c(1,1,2,1)))))
    d <- lesson_probability_draws(f)
    expect_equal(d[,1,1],plogis(c(-1.5,-.5)))
    expect_equal(d[,2,1],plogis(c(.5,3.5)))
    expect_equal(dim(d),c(2L,2L,1L))
  })
}

test_that("fit identity verification checks the data actually fitted", {
  expect_true(exists("validate_lesson_fit_identity",mode="function"))
})
if(exists("validate_lesson_fit_identity",mode="function")) {
  test_that("a changed binary outcome is rejected even with correct labels", {
    z <- matrix(c(1,0,0,1),2,2,dimnames=list(c("a","b"),c("sp1","sp2")))
    b <- list(sim=list(true_params=list(z_true=z)))
    f <- list(infos=list(model="binary",siteNames=c("b","a"),
                         speciesNames=c("sp2","sp1"),OTU=z[c("b","a"),c("sp2","sp1")]))
    expect_silent(validate_lesson_fit_identity(f,b))
    f$infos$OTU[1,1] <- 1-f$infos$OTU[1,1]
    expect_error(validate_lesson_fit_identity(f,b),"fitted observations")
  })
  test_that("fitter-added identity columns do not change the original inputs", {
    info <- data.frame(Site=c(2,1),Sample=c(2,1),Primer=c(1,1),X_theta=c(.8,.2))
    reads <- matrix(c(7,0),2,1,dimnames=list(NULL,"sp1"))
    b <- list(sim=list(data_list=list(info=info,OTU=reads)))
    fitted_info <- info[2:1,]; fitted_info$SiteSample <- c("1-1","2-2")
    f <- list(infos=list(model="two_stage",speciesNames="sp1",data_info=fitted_info,
                         OTU=reads[2:1,,drop=FALSE]))
    expect_silent(validate_lesson_fit_identity(f,b))
    f$infos$data_info$X_theta[1] <- 99
    expect_error(validate_lesson_fit_identity(f,b),"identities")
  })
}

bundle_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
if(file.exists(bundle_path)) {
  test_that("the shipped bundle retains software and summary provenance", {
    b <- readRDS(bundle_path)
    expect_true(all(vapply(b$manifests,function(m)!is.null(m$session$R.version),logical(1))))
    expect_true(length(b$summary_source_hashes)>=2)
  })
}
