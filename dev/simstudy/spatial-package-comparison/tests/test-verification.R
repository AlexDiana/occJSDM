library(testthat)
source("../adapter-contract.R")
source("../verify.R")

test_that("independent draw arithmetic checks intervals, contrasts and Monte Carlo errors", {
  set.seed(914)
  ids <- list(site_id=paste0("s",1:30),species_id=c("a","b"))
  field <- array(rnorm(30*2*120*4),c(30,2,120,4))
  draws <- list(p_full=plogis(field),p_space_zero=plogis(field*.3),
    p_env_zero=plogis(field*.8),p_both_zero=field*0+.5,field=field)
  ex <- function(idx) lapply(draws,function(x)x[idx,,,,drop=FALSE])
  p <- summarize_draw_blocks(ids,ex,"logit","fixture")
  keep <- c(2L,29L)
  audit <- list(site_id=ids$site_id[keep],species_id=ids$species_id,
    draws=lapply(draws,function(x)x[keep,,,,drop=FALSE]))
  expect_lte(verify_draw_summaries(audit,p),1e-12)
  for(group in c("intervals","mcse")) {
    bad <- p; bad[[group]][[1]][2,1] <- bad[[group]][[1]][2,1]+.01
    expect_error(verify_draw_summaries(audit,bad),"summary mismatch")
  }
  for(field in c("delta_p","field_mean","field_median")) {
    bad <- p; bad[[field]][2,1] <- bad[[field]][2,1]+.01
    expect_error(verify_draw_summaries(audit,bad),"summary mismatch")
  }
  bad <- audit; bad$site_id[1] <- "absent"
  expect_error(verify_draw_summaries(bad,p),"IDs")
  bad <- p; bad$mcse$field[2,1] <- NA_real_
  expect_error(verify_draw_summaries(audit,bad),"summary mismatch")
  # This comparison uses all retained draws at once, across a 25-site block boundary.
  expect_equal(p$field_median,apply(draws$field,c(1,2),median))
  expect_equal(p$intervals$probability_lower,apply(draws$p_full,c(1,2),quantile,.025))
})
