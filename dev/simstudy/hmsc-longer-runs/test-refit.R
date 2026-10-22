# Standalone tests for the Hmsc longer-runs check. Run with
# `Rscript test-refit.R` from this directory. Tests that need the extension
# archive, the job selection or the study library are skipped when absent;
# they read those files and write only to temporary directories.
library(testthat)
source('common.R')
source('launch.R')
source('compare.R')

archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
root <- file.path(archives,'lesson-4-extension-20260923')
study <- file.path(archives,'hmsc-longer-runs-20261013')
have_root <- file.exists(file.path(root,'status.csv')) && dir.exists(file.path(root,'jobs'))
have_jobs <- file.exists(file.path(SCRIPT_REPO,JOBS_FILE))
have_hmsc <- dir.exists(file.path(study,'library/Hmsc'))
need <- function(ok,why) if(!ok) skip(why)

test_that('the schedule is four times the longer attempt with the same retained draws, and seeds are new', {
  expect_identical(iterations_per_chain(SCHEDULE),4L*iterations_per_chain(ORIGINAL_LONGER))
  expect_identical(SCHEDULE$draws,ORIGINAL_LONGER$draws)
  expect_identical(SCHEDULE$chains,4L)
  expect_identical(iterations_per_chain(SCHEDULE),96000L)
  expect_identical(job_seed(1L),26100103L)
  expect_false(any(job_seed(1:640) %in% c(26100000L+100L*(1:640)+1L,26100000L+100L*(1:640)+2L)))
})

test_that('the job selection is the frozen one and matches the archive', {
  need(have_root && have_jobs,'extension archive or job selection not available')
  jobs <- read_jobs(SCRIPT_REPO)
  expect_identical(nrow(jobs),46L)
  expect_identical(anyDuplicated(jobs$job_number),0L)
  expect_true(all(jobs$coef_fail>0))
  expect_identical(sum(jobs$badly_unmixed),14L)
  expect_equal(max(jobs$max_coef_rhat),1.293,tolerance=.001)
  status <- utils::read.csv(file.path(root,'status.csv'),stringsAsFactors=FALSE)
  m <- status[match(jobs$job,status$job),]
  expect_true(all(m$package=='Hmsc' & m$fit_ok & !m$diagnostic_pass))
  expect_identical(m$input_md5,jobs$input_md5)
  expect_identical(m$job_number,jobs$job_number)
  for(j in jobs$job[c(1,46)]) {
    t <- utils::read.csv(file.path(root,'jobs',j,'attempt-2-diagnostics.csv'),stringsAsFactors=FALSE)
    expect_identical(sum(!t$pass & t$block=='coefficient'),jobs$coef_fail[jobs$job==j])
  }
  q <- parse_jobs('all',jobs);expect_identical(q$job_number,sort(jobs$job_number));expect_identical(q$label[1],sprintf('job-%03d',q$job_number[1]))
  two <- parse_jobs(paste(jobs$job_number[2:1],collapse=','),jobs);expect_identical(two$job_number,sort(jobs$job_number[1:2]))
  expect_error(parse_jobs('999999',jobs),'not in the selection')
})

test_that('the extension code and the frozen worker are the recorded ones', {
  need(have_root,'extension archive not available')
  expect_silent(check_fit_job(root))
  code <- load_code(SCRIPT_REPO)
  expect_true(all(c('extract_bayesian','parameter_diagnostics','probability_intervals','marginal_truth','fixed_grid','basis','checked_posterior_mean') %in% ls(code)))
})

test_that('a tiny Hmsc fit through the frozen statements is reproducible, extracts and compares', {
  need(have_root && have_jobs && have_hmsc,'archive, job selection or study library not available')
  load_hmsc(study);code <- load_code(SCRIPT_REPO)
  jobs <- read_jobs(SCRIPT_REPO);j <- jobs[jobs$n_sites==100 & !jobs$use_traits,][1,]
  input <- job_input(root,j);S <- as.integer(j$n_species)
  tiny <- list(warmup=10L,draws=12L,thin=1L,chains=4L)
  fit <- fit_hmsc(input,tiny,job_seed(j$job_number))
  again <- fit_hmsc(input,tiny,job_seed(j$job_number))
  expect_identical(fit$postList[[1]][[12]]$Beta,again$postList[[1]][[12]]$Beta)
  other <- fit_hmsc(input,tiny,job_seed(j$job_number)+1L)
  expect_false(identical(fit$postList[[1]][[12]]$Beta,other$postList[[1]][[12]]$Beta))
  p <- code$extract_bayesian(fit,'Hmsc')
  expect_identical(dim(p$beta),c(3L,S,12L,4L));expect_identical(p$link,'probit')
  d <- code$parameter_diagnostics(p,input)
  expect_true(all(c('coefficient','covariance','new_site_probability') %in% d$table$block))
  expect_identical(sum(d$table$block=='coefficient'),3L*S)
  grid <- code$fixed_grid();raw <- code$basis(grid,j$response)
  x <- as.data.frame(sweep(sweep(as.matrix(raw),2,input$centre,'-'),2,input$spread,'/'))
  pi <- code$probability_intervals(p,x);expect_identical(nrow(pi),5L*S);expect_true(all(pi$lower<=pi$estimate & pi$estimate<=pi$upper))
  cs <- coefficient_summary(p$beta);expect_identical(dim(cs$mean),c(3L,S));expect_true(all(cs$lower<=cs$upper))
  est <- code$checked_posterior_mean(p,input$test_x,input$x[c(1,21,41,61,81),,drop=FALSE],
    matrix(d$table$mean[d$table$block=='new_site_probability'],nrow=5))$mean
  expect_identical(dim(est),c(300L,S))
})
