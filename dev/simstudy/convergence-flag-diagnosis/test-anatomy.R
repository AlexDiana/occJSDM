# Standalone research tests for anatomy.R (Task 1 of PLAN.md in this
# directory). Run with `Rscript test-anatomy.R` from this directory. Tests on
# the saved pr11 fit of design-qfar_K6-sites300-05 are skipped when the
# archives are absent; they read the fit and its input, never write them.
library(testthat)
source('anatomy.R')

archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
have_archives <- all(dir.exists(file.path(archives,c('pr11-current-20260927','intercept-prior-inputs'))))
need_archives <- function() if(!have_archives) skip('saved study archives not available')

# Synthetic anatomy rows for one quantity of one species, from an
# iterations x chains matrix of draws.
synthetic <- function(draws,quantity='B0',species=1L,key='synthetic')
  anatomy_rows(draws,key=key,species=species,quantity=quantity,truth=NA_real_)
normal_chains <- function(means,sd=1,n=4000L) sapply(means,function(m) rnorm(n,m,sd))

test_that('(a) two chains centred at 0 and two at 5 are separated', {
  set.seed(101)
  a <- rbind(synthetic(normal_chains(c(0,0,5,5))),
    synthetic(normal_chains(c(.2,.2,.7,.7),sd=.01),'mean_psi_original_sites'))
  s <- chain_separation(a)$separation
  row <- s[s$quantity=='B0',]
  expect_identical(row$label,'separated')
  expect_gt(row$separation,3)
  expect_lte(row$max_split_rhat_within_chain,1.05)
  expect_gt(row$rhat,1.05)
})

test_that('(b) four chains drifting linearly from 0 to 5 are drifting', {
  set.seed(102)
  drift <- sapply(1:4,function(i) seq(0,5,length.out=4000L)+rnorm(4000L))
  s <- chain_separation(synthetic(drift))$separation
  expect_identical(s$label,'drifting')
  expect_gt(s$max_split_rhat_within_chain,1.05)
})

test_that('(c) four identical-distribution chains agree', {
  set.seed(103)
  s <- chain_separation(synthetic(normal_chains(c(0,0,0,0))))$separation
  expect_identical(s$label,'agrees')
  expect_lt(s$separation,3)
  expect_lt(s$rhat,1.01)
  expect_gt(s$ess_bulk,10000)
})

test_that('per-chain rows hold the chain summaries and the all-chain diagnostics', {
  set.seed(104)
  x <- normal_chains(c(0,1,2),sd=c(1,2,3))
  a <- synthetic(x,'theta0',species=3L,key='k')
  expect_identical(a$chain,1:3)
  expect_equal(a$chain_mean,unname(colMeans(x)))
  expect_equal(a$chain_sd,unname(apply(x,2,sd)))
  expect_equal(a$split_rhat_within_chain,
    unname(apply(x,2,function(v) posterior::rhat(matrix(v,ncol=1L)))))
  expect_equal(unique(a$rhat_all_chains),posterior::rhat(x))
  expect_equal(unique(a$ess_bulk_all_chains),posterior::ess_bulk(x))
  expect_true(all(c('key','species','chain','quantity','chain_mean','chain_sd',
    'split_rhat_within_chain','truth') %in% names(a)))
})

test_that('separation is the chain-mean range over the pooled within-chain SD', {
  set.seed(105)
  x <- normal_chains(c(0,1,2,6),sd=c(1,2,1,2))
  s <- chain_separation(rbind(synthetic(x),
    synthetic(normal_chains(c(.2,.2,.7,.7),sd=.01),'mean_psi_original_sites')))$separation
  s <- s[s$quantity=='B0',]
  expect_equal(s$gap,diff(range(colMeans(x))))
  expect_equal(s$pooled_within_sd,sqrt(mean(apply(x,2,sd)^2)))
  expect_equal(s$separation,s$gap/s$pooled_within_sd)
})

test_that('chain groups split mean_psi_original_sites chain means at the largest gap (R2)', {
  set.seed(106)
  a <- rbind(synthetic(normal_chains(c(0,5,.1,5.1)),'theta0',species=6L),
    synthetic(normal_chains(c(.2,.7,.21,.69),sd=.01),'mean_psi_original_sites',species=6L),
    synthetic(normal_chains(c(0,0,0,0)),'theta0',species=2L),
    synthetic(normal_chains(c(.4,.4,.4,.4),sd=.01),'mean_psi_original_sites',species=2L))
  g <- chain_separation(a)$chain_groups
  expect_identical(unique(g$species),6L)
  g <- g[order(g$chain),]
  expect_identical(g$chain,1:4)
  expect_identical(g$chain_group,c(1L,2L,1L,2L))
  expect_identical(unique(g$species_label),'separated')
})

test_that('loadings alone do not make a species separated (R6)', {
  set.seed(107)
  a <- rbind(synthetic(normal_chains(c(-1,-1,1,1),sd=.1),'L1'),
    synthetic(normal_chains(c(0,0,0,0)),'B0'),
    synthetic(normal_chains(c(.3,.3,.3,.3),sd=.01),'mean_psi_original_sites'))
  s <- chain_separation(a)
  expect_identical(s$separation$label[s$separation$quantity=='L1'],'separated')
  expect_identical(s$species$label,'agrees')
  expect_identical(s$species$separated_loadings,'L1')
  expect_identical(s$species$separated_quantities,'')
})

test_that('constant draws (loadings fixed by the constraint) agree and are marked fixed', {
  a <- synthetic(matrix(1,100L,4L),'L1')
  s <- chain_separation(a)$separation
  expect_identical(s$label,'agrees')
  expect_true(s$fixed)
  expect_identical(s$separation,0)
})

test_that('the selected fit and remapped input are found for a key', {
  need_archives()
  f <- selected_fit('design-qfar_K6-sites300-05')
  expect_identical(basename(f),'design-qfar_K6-sites300-05-fit.rds')
  expect_identical(basename(dirname(f)),'long')
  i <- selected_input('design-qfar_K6-sites300-05')
  expect_false(startsWith(i,'/Users/douglasyu/Documents'))
  expect_identical(unname(tools::md5sum(i)),'574ed4df46b9c1bd227c79d2185a7fcb')
})

test_that('an input whose md5 differs from the fit job record is refused', {
  need_archives()
  expect_error(chain_anatomy(selected_fit('design-qfar_K6-sites300-05'),
    selected_input('design-qfar_K6-sites300-01')),'md5')
})

test_that('(d) reconstructed occupancy draws average to the stored psi_output', {
  need_archives()
  fit <- selected_fit('design-qfar_K6-sites300-05');input_path <- selected_input('design-qfar_K6-sites300-05')
  a <- chain_anatomy(fit,input_path)
  stored <- readRDS(fit)$fit$results_output$psi_output
  difference <- max(abs(attr(a,'psi_mean')-stored))
  cat(sprintf('\n(d) max |reconstructed - stored psi_output| = %.3e\n',difference))
  expect_lt(difference,1e-10)
  expect_setequal(unique(a$quantity),ANATOMY_QUANTITIES)
  expect_identical(nrow(a),length(ANATOMY_QUANTITIES)*10L*4L)
  expect_identical(sort(unique(a$chain)),1:4)
  expect_identical(unique(a$key),'design-qfar_K6-sites300-05')
  # Truths taken independently of anatomy.R from the input.
  input <- readRDS(input_path);tp <- input$sim$true_params
  first <- function(q) {x <- a[a$quantity==q & a$chain==1L,];x$truth[order(x$species)]}
  expect_equal(first('mean_psi_original_sites'),unname(colMeans(plogis(tp$jsdmParams_true$eta[1:100,]))),tolerance=1e-12)
  expect_equal(first('theta0'),input$truth$params$theta0)
  expect_equal(first('B0'),tp$jsdmParams_true$B0)
  expect_equal(first('B_slope2'),tp$jsdmParams_true$B[2,])
  expect_equal(first('L2'),tp$jsdmParams_true$L[2,])
})

test_that('several fit files of one key concatenate their chains in order (R4)', {
  need_archives()
  fit <- selected_fit('design-qfar_K6-sites300-05')
  a <- chain_anatomy(c(fit,fit),selected_input('design-qfar_K6-sites300-05'))
  expect_identical(sort(unique(a$chain)),1:8)
  first <- a[a$chain<=4L,];second <- a[a$chain>4L,]
  first <- first[order(first$quantity,first$species,first$chain),]
  second <- second[order(second$quantity,second$species,second$chain),]
  expect_identical(second$chain-4L,first$chain)
  expect_identical(second$file_chain,first$file_chain)
  expect_equal(second$chain_mean,first$chain_mean)
  expect_equal(second$split_rhat_within_chain,first$split_rhat_within_chain)
  stored <- readRDS(fit)$fit$results_output$psi_output
  expect_lt(max(abs(attr(a,'psi_mean')-stored)),1e-10)
  expect_identical(nrow(attr(a,'psi_check')),2L)
})
