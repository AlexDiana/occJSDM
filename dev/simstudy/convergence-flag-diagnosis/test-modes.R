# Standalone research tests for modes.R (Task 3 of PLAN.md in this directory).
# Run with `Rscript test-modes.R` from this directory. Synthetic draws only:
# no fit or archive is read.
library(testthat)
source('modes.R')

# Species-6-like draws on the scales of the pr11 fit of community 5 (Task 1:
# chains {1,3} theta0 about 0.035, B0 -0.49, mean_psi 0.475; chains {2,4}
# theta0 0.25, B0 -0.02, mean_psi 0.45). `labels` is an iterations x chains
# matrix of the generating cluster (1 = low theta0, 2 = high theta0).
cluster_draws <- function(labels,seed) {
  set.seed(seed)
  n <- length(labels);low <- as.vector(labels)==1L
  shape <- function(x) matrix(x,nrow(labels),ncol(labels))
  list(theta0=shape(ifelse(low,plogis(rnorm(n,qlogis(.035),.25)),plogis(rnorm(n,qlogis(.25),.15)))),
    B0=shape(ifelse(low,rnorm(n,-.49,.35),rnorm(n,-.02,.35))),
    mean_psi_original_sites=shape(ifelse(low,plogis(rnorm(n,qlogis(.475),.12)),plogis(rnorm(n,qlogis(.45),.12)))))
}

test_that('a known 70/30 split mixed within chains is recovered within 2 percentage points', {
  set.seed(301)
  labels <- matrix(sample(rep(1:2,c(11200L,4800L))),2000L,8L)
  r <- assign_modes(cluster_draws(labels,302),MODE_QUANTITIES)
  expect_identical(r$n_modes,2L)
  expect_identical(dim(r$labels),dim(labels))
  overall <- r$overall
  expect_lt(abs(overall$share[overall$mode==1L]-.7),.02)
  expect_lt(abs(overall$share[overall$mode==2L]-.3),.02)
  expect_gt(mean(r$labels==labels),.98)
  # Per-chain shares follow each chain's own realised split.
  truth <- colMeans(labels==1L);got <- r$share$share[r$share$mode==1L]
  expect_lt(max(abs(got-truth)),.02)
})

test_that('a 70/30 split between chains is recovered, each chain in one mode', {
  labels <- matrix(rep(c(1L,2L,1L,2L,1L,1L,2L,1L,1L,1L),each=1500L),1500L,10L)
  r <- assign_modes(cluster_draws(labels,311),MODE_QUANTITIES)
  expect_identical(r$n_modes,2L)
  expect_lt(abs(r$overall$share[r$overall$mode==1L]-.7),.02)
  regions <- chain_regions(r)
  expect_identical(regions$region,ifelse(labels[1,]==1L,'1','2'))
  expect_false(any(regions$visits_both))
  # Mode 1 is the low-theta0 mode.
  expect_lt(r$components$mean_theta0[1],r$components$mean_theta0[2])
})

test_that('unimodal Gaussian draws give a single mode', {
  set.seed(321)
  n <- 2000L;K <- 8L
  d <- list(theta0=matrix(plogis(rnorm(n*K,qlogis(.1),.3)),n,K),B0=matrix(rnorm(n*K,-.3,.4),n,K),
    mean_psi_original_sites=matrix(plogis(rnorm(n*K,0,.1)),n,K))
  r <- assign_modes(d,MODE_QUANTITIES)
  expect_identical(r$n_modes,1L)
  expect_true(all(r$labels==1L))
  expect_equal(r$share$share,rep(1,K))
  expect_true(all(chain_regions(r)$region=='1'))
})

test_that('skewed unimodal and correlated draws give a single mode', {
  set.seed(331)
  n <- 3000L;K <- 4L
  z <- rnorm(n*K)
  d <- list(theta0=matrix(rbeta(n*K,2,40),n,K),B0=matrix(z+rnorm(n*K,0,.3),n,K),
    mean_psi_original_sites=matrix(plogis(.5*z+rnorm(n*K,0,.2)),n,K))
  expect_identical(assign_modes(d,MODE_QUANTITIES)$n_modes,1L)
})

test_that('two overlapping components (a unimodal mixture) give a single mode', {
  set.seed(341)
  n <- 4000L;K <- 4L;lab <- sample(1:2,n*K,TRUE)
  shift <- ifelse(lab==1L,-.5,.5)
  d <- list(theta0=matrix(plogis(qlogis(.1)+.3*(shift+rnorm(n*K))),n,K),
    B0=matrix(rnorm(n*K),n,K),mean_psi_original_sites=matrix(plogis(rnorm(n*K,0,.1)),n,K))
  expect_identical(assign_modes(d,MODE_QUANTITIES)$n_modes,1L)
})

test_that('a constant quantity is dropped and reported, not clustered', {
  labels <- matrix(rep(1:2,c(6000L,2000L)),2000L,4L)
  d <- cluster_draws(labels,351)
  d$theta0[] <- .038
  r <- assign_modes(d,MODE_QUANTITIES)
  expect_identical(r$dropped_constant,'theta0')
  expect_identical(r$quantities_used,c('B0','mean_psi_original_sites'))
  d$B0[] <- 0;d$mean_psi_original_sites[] <- .5
  r <- assign_modes(d,MODE_QUANTITIES)
  expect_identical(r$n_modes,1L)
  expect_identical(r$quantities_used,character())
})

test_that('share tables have one row per chain and mode and sum to one per chain', {
  labels <- matrix(sample(rep(1:2,c(2800L,1200L))),1000L,4L)
  r <- assign_modes(cluster_draws(labels,361),MODE_QUANTITIES)
  expect_identical(names(r$share),c('chain','mode','draws','share'))
  expect_identical(nrow(r$share),8L)
  expect_equal(as.vector(tapply(r$share$share,r$share$chain,sum)),rep(1,4))
  expect_identical(sum(r$share$draws),4000L)
  m <- mode_mass_table(r,run='synthetic')
  expect_true(all(c('run','chain','mode','draws','share','n_modes','visits_both','region') %in% names(m)))
})

test_that('malformed input is refused', {
  d <- cluster_draws(matrix(1L,100L,2L),371)
  expect_error(assign_modes(d,c('theta0','nope')),'not in draws')
  d2 <- d;d2$B0 <- d2$B0[1:50,]
  expect_error(assign_modes(d2,MODE_QUANTITIES),'same iterations x chains')
  d3 <- d;d3$theta0[1] <- 1.5
  expect_error(assign_modes(d3,MODE_QUANTITIES),'outside \\[0, 1\\]')
})

test_that('the chain strip figure runs for 16 chains, with and without modes', {
  labels <- matrix(rep(rep(1:2,8),each=500L),500L,16L)
  d <- cluster_draws(labels,381)
  r <- assign_modes(d,MODE_QUANTITIES)
  f1 <- tempfile(fileext='.png');f2 <- tempfile(fileext='.png')
  expect_silent(plot_chain_strips(d,f1,truth=c(theta0=.038,B0=-.369,mean_psi_original_sites=.51),modes=r,
    title='synthetic 16-chain check'))
  expect_silent(plot_chain_strips(d,f2))
  expect_gt(file.size(f1),5000);expect_gt(file.size(f2),5000)
  unlink(c(f1,f2))
})

test_that('a split in theta0 alone is found when the other quantities are correlated noise', {
  # The first principal component of the standardised draws follows the
  # correlated B0 and occupancy noise, so a start from it alone misses the
  # theta0 split; the per-quantity starts find it.
  set.seed(391)
  n <- 3000L;K <- 4L;labels <- matrix(rep(c(1L,2L,1L,2L),each=n),n,K)
  low <- as.vector(labels)==1L;shared <- rnorm(n*K)
  d <- list(theta0=matrix(ifelse(low,rbeta(n*K,3,80),rbeta(n*K,25,75)),n,K),
    B0=matrix(-.3+.6*shared+rnorm(n*K,0,.15),n,K),
    mean_psi_original_sites=matrix(plogis(.5*shared+rnorm(n*K,0,.1)),n,K))
  r <- assign_modes(d,MODE_QUANTITIES,c(theta0='identity',B0='identity',mean_psi_original_sites='identity'))
  expect_identical(r$n_modes,2L)
  expect_gt(mean(r$labels==labels),.99)
})

test_that('draws of a long-tailed low mode and a compact high mode are not mislabelled in the tails', {
  # Shaped like species 6 in the pr11 fit (Task 1 data): the low mode mixes
  # near-zero theta0 draws with a main body (quantiles 0.1%, 50%, 99% about
  # 2e-5, 0.027, 0.10 against 1e-4, 0.031, 0.106 in the fit), the high mode
  # is centred at 0.25 (1% quantile 0.16), theta0 and B0 are correlated with
  # opposite signs in the two modes, and B0 and occupancy are correlated 0.96.
  # On the logit scale with a single start, 0.3% to 0.7% of a high-mode
  # chain's draws were labelled mode 1 (1.4% in the pr11 fit itself).
  set.seed(401)
  n <- 6000L;K <- 4L;labels <- matrix(rep(c(1L,2L,1L,2L),each=n),n,K)
  low <- as.vector(labels)==1L;N <- n*K
  t_low <- ifelse(runif(N)<.25,rbeta(N,1,200),rbeta(N,3,75));t_high <- rbeta(N,28,84)
  theta0 <- ifelse(low,t_low,t_high)
  zt <- ifelse(low,(t_low-mean(t_low))/sd(t_low),(t_high-mean(t_high))/sd(t_high));rho <- ifelse(low,-.55,.45)
  zB <- rho*zt+sqrt(1-rho^2)*rnorm(N);B0 <- ifelse(low,-.5,-.01)+.65*zB
  psi <- plogis(ifelse(low,qlogis(.47),qlogis(.45))+.4*(.96*zB+.28*rnorm(N)))
  r <- assign_modes(list(theta0=matrix(theta0,n,K),B0=matrix(B0,n,K),mean_psi_original_sites=matrix(psi,n,K)),MODE_QUANTITIES)
  expect_identical(r$n_modes,2L)
  regions <- chain_regions(r)
  expect_identical(regions$region,c('1','2','1','2'))
  expect_lt(max(pmin(regions$share_mode1,regions$share_mode2)),.002)
  expect_gt(mean(r$labels==labels),.999)
})
