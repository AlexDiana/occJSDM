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
  f3 <- tempfile(fileext='.png')
  expect_silent(plot_chain_strips(d,f3,modes=named_result('anchored',r$labels)))
  expect_gt(file.size(f1),5000);expect_gt(file.size(f2),5000);expect_gt(file.size(f3),5000)
  unlink(c(f1,f2,f3))
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

# ---- Anchored assignment (ruling R17) -----------------------------------------------
# The primary assignment classifies each draw with a frozen two-component
# classifier calibrated on the saved pr11 fit (chains 1 and 3 near-truth,
# chains 2 and 4 mirror) from quantities that stay free in every variant.
# Pseudo-chains are blocks of that fit's species-6 draws.

anchor_quantities <- c('B_slope1','B_slope2','beta_theta_intercept','mean_psi_original_sites')
archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
have_pr11 <- all(dir.exists(file.path(archives,c('pr11-current-20260927','intercept-prior-inputs'))))
need_pr11 <- function() if(!have_pr11) skip('saved pr11 fit not available')
pr11_cache <- new.env()
pr11_draws <- function() {
  if(is.null(pr11_cache$d)) {
    anatomy <- new.env(parent=globalenv());source('anatomy.R',local=anatomy)
    key <- 'design-qfar_K6-sites300-05'
    a <- anatomy$chain_anatomy(anatomy$selected_fit(key),anatomy$selected_input(key),keep_draws=TRUE)
    pr11_cache$d <- species_mode_draws(a,6L,unique(c('theta0','B0','mean_psi_original_sites',anchor_quantities)))
  }
  pr11_cache$d
}
# nl near-truth and nh mirror pseudo-chains: each of chains {1,3} and {2,4}
# cut into nb blocks.
pseudo_chains <- function(nl,nh,nb) {
  d <- pr11_draws()
  blocks <- function(ch) lapply(d,function(m) do.call(cbind,lapply(ch,function(j) matrix(m[,j],nrow(m)/nb,nb))))
  lo <- blocks(c(1,3));hi <- blocks(c(2,4))
  out <- lapply(names(d),function(q) cbind(lo[[q]][,seq_len(nl),drop=FALSE],hi[[q]][,seq_len(nh),drop=FALSE]))
  stats::setNames(out,names(d))
}
expect_split <- function(r,nl,nh) {
  reg <- anchored_regions(r)
  expect_identical(reg$region,c(rep('near-truth',nl),rep('mirror',nh)))
  expect_lte(max(pmin(reg$share_near_truth,reg$share_mirror)),.01)
}

test_that('the frozen anchor classifier is the one calibrated on the saved pr11 fit, without theta0 or B0', {
  need_pr11()
  anchor <- read_anchor(ANCHOR_FILE)
  expect_identical(ANCHOR_QUANTITIES,anchor_quantities)
  expect_identical(anchor$quantities,anchor_quantities)
  expect_false(any(c('theta0','B0') %in% anchor$quantities))
  expect_identical(anchor$modes,c('near-truth','mirror'))
  expect_identical(anchor$weight,c(.5,.5))
  fresh <- calibrate_anchor(pr11_draws(),near_chains=c(1L,3L),mirror_chains=c(2L,4L))
  expect_identical(fresh$mean,anchor$mean);expect_identical(fresh$cov,anchor$cov)
})

test_that('imbalanced splits of 16 pseudo-chains are assigned chain by chain with at most 1% off-region', {
  need_pr11()
  anchor <- read_anchor(ANCHOR_FILE)
  for(s in list(c(15,1),c(14,2),c(2,14),c(1,15))) expect_split(assign_anchored(pseudo_chains(s[1],s[2],8),anchor),s[1],s[2])
})

test_that('a 7 and 1 split of 8 pseudo-chains is assigned correctly, either way round', {
  need_pr11()
  anchor <- read_anchor(ANCHOR_FILE)
  expect_split(assign_anchored(pseudo_chains(7,1,4),anchor),7,1)
  expect_split(assign_anchored(pseudo_chains(1,7,4),anchor),1,7)
})

test_that('theta0 held constant, as in variant (b), leaves the anchored assignment unchanged', {
  need_pr11()
  anchor <- read_anchor(ANCHOR_FILE)
  d <- pseudo_chains(15,1,8);held <- d;held$theta0[] <- .038004923556
  a <- assign_anchored(d,anchor);b <- assign_anchored(held,anchor)
  expect_identical(b$labels,a$labels)
  expect_split(b,15,1)
  expect_identical(assign_anchored(held[setdiff(names(held),c('theta0','B0'))],anchor)$labels,a$labels)
})

test_that('a classifier calibrated on chains 1 and 2 assigns held-out chains 3 and 4', {
  need_pr11()
  d <- pr11_draws()
  held_out <- calibrate_anchor(d,near_chains=1L,mirror_chains=2L)
  r <- assign_anchored(lapply(d,function(m) m[,3:4,drop=FALSE]),held_out)
  expect_split(r,1,1)
})

test_that('the theta0 cut at 0.135 assigns the pr11 pseudo-chains correctly', {
  need_pr11()
  expect_identical(THETA0_CUT,.135)
  expect_split(theta0_cut_assignment(pseudo_chains(14,2,8)),14,2)
  expect_split(theta0_cut_assignment(pseudo_chains(1,7,4)),1,7)
})

test_that('the per-chain comparison reports where the secondary assignments disagree with the primary', {
  need_pr11()
  anchor <- read_anchor(ANCHOR_FILE)
  d <- pseudo_chains(14,2,8)
  primary <- assign_anchored(d,anchor);cut <- theta0_cut_assignment(d)
  refit <- assign_modes(d[MODE_QUANTITIES])
  cmp <- compare_assignments(primary,cut,refit)
  expect_true(all(c('chain','primary_region','cut_region','refit_region','agree_cut','agree_refit') %in% names(cmp)))
  expect_true(all(cmp$agree_cut))
  # The refit with chain-partition starts recovers a 14 and 2 split.
  expect_identical(refit$n_modes,2L);expect_true(all(cmp$agree_refit))
  # A refit that put chain 1 in the other mode is reported as disagreeing there only.
  wrong <- refit;wrong$labels[,1] <- 3L-wrong$labels[,1]
  wrong$share <- chain_share(wrong$labels,2L)
  cmp2 <- compare_assignments(primary,cut,wrong)
  expect_identical(which(!cmp2$agree_refit),1L)
})

test_that('anchored assignment works on synthetic draws and flags draws far from both components', {
  set.seed(421)
  n <- 2000L;mk <- function(m,s) matrix(rnorm(n*2,m,s),n,2)
  near <- list(x=mk(0,1),y=mk(0,1));far <- list(x=mk(6,1),y=mk(-6,1))
  anchor <- calibrate_anchor(list(x=cbind(near$x,far$x),y=cbind(near$y,far$y)),near_chains=1:2,mirror_chains=3:4,
    quantities=c('x','y'))
  test <- list(x=cbind(rnorm(n),rnorm(n,6),rnorm(n,30)),y=cbind(rnorm(n),rnorm(n,-6),rnorm(n,30)))
  r <- assign_anchored(test,anchor)
  expect_identical(anchored_regions(r)$region[1:2],c('near-truth','mirror'))
  expect_lt(r$atypical$share[1],.01);expect_gt(r$atypical$share[3],.99)
  f <- tempfile(fileext='.csv');write_anchor(anchor,f)
  expect_error(read_anchor(f),'md5 mismatch')
  back <- read_anchor(f,md5=NULL);expect_identical(back$mean,anchor$mean);expect_identical(back$cov,anchor$cov)
  unlink(f)
})

test_that('region patterns follow AMENDMENT-1', {
  reg <- function(region) data.frame(chain=seq_along(region),region=region,stringsAsFactors=FALSE)
  expect_identical(region_pattern(reg(rep('near-truth',16))),'one mode only: near-truth')
  expect_identical(region_pattern(reg(rep('mirror',16))),'one mode only: mirror')
  expect_identical(region_pattern(reg(rep(c('near-truth','mirror'),8))),'each chain in one region')
  expect_identical(region_pattern(reg(rep('both',16))),'every chain visits both')
  expect_identical(region_pattern(reg(c(rep('both',15),'mirror'))),'other')
})
