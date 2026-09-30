# Standalone research tests for data-features.R (Task 2 of PLAN.md in this
# directory). Run with `Rscript test-data-features.R` from this directory.
# Tests on the saved pr11 fit and input of design-qfar_K6-sites300-05 are
# skipped when the archives are absent; they read them, never write them.
library(testthat)
source('data-features.R')

archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
have_archives <- all(dir.exists(file.path(archives,c('pr11-current-20260927','intercept-prior-inputs'))))
need_archives <- function() if(!have_archives) skip('saved study archives not available')

KEY <- 'design-qfar_K6-sites300-05'
real_input <- local({cache <- NULL;function() {
  if(is.null(cache)) cache <<- readRDS(selected_input(KEY))
  cache
}})

# A tiny input built by hand. Four sites, two field samples per site, two
# primers, two PCRs per sample and primer (so 32 rows), two species. Rows run
# in Site, Sample, Primer, PCR order; reads are listed per sample as
# primer 1 (PCR 1, PCR 2) then primer 2 (PCR 1, PCR 2).
tiny_input <- function() {
  info <- data.frame(Site=rep(1:4,each=8),Sample=rep(1:8,each=4),Primer=rep(rep(1:2,each=2),8))
  reads1 <- c(
    5,3,0,2,   # sample 1, site 1 (occupied), collected: 3 positives
    0,1,0,0,   # sample 2, site 1, not collected: 1 positive (reads exactly 1)
    0,0,0,0,   # sample 3, site 2 (occupied), not collected
    0,0,0,0,   # sample 4, site 2, not collected
    7,0,9,9,   # sample 5, site 3 (unoccupied), collected: 3 positives
    0,0,0,0,   # sample 6, site 3, not collected
    2,0,0,0,   # sample 7, site 4 (unoccupied), not collected: 1 positive
    0,0,0,3)   # sample 8, site 4, not collected: 1 positive
  OTU <- cbind(sp1=reads1,sp2=rep(0,32))
  z <- cbind(c(1L,1L,0L,0L),c(0L,1L,1L,0L))
  w <- cbind(c(1,0,0,0,1,0,0,0),rep(0,8))
  eta <- cbind(c(2,1,-1,-2),c(0,0,0,0))
  list(sim=list(data_list=list(info=info,OTU=OTU),
      true_params=list(z_true=z,w_true=w,jsdmParams_true=list(eta=eta))),
    design=list(original_sites=2L))
}

test_that('per-site positive counts equal a hand count of the observation array', {
  need_archives()
  # Read off the raw reads of OTU_6 for sites 1, 2 and 7 (printed with the
  # site, sample and primer of each row) and counted by eye: a PCR is positive
  # if it has at least one read. Site 1: primer 1 has reads 1 and 18 in
  # sample 1 and 5 in sample 2, primer 2 has 4 in sample 2. Site 2: primer 1
  # has 221, 52, 20 in sample 3 and 9 in sample 4, primer 2 has 73, 96, 231 in
  # sample 3 and nothing in sample 4. Site 7: primer 1 has 19 and 17 in sample
  # 13 and 1 in sample 14, primer 2 has 10 in sample 13.
  f <- site_features(real_input(),6L)
  expect_identical(nrow(f),300L)
  s1 <- f[f$site==1L,];s2 <- f[f$site==2L,];s7 <- f[f$site==7L,]
  expect_identical(c(s1$positive_pcr_primer1,s1$positive_pcr_primer2,s1$positive_pcr_total,s1$positive_samples),c(3L,1L,4L,2L))
  expect_identical(c(s2$positive_pcr_primer1,s2$positive_pcr_primer2,s2$positive_pcr_total,s2$positive_samples),c(4L,3L,7L,2L))
  expect_identical(c(s7$positive_pcr_primer1,s7$positive_pcr_primer2,s7$positive_pcr_total,s7$positive_samples),c(3L,1L,4L,2L))
  # The per-sample split of site 2 (samples 3 and 4).
  expect_identical(c(s2$positives_sample1_primer1,s2$positives_sample1_primer2,
    s2$positives_sample2_primer1,s2$positives_sample2_primer2),c(3L,3L,1L,0L))
  # Generating truth of the same sites: unoccupied, one collected sample; and
  # occupied with no collected sample.
  expect_identical(c(s2$true_occupied,s2$collected_sample1,s2$collected_sample2),c(0L,1L,0L))
  expect_identical(c(s7$true_occupied,s7$collected_sample1,s7$collected_sample2),c(1L,0L,0L))
})

test_that('a read count of exactly the threshold is positive, and the threshold is honoured', {
  need_archives()
  f1 <- site_features(real_input(),6L)
  f2 <- site_features(real_input(),6L,threshold=2)
  # Site 1 has a single read in one primer 1 PCR of sample 1: positive at
  # threshold 1 (the fitted threshold), not at threshold 2.
  expect_identical(f1$positive_pcr_primer1[f1$site==1L],3L)
  expect_identical(f2$positive_pcr_primer1[f2$site==1L],2L)
  expect_identical(READ_THRESHOLD,1)
})

test_that('positive counts of every site agree with an independent recomputation', {
  need_archives()
  input <- real_input();info <- input$sim$data_list$info;Y <- input$sim$data_list$OTU
  for(sp in c(6L,3L,4L)) {
    f <- site_features(input,sp)
    for(s in c(1L,5L,17L,100L,101L,299L,300L)) {
      rows <- which(info$Site==s)
      hand <- vapply(1:2,function(p) {
        n <- 0L
        for(r in rows) if(info$Primer[r]==p && Y[r,sp]>=1) n <- n+1L
        n
      },integer(1))
      expect_identical(c(f$positive_pcr_primer1[s],f$positive_pcr_primer2[s]),hand)
      expect_identical(f$positive_pcr_total[s],sum(hand))
      by_sample <- tapply(Y[rows,sp]>=1,info$Sample[rows],any)
      expect_identical(f$positive_samples[s],sum(by_sample))
    }
    expect_identical(sum(f$positive_pcr_total),sum(Y[,sp]>=1))
  }
})

test_that('true occupancy, occupancy probability and collection states come from the input', {
  need_archives()
  input <- real_input();tp <- input$sim$true_params
  f <- site_features(input,6L)
  expect_identical(f$true_occupied,as.integer(tp$z_true[,6]))
  expect_equal(f$true_psi,plogis(tp$jsdmParams_true$eta[,6]))
  # Samples 2s-1 and 2s belong to site s.
  expect_identical(f$collected_sample1,as.integer(tp$w_true[seq(1,600,2),6]))
  expect_identical(f$collected_sample2,as.integer(tp$w_true[seq(2,600,2),6]))
  expect_identical(f$original_site,f$site<=100L)
  expect_identical(f$species_name,rep('OTU_6',300L))
})

test_that('the input holds the reads the pr11 fit used', {
  need_archives()
  saved <- readRDS(selected_fit(KEY))
  expect_identical(unname(saved$fit$infos$OTU),unname(real_input()$sim$data_list$OTU))
  expect_identical(saved$fit$infos$n,300L)
})

test_that('site features of the hand-built input match a hand count', {
  f <- site_features(tiny_input(),1L)
  expect_identical(f$site,1:4)
  expect_identical(f$positive_pcr_primer1,c(3L,0L,1L,1L))
  expect_identical(f$positive_pcr_primer2,c(1L,0L,2L,1L))
  expect_identical(f$positive_pcr_total,c(4L,0L,3L,2L))
  expect_identical(f$positive_samples,c(2L,0L,1L,2L))
  expect_identical(f$true_occupied,c(1L,1L,0L,0L))
  expect_identical(f$collected_sample1,c(1L,0L,1L,0L))
  expect_identical(f$collected_sample2,c(0L,0L,0L,0L))
  expect_identical(f$original_site,c(TRUE,TRUE,FALSE,FALSE))
  expect_equal(f$true_psi,plogis(c(2,1,-1,-2)))
  expect_true(all(site_features(tiny_input(),2L)$positive_pcr_total==0L))
})

test_that('positive routes split positives by true occupancy and by true collection', {
  r <- positive_routes(site_features(tiny_input(),1L))
  expect_identical(r$n_sites,4L)
  expect_identical(r$n_occupied_sites,2L)
  expect_identical(r$n_unoccupied_sites,2L)
  expect_identical(r$n_pcr,32L)
  expect_identical(r$n_positive_pcr,9L)
  expect_identical(r$positive_pcr_at_occupied_sites,4L)
  expect_identical(r$positive_pcr_at_unoccupied_sites,5L)
  expect_equal(r$share_positive_pcr_at_unoccupied_sites,5/9)
  # Unoccupied sites: sample 5 is collected (field false positive), samples 7
  # and 8 are not (laboratory false positive). Occupied sites: sample 1 is a
  # true detection, sample 2 a laboratory false positive.
  expect_identical(c(r$positive_pcr_unoccupied_collected,r$positive_pcr_unoccupied_uncollected),c(3L,2L))
  expect_identical(c(r$positive_pcr_occupied_collected,r$positive_pcr_occupied_uncollected),c(3L,1L))
  expect_identical(r$occupied_sites_without_positive,1L)
  expect_identical(r$unoccupied_sites_with_positive,2L)
  expect_identical(r$unoccupied_sites_with_collected_sample,1L)
  expect_identical(r$occupied_sites_without_collected_sample,1L)
  expect_equal(r$naive_detection_rate_sites,1/2)
  expect_equal(r$naive_false_positive_rate_sites,1)
  expect_equal(r$naive_positive_rate_sites,3/4)
  expect_equal(r$positive_pcr_rate,9/32)
  # 4 collected-sample PCR blocks of 4 PCRs: samples 1 and 5 hold 6 positives
  # in 8 PCRs; the other 6 samples hold 3 in 24.
  expect_equal(r$positive_pcr_rate_collected,6/8)
  expect_equal(r$positive_pcr_rate_uncollected,3/24)
  # Scope: only the first two (original) sites.
  o <- positive_routes(site_features(tiny_input(),1L),scope='original')
  expect_identical(c(o$n_sites,o$n_positive_pcr,o$positive_pcr_at_unoccupied_sites,o$occupied_sites_without_positive),c(2L,4L,0L,1L))
})

test_that('positive routes of species 6 agree with a direct count on the real input', {
  need_archives()
  input <- real_input();info <- input$sim$data_list$info;Y <- input$sim$data_list$OTU
  tp <- input$sim$true_params
  r <- positive_routes(site_features(input,6L))
  pos <- Y[,6]>=1;z_row <- tp$z_true[info$Site,6];w_row <- tp$w_true[info$Sample,6]
  expect_identical(r$n_positive_pcr,sum(pos))
  expect_identical(r$positive_pcr_at_unoccupied_sites,sum(pos & z_row==0))
  expect_identical(r$positive_pcr_at_occupied_sites,sum(pos & z_row==1))
  expect_identical(r$positive_pcr_unoccupied_collected,sum(pos & z_row==0 & w_row==1))
  expect_identical(r$positive_pcr_occupied_uncollected,sum(pos & z_row==1 & w_row==0))
  site_pos <- tapply(pos,info$Site,any)
  expect_identical(r$occupied_sites_without_positive,sum(tp$z_true[,6]==1 & !site_pos))
  expect_identical(r$n_occupied_sites,sum(tp$z_true[,6]))
  expect_equal(r$naive_detection_rate_sites,mean(site_pos[tp$z_true[,6]==1]))
})

test_that('comparison species are those labelled agrees with the closest prevalence', {
  labels <- data.frame(species=1:6,label=c('agrees','separated','agrees','drifting','agrees','agrees'),
    stringsAsFactors=FALSE)
  prev <- c(.51,.46,.44,.455,.30,.47)
  # Species 2 is separated, 4 drifting and 6 is the target, so the candidates
  # are 1 (0.04 from the target), 3 (0.03) and 5 (0.17), closest first.
  expect_identical(choose_comparison_species(prev,labels,target=6L,n=2L),c(3L,1L))
  prev[3] <- .40
  expect_identical(choose_comparison_species(prev,labels,target=6L,n=2L),c(1L,3L))
  expect_error(choose_comparison_species(prev,labels[1:5,],target=6L,n=2L),'target')
  expect_error(choose_comparison_species(prev,labels,target=6L,n=4L),'candidate')
})

test_that('site occupancy errors are summarised as MAE, signed error, RMSE and correlation', {
  e <- site_occupancy_errors(psi=c(.2,.6,.9,.5),truth=c(.1,.7,.5,.5),sites=1:3)
  expect_equal(e$mae,(.1+.1+.4)/3)
  expect_equal(e$signed_error,(.1-.1+.4)/3)
  expect_equal(e$rmse,sqrt((.01+.01+.16)/3))
  expect_equal(e$correlation,cor(c(.2,.6,.9),c(.1,.7,.5)))
  expect_identical(e$n_sites,3L)
})

test_that('a fit is subset to chains without changing the draws it keeps', {
  need_archives()
  fit <- readRDS(selected_fit(KEY))$fit
  sub <- subset_fit_chains(fit,c(2L,4L))
  jo <- fit$results_output$jsdm_output;js <- sub$results_output$jsdm_output
  expect_identical(dim(js$B0_output),c(10L,12000L,2L))
  expect_identical(js$B0_output[,,2],jo$B0_output[,,4])
  expect_identical(dim(js$U_output),c(300L,2L,12000L,2L))
  expect_identical(js$U_output[,,,1],jo$U_output[,,,2])
  expect_identical(sub$results_output$theta0_output[,,1],fit$results_output$theta0_output[,,2])
  expect_identical(sub$results_output$beta_theta_output[,,,2],fit$results_output$beta_theta_output[,,,4])
  expect_error(subset_fit_chains(fit,5L),'chain')
})

test_that('per-chain occupancy reconstructions average to the stored psi_output and to Task 1 chain means', {
  need_archives()
  fit <- readRDS(selected_fit(KEY))$fit;input <- real_input()
  truth <- anatomy_truth(fit,input)
  recon <- chain_reconstruction(fit,truth$original_sites)
  expect_length(recon,4L)
  total <- Reduce(`+`,lapply(recon,`[[`,'psi_sum'))/(recon[[1]]$ni*4)
  difference <- max(abs(total-fit$results_output$psi_output))
  cat(sprintf('\nmax |chain reconstructions - stored psi_output| = %.3e\n',difference))
  expect_lt(difference,1e-10)
  summary <- read.csv(file.path(ANATOMY_DIR,'results/anatomy/chain-summary.csv'),stringsAsFactors=FALSE)
  summary <- summary[summary$key==KEY & summary$species==6L,]
  for(q in c('theta0','B_slope1','mean_psi_original_sites','L2')) for(ch in 1:4) {
    expect_equal(mean(recon[[ch]]$draws[[q]][6,,1]),summary$chain_mean[summary$quantity==q & summary$chain==ch],
      tolerance=1e-5)
  }
})

test_that('group rows carry the pooled mean, the generating value and the absolute error', {
  set.seed(201)
  make <- function(m) array(rnorm(2*200,m,.1),c(2L,200L,1L))
  recon <- list(list(draws=list(B0=make(1),theta0=make(.2)),psi_sum=NULL,ni=200L),
    list(draws=list(B0=make(3),theta0=make(.6)),psi_sum=NULL,ni=200L),
    list(draws=list(B0=make(1.1),theta0=make(.22)),psi_sum=NULL,ni=200L))
  # Draws are species x iteration x chain; species 2 is used.
  truth <- c(B0=1.05,theta0=.5)
  g <- group_parameter_rows(recon,chains=c(1L,3L),species=2L,truth=truth,
    key='k',chain_group='1')
  expect_identical(g$quantity,c('B0','theta0'))
  expect_identical(unique(g$chains),'1;3')
  pooled <- mean(c(recon[[1]]$draws$B0[2,,1],recon[[3]]$draws$B0[2,,1]))
  expect_equal(g$group_mean[g$quantity=='B0'],pooled)
  expect_equal(g$abs_error[g$quantity=='B0'],abs(pooled-1.05))
  expect_equal(g$error[g$quantity=='theta0'],mean(c(recon[[1]]$draws$theta0[2,,1],
    recon[[3]]$draws$theta0[2,,1]))-.5)
  expect_true(g$truth_in_interval[g$quantity=='B0'])
  expect_false(g$truth_in_interval[g$quantity=='theta0'])
  expect_true(all(g$q025<g$group_mean & g$group_mean<g$q975))
})
