# Standalone research tests for impact.R (Task 4 of PLAN.md in this
# directory). Run with `Rscript test-impact.R` from this directory. The last
# test scores the saved pr11 fit of community 5 with the archived scorers and
# takes a few minutes; set IMPACT_FAST=1 to skip it while iterating.
library(testthat)
source('impact.R')

ctx <- impact_setup()
scores <- read_recheck('community-scores.csv')
key5 <- TARGET_KEY

# ---- Current-main-recheck definitions ------------------------------------------------

test_that('original-site groups follow the design scorer bands and error definitions', {
  # Three species, four fitted sites, the first two original. Truth sits on both
  # band edges: 0.2 and 0.8 are in the middle band (design_helpers.R:124).
  truth <- matrix(c(.1,.2,.5,.5, .8,.9,.5,.5, .3,.7,.5,.5),4,3)
  estimate <- truth+matrix(c(.05,-.1,9,9, .1,-.2,9,9, 0,.3,9,9),4,3)
  g <- recheck_groups(truth,estimate,2L,ctx$sc$error_metrics)
  expect_identical(g$group,c('all','low','medium','high'))
  expect_identical(g$metric,rep('occupancy_original_sites',4))
  expect_identical(as.integer(g$n_elements),c(6L,1L,4L,1L))
  e <- c(.05,-.1,.1,-.2,0,.3)
  expect_equal(g$bias[1],mean(e),tolerance=1e-15)
  expect_equal(g$mae[1],mean(abs(e)),tolerance=1e-15)
  expect_equal(g$rmse[1],sqrt(mean(e^2)),tolerance=1e-15)
  expect_equal(g$bias[g$group=='low'],.05,tolerance=1e-15)
  expect_equal(g$bias[g$group=='high'],-.2,tolerance=1e-15)
  expect_equal(g$mae[g$group=='medium'],mean(abs(c(-.1,.1,0,.3))),tolerance=1e-15)
  # The non-original sites (estimate error 9) never enter.
  expect_true(all(abs(g$bias)<1))
  # An empty band is left out, as add_group is only called for a non-empty band.
  g2 <- recheck_groups(matrix(.5,2,2),matrix(.6,2,2),2L,ctx$sc$error_metrics)
  expect_identical(g2$group,c('all','medium'))
})

test_that('substitution replaces only the target community rows of one version', {
  new <- subset(scores,key==key5 & version=='current' & metric=='occupancy_original_sites')
  new <- new[c('metric','group','n_elements','truth','estimate','bias','mae','rmse')]
  new$estimate <- new$estimate+.01;new$bias <- new$bias+.01;new$mae <- new$mae+.02;new$rmse <- new$rmse+.03
  s <- substitute_scores(scores,key5,'current',new)
  hit <- s$key==key5 & s$version=='current' & s$metric=='occupancy_original_sites'
  expect_identical(sum(hit),4L)
  expect_identical(s[!hit,],scores[!hit,])
  expect_equal(s$mae[hit],scores$mae[hit]+.02)
  expect_equal(s$bias[hit],scores$bias[hit]+.01)
  # Chain diagnostics of a substituted estimate are not defined.
  expect_true(all(is.na(s$rhat[hit])))
  # A substitute must keep the cell sets and truth of the rows it replaces.
  bad <- new;bad$truth[1] <- bad$truth[1]+1e-6
  expect_error(substitute_scores(scores,key5,'current',bad),'truth')
  bad <- new;bad$n_elements[2] <- bad$n_elements[2]+1L
  expect_error(substitute_scores(scores,key5,'current',bad),'cells')
  expect_error(substitute_scores(scores,key5,'current',new[1:3,]),'groups')
})

test_that('arm means, design contrasts and code-version changes reproduce the published summaries', {
  mean_ci <- ctx$recheck$mean_ci
  sens <- read_recheck('chain-sensitivity-summary.csv')
  design <- read_recheck('design-paired-summary.csv')
  paired <- read_recheck('paired-summary.csv')
  gap <- 0
  for(i in which(sens$family=='design')) {
    a <- arm_summary(scores,'current',sens$scenario[i],sens$arm[i],sens$group[i])
    gap <- max(gap,abs(a[['mae']]-sens$selected_mae[i]),abs(a[['bias']]-sens$selected_bias[i]))
  }
  expect_lt(gap,1e-10)
  gap <- 0
  for(i in which(design$family=='design')) {
    d <- design_contrast(scores,design$version[i],design$scenario[i],design$reference[i],design$alternative[i],design$group[i],mean_ci)
    expect_identical(d$improved,as.integer(design$improved[i]))
    gap <- max(gap,abs(d$reduction_mae-design$reduction_mae[i]),abs(d$lower-design$lower[i]),abs(d$upper-design$upper[i]))
  }
  expect_lt(gap,1e-10)
  gap <- 0
  rows <- which(paired$comparison=='selected' & paired$family=='design' & paired$metric=='occupancy_original_sites')
  expect_identical(length(rows),32L)
  for(i in rows) {
    v <- version_change(scores,paired$scenario[i],paired$arm[i],paired$group[i],mean_ci)
    gap <- max(gap,abs(v$delta_mae_mean-paired$delta_mae_mean[i]),abs(v$delta_mae_lower-paired$delta_mae_lower[i]),
      abs(v$delta_mae_upper-paired$delta_mae_upper[i]),abs(v$current_mae-paired$current_mae[i]))
  }
  expect_lt(gap,1e-10)
})

test_that('the chain-sensitivity ranges are read from the published tables', {
  r <- sensitivity_range('arm','qfar_K6','sites300','all','mae')
  expect_equal(unname(r),c(0.166824291212352,0.173126808754637),tolerance=1e-14)
  r <- sensitivity_range('contrast','qfar_K6','field4:sites300','all','reduction_mae')
  expect_equal(unname(r),c(-0.00115873500067254,0.00612184833368518),tolerance=1e-14)
  r <- sensitivity_range('arm','qfar_K6','sites300','low','bias')
  expect_equal(unname(r),c(0.147193381754557,0.158933968780993),tolerance=1e-14)
})

test_that('interval columns are those interval_cells takes from score_draw_block', {
  set.seed(1)
  draws <- array(runif(5*40*2),c(5,40,2));truth <- runif(5)
  a <- ctx$ip$interval_cells(draws,truth,ctx$sc)
  b <- iv_from_block(ctx$sc$score_draw_block(draws,truth,'cell')$elements)
  expect_identical(b,a)
})

# ---- Phase B of the intercept-prior study ----------------------------------------------

test_that('the phase B gate is reproduced from the archived per-community tables', {
  g <- phase_b_groups()
  conv <- phase_b_convergence()
  gate <- phase_b_gate(g,conv,ctx$ip)
  published <- read.csv(file.path(IP_RESULTS,'gate-detail-B.csv'),stringsAsFactors=FALSE)
  d <- gate$detail
  expect_identical(nrow(d),nrow(published))
  expect_identical(paste(d$sd,d$criterion,d$component,d$stratum),paste(published$sd,published$criterion,published$component,published$stratum))
  expect_identical(d$pass,published$pass)
  num <- c('communities','needed','improved','value_new','value_control','statistic','threshold','missing_new')
  expect_lt(max(abs(as.matrix(d[num])-as.matrix(published[num])),na.rm=TRUE),1e-12)
  expect_identical(is.na(as.matrix(d[num])),is.na(as.matrix(published[num])))
  rows <- read.csv(file.path(IP_RESULTS,'gate-B.csv'),stringsAsFactors=FALSE)
  expect_identical(gate$rows$result,rows$result)
  expect_identical(gate$rows$result,rep('FAIL',3))
})

test_that('phase B substitution replaces only the control rows of one community', {
  g <- phase_b_groups()
  hit <- g$sd==1 & g$key==key5
  expect_identical(sum(hit),4L)
  new <- g[hit,c('scope','group','cells','truth_mean','estimate_mean','signed_error','abs_signed_error','mean_abs_cell_error','coverage')]
  new$signed_error <- new$signed_error-1;new$abs_signed_error <- abs(new$signed_error);new$coverage <- new$coverage/2
  s <- substitute_phase_b(g,key5,new)
  expect_identical(s[!hit,],g[!hit,])
  expect_equal(s$coverage[hit],g$coverage[hit]/2)
  bad <- new;bad$cells[1] <- bad$cells[1]+1L
  expect_error(substitute_phase_b(g,key5,bad),'cells')
  # One fewer flagged control fit at qfar, nothing else changed.
  conv <- phase_b_convergence();c2 <- converged_convergence(conv,'qfar')
  at <- c2$sd==1 & c2$stratum=='qfar'
  expect_identical(c2$selected_flagged[at],conv$selected_flagged[at]-1L)
  expect_identical(c2[!at,],conv[!at,])
})

# ---- The other four flagged fits ------------------------------------------------------

test_that('flagged diagnostics are attributed to species by the scorer element layouts', {
  expect_identical(element_species('additional_diagnostics','original_probability',701L),8L)
  expect_identical(element_species('additional_diagnostics','original_probability',800L),8L)
  expect_identical(element_species('additional_diagnostics','original_probability',912L),10L)
  expect_identical(element_species('additional_diagnostics','environment_slope',16L),8L)
  expect_identical(element_species('additional_diagnostics','environment_slope',20L),10L)
  expect_identical(element_species('elements','theta0',8L),8L)
  expect_identical(element_species('elements','p',3L),2L)
  expect_identical(element_species('additional_diagnostics','sigma_h',1L),NA_integer_)
  expect_identical(element_species('groups','occupancy',NA_integer_),NA_integer_)
})

test_that('the four other flagged fits are reported with their Task 1 labels and flag sources', {
  o <- other_flagged_fits()
  expect_identical(nrow(o),40L)
  expect_setequal(unique(o$key),setdiff(FLAGGED_KEYS,key5))
  expect_false(any(o$task1_label=='separated'))
  expect_true(all(o$max_separation<1))
  expect_false(any(o$extended_run_needed))
  q2 <- o[o$key=='design-qnear_K6-sites300-02',]
  expect_identical(q2$task1_label[q2$species==8],'drifting')
  expect_identical(q2$flagged_diagnostics[q2$species==8],35L)
  expect_identical(sum(q2$flagged_diagnostics),35L)
  f4 <- o[o$key=='design-qfar_K6-field4-09',]
  expect_identical(sum(f4$flagged_diagnostics),0L)
  expect_true(all(grepl('B_output',f4$package_warnings)))
})

# ---- The saved pr11 fit of community 5 (slow) -----------------------------------------

if(!identical(Sys.getenv('IMPACT_FAST'),'1')) test_that('the pr11 fit reproduces the published community 5 scores and chain means', {
  p <- score_pr11(ctx)
  published <- subset(scores,key==key5 & version=='current')
  # Every group row of the archived design scorer, as summarise.R exported it.
  m <- merge(published,p$recheck_rows,by=c('metric','group'),suffixes=c('_published','_fit'))
  expect_identical(nrow(m),nrow(published))
  num <- c('n_elements','truth','estimate','bias','mae','rmse','rhat','ess_mean','mcse','chain_gap')
  gap <- max(vapply(num,function(k) max(abs(m[[paste0(k,'_published')]]-m[[paste0(k,'_fit')]])),numeric(1)))
  expect_lt(gap,1e-10)
  band <- recheck_groups(p$truth,p$sources$pr11_pooled$estimate,p$n0,ctx$sc$error_metrics)
  o <- published[published$metric=='occupancy_original_sites',];o <- o[match(band$group,o$group),]
  expect_lt(max(abs(c(band$bias-o$bias,band$mae-o$mae,band$rmse-o$rmse))),1e-10)
  # The phase B record of the same fit (intercept-prior control).
  rec <- readRDS(file.path(IP_ARCHIVE,'summary/B/scores/sd1',paste0(key5,'-long.rds')))
  pb <- source_scores(p$sources$pr11_pooled,p,ctx)$phase_b
  cols <- c('cells','truth_mean','estimate_mean','signed_error','abs_signed_error','mean_abs_cell_error','coverage')
  expect_identical(pb$group,rec$groups$group)
  expect_lt(max(abs(as.matrix(pb[cols])-as.matrix(rec$groups[cols]))),1e-10)
  # Chain means reproduce the published observed-chain scores (chain 0 is pooled).
  ch <- read_recheck('flagged-chain-scores.csv');ch <- ch[ch$key==key5,]
  gap <- 0
  for(k in 0:4) {
    est <- if(k==0L) p$sources$pr11_pooled$estimate else matrix(p$chain_means[,k],nrow(p$truth))
    g <- recheck_groups(p$truth,est,p$n0,ctx$sc$error_metrics)
    z <- ch[ch$chain==k,];z <- z[match(g$group,z$group),]
    gap <- max(gap,abs(g$bias-z$bias),abs(g$mae-z$mae))
  }
  expect_lt(gap,1e-10)
  # The chain-group estimates are the means of their chains' draws.
  expect_lt(max(abs(p$sources$near_truth$estimate-matrix(rowMeans(p$chain_means[,c(1,3)]),nrow(p$truth)))),1e-15)
  expect_lt(max(abs(p$sources$mirror$estimate-matrix(rowMeans(p$chain_means[,c(2,4)]),nrow(p$truth)))),1e-15)
  expect_identical(p$chain_groups,list(near_truth=c(1L,3L),mirror=c(2L,4L)))
})
