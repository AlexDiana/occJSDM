#!/usr/bin/env Rscript
# Research tests for the scorer: analysis.R (phase A per-community outcomes),
# analysis-b.R (phase B), summarise.R (tables, the gates and the final
# decision) and verify.R (independent audit). Run from this directory with
# `Rscript test-analysis.R`; the exit status is 1 if any expectation fails or
# any block errors.
#
# Blind development (ruling R28): every metric and gate criterion is checked on
# hand-computed synthetic examples. Real fits enter only as mechanics fixtures:
# archived control fits (whose archived results the scorer must reproduce) and
# the phase A sd 3 pilot fits (pilot-schedule fits that are never scored; only
# their loading path is exercised and no outcome of them is printed or
# compared). The phase B new-arm loading path is exercised on a synthetic
# new-arm copy of a control fit, written to a temporary directory. No new-arm
# fit is read, and nothing under the study's fits/B is read. Tests that need
# the archives skip when absent.
if(!identical(Sys.getenv('OCCJSDM_TEST_ANALYSIS_INNER'),'1')) {
  Sys.setenv(OCCJSDM_TEST_ANALYSIS_INNER='1')
  self <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]))
  res <- testthat::test_file(self,reporter=testthat::ProgressReporter$new(show_praise=FALSE))
  d <- as.data.frame(res)
  cat(sprintf('\n%d test blocks, %d expectations, %d failed, %d blocks with errors, %d skipped\n',
    nrow(d),sum(d$nb),sum(d$failed),sum(d$error),sum(d$skipped)))
  quit(save='no',status=if(sum(d$failed)+sum(d$error)) 1L else 0L)
}

library(testthat)
for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R','analysis-b.R','supplement.R','summarise.R')) source(f)
here <- normalizePath('.')
repo <- normalizePath('../../..')
archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
inputs_root <- file.path(archives,'intercept-prior-inputs')
study <- file.path(archives,'intercept-prior-20260929')
have_archives <- all(dir.exists(file.path(archives,c(unname(ARCHIVES),'intercept-prior-inputs'))))
need_archives <- function() if(!have_archives) skip('saved study archives not available')
rscript <- file.path(R.home('bin'),'Rscript')
run_cli <- function(script,args) {
  # system2 warns only of the non-zero exit status that the refusal tests expect.
  out <- suppressWarnings(system2(rscript,c(shQuote(file.path(here,script)),shQuote(args)),stdout=TRUE,stderr=TRUE))
  list(status=attr(out,'status') %||% 0L,output=paste(out,collapse='\n'))
}
# The archived scorer definitions, md5-checked and loaded on first use.
delayedAssign('sc',if(have_archives) load_metric_scorers(repo,archives) else NULL)
delayedAssign('bsc',if(have_archives) load_b_scorers(repo,archives) else NULL)
# Type-7 sample quantile written out by hand: h = (N - 1) p + 1.
q7 <- function(x,p) {x <- sort(x);h <- (length(x)-1)*p+1;lo <- floor(h);x[lo]+(h-lo)*(x[min(lo+1,length(x))]-x[lo])}
row_of <- function(g,scope,group) {r <- g[g$scope==scope & g$group==group,,drop=FALSE];stopifnot(nrow(r)==1L);r}

# ---------------------------------------------------------------------------
# Metrics on hand-computed three-species examples
# ---------------------------------------------------------------------------

# Three species at four sites, of which sites 1 and 2 are the original sites.
# Truth (columns are species):
#   sp1 .05 .80 .10 .10   all-site mean .2625 (not rare)
#   sp2 .85 .20 .90 .95   all-site mean .725
#   sp3 .25 .15 .05 .05   all-site mean .125 (rare: below .2 over all fitted sites)
# Estimate - truth:
#   sp1 +.05 -.10 +.02 -.02
#   sp2 -.05 +.10 -.05 -.05
#   sp3 +.05 +.05 +.05 -.03
hand_a1 <- function() {
  truth <- cbind(c(.05,.80,.10,.10),c(.85,.20,.90,.95),c(.25,.15,.05,.05))
  err <- cbind(c(.05,-.10,.02,-.02),c(-.05,.10,-.05,-.05),c(.05,.05,.05,-.03))
  covered <- cbind(c(TRUE,FALSE,TRUE,TRUE),c(TRUE,FALSE,TRUE,FALSE),c(FALSE,TRUE,TRUE,TRUE))
  make_cells('A1',truth=truth,estimate=truth+err,covered=covered,n_original=2L)
}

test_that('A1 signed error, MAE and coverage by band on the original sites match the hand computation', {
  need_archives()
  g <- group_table(hand_a1(),sc)
  # low: sp1 site 1 (+.05) and sp3 site 2 (+.05)
  low <- row_of(g,'primary','low')
  expect_identical(low$cells,2L);expect_equal(low$signed_error,5,tolerance=1e-12)
  expect_equal(low$abs_signed_error,5,tolerance=1e-12);expect_equal(low$mean_abs_cell_error,5,tolerance=1e-12)
  expect_equal(low$coverage,1)
  # middle, inclusive at .2 and .8: sp1 site 2 (-.10), sp2 site 2 (+.10), sp3 site 1 (+.05)
  mid <- row_of(g,'primary','middle')
  expect_identical(mid$cells,3L);expect_equal(mid$signed_error,100*(-.10+.10+.05)/3,tolerance=1e-12)
  expect_equal(mid$mean_abs_cell_error,100*(.10+.10+.05)/3,tolerance=1e-12);expect_equal(mid$coverage,0)
  # high: sp2 site 1 (-.05); a negative signed error has a positive absolute value
  high <- row_of(g,'primary','high')
  expect_identical(high$cells,1L);expect_equal(high$signed_error,-5,tolerance=1e-12)
  expect_equal(high$abs_signed_error,5,tolerance=1e-12);expect_equal(high$coverage,1)
  # overall MAE over the six original cells
  all <- row_of(g,'primary','all')
  expect_identical(all$cells,6L);expect_equal(all$signed_error,100*.10/6,tolerance=1e-12)
  expect_equal(all$mean_abs_cell_error,100*.40/6,tolerance=1e-12);expect_equal(all$coverage,3/6)
  expect_equal(all$truth_mean,mean(c(.05,.80,.85,.20,.25,.15)),tolerance=1e-12)
})

test_that('the A1 rare group is the species below 0.2 mean truth over all fitted sites, with all its fitted sites', {
  need_archives()
  g <- group_table(hand_a1(),sc)
  rare <- row_of(g,'primary','rare_below_20pct')
  # sp3 only (sp1 has original-site mean .425 and all-site mean .2625), all four sites
  expect_identical(rare$cells,4L);expect_equal(rare$signed_error,100*(.05+.05+.05-.03)/4,tolerance=1e-12)
  expect_equal(rare$mean_abs_cell_error,100*.18/4,tolerance=1e-12);expect_equal(rare$coverage,.75)
  # no rare species: the group is absent rather than empty
  cells <- hand_a1();cells$truth[,3] <- c(.25,.15,.30,.30);cells$estimate[,3] <- cells$truth[,3]
  expect_false('rare_below_20pct' %in% group_table(cells,sc)$group)
})

test_that('A1 fits with more than the original sites are also scored on all sites, descriptively', {
  need_archives()
  g <- group_table(hand_a1(),sc)
  expect_setequal(unique(g$scope),c('primary','allsites'))
  low <- row_of(g,'allsites','low')
  expect_identical(low$cells,6L);expect_equal(low$signed_error,100*(.05+.02-.02+.05+.05-.03)/6,tolerance=1e-12)
  all <- row_of(g,'allsites','all')
  expect_identical(all$cells,12L);expect_equal(all$mean_abs_cell_error,100*(.19+.25+.18)/12,tolerance=1e-12)
  expect_equal(row_of(g,'allsites','high')$signed_error,-5,tolerance=1e-12)
  # a fit whose fitted sites are the original sites has the primary scope only
  cells <- hand_a1();cells$n_original <- 4L
  expect_identical(unique(group_table(cells,sc)$scope),'primary')
})

# Four species at four sites with designed prevalences 1%, 5%, 25% and 75%.
hand_a2 <- function() {
  truth <- cbind(c(.005,.010,.012,.013),c(.02,.03,.07,.08),c(.10,.20,.30,.40),c(.50,.70,.90,.90))
  estimate <- matrix(rep(c(.0125,.0475,.25875,.77),each=4),4,4)
  covered <- cbind(c(FALSE,TRUE,TRUE,TRUE),rep(FALSE,4),c(FALSE,TRUE,TRUE,FALSE),rep(FALSE,4))
  make_cells('A2',truth=truth,estimate=estimate,covered=covered,target_prevalence=c(.01,.05,.25,.75))
}

test_that('A2 bands, prevalence groups and the gated 1% and 5% rare group match the hand computation', {
  need_archives()
  g <- group_table(hand_a2(),sc)
  expect_identical(unique(g$scope),'primary')
  e <- hand_a2();err <- e$estimate-e$truth
  low <- row_of(g,'primary','low')
  expect_identical(low$cells,9L)
  expect_equal(low$signed_error,100*mean(c(.0075,.0025,.0005,-.0005,.0275,.0175,-.0225,-.0325,.15875)),tolerance=1e-12)
  expect_equal(low$coverage,3/9)
  expect_identical(row_of(g,'primary','middle')$cells,5L);expect_identical(row_of(g,'primary','high')$cells,2L)
  expect_equal(row_of(g,'primary','high')$signed_error,100*mean(c(-.13,-.13)),tolerance=1e-12)
  rare <- row_of(g,'primary','rare_1_5pct')
  expect_identical(rare$cells,8L);expect_equal(rare$signed_error,0,tolerance=1e-12)
  expect_equal(rare$mean_abs_cell_error,100*.111/8,tolerance=1e-12);expect_equal(rare$coverage,3/8)
  p1 <- row_of(g,'primary','prevalence_1pct')
  expect_identical(p1$cells,4L);expect_equal(p1$signed_error,100*.01/4,tolerance=1e-12);expect_equal(p1$coverage,.75)
  expect_equal(row_of(g,'primary','prevalence_75pct')$signed_error,100*mean(err[,4]),tolerance=1e-12)
  expect_equal(row_of(g,'primary','all')$mean_abs_cell_error,100*mean(abs(err)),tolerance=1e-12)
  expect_false('rare_below_20pct' %in% g$group)
})

test_that('coverage uses the type-7 0.025 and 0.975 quantiles of pooled draws, inclusive at both ends', {
  need_archives()
  # 41 draws in one chain: h = 40 * 0.025 + 1 = 2, so the lower bound is the second draw exactly.
  x <- array(c(seq(0,.4,by=.01),seq(0,.4,by=.01)),c(41,1,2));x <- aperm(x,c(3,1,2))
  truth <- c(.01,.011)
  b <- interval_cells(x,truth,sc)
  expect_equal(b$lower,rep(.01,2));expect_equal(b$upper,rep(.39,2))
  expect_identical(b$covered,c(TRUE,TRUE))
  expect_identical(interval_cells(x,c(.0099,.3901),sc)$covered,c(FALSE,FALSE))
  # pooled over chains: two chains of four draws each
  y <- array(c(1:4,5:8)/10,c(1,4,2))
  b <- interval_cells(y,.15,sc)
  expect_equal(b$lower,q7((1:8)/10,.025));expect_equal(b$upper,q7((1:8)/10,.975));expect_equal(b$estimate,.45)
  expect_identical(b$covered,TRUE)
})

test_that('B0 bias and coverage average over species on the logit scale', {
  # intervals contain the truth for species 1 and 2; species 3's truth 2 is above its upper bound 1.8
  b0 <- data.frame(species=1:3,truth=c(-1,0,2),estimate=c(-.5,.25,1),lower=c(-2,-1,.5),upper=c(0,1,1.8))
  b0$bias <- b0$estimate-b0$truth;b0$covered <- b0$lower<=b0$truth & b0$truth<=b0$upper
  s <- b0_summary(b0)
  expect_identical(s$species,3L);expect_equal(s$b0_bias,mean(c(.5,.25,-1)),tolerance=1e-12)
  expect_equal(s$b0_abs_bias,mean(c(.5,.25,1)),tolerance=1e-12);expect_equal(s$b0_coverage,2/3)
})

# ---------------------------------------------------------------------------
# Cells from saved draws, on synthetic fits
# ---------------------------------------------------------------------------

probs <- list(c(.05,.06,.07,.08,.09,.10,.11,.12),c(.40,.45,.50,.55,.60,.42,.48,.52),c(.95,.90,.85,.80,.97,.93,.88,.83))
# Two covariates and two site factors, all zero, as in the real A1 fits (one
# of either would make the archived scorer's array slices drop a dimension).
a1_fixture <- function(n=100L) {
  S <- 3L;ni <- 4L;nc <- 2L
  B0 <- array(NA_real_,c(S,ni,nc));for(s in 1:S) B0[s,,] <- matrix(qlogis(probs[[s]]),ni,nc)
  X <- cbind(seq(-1,1,length.out=n),rep(c(-1,1),length.out=n))
  js <- list(B0_output=B0,B_output=array(0,c(2,S,ni,nc)),U_output=array(0,c(n,2,ni,nc)),
    L_output=array(0,c(2,S,ni,nc)),sigmah_output=matrix(1,ni,nc))
  fit <- list(results_output=list(jsdm_output=js),infos=list(model='binary',ps=0L),X_psi=X)
  b0 <- qlogis(c(.04,.50,.90));B <- matrix(0,2,S);U <- matrix(0,n,2);L <- matrix(0,2,S)
  input <- list(n=n,data=list(OTU=matrix(0L,n,S,dimnames=list(NULL,paste0('sp',1:S)))),
    truth=list(B0=b0,B=B,U=U,L=L,eta=sweep(X%*%B+U%*%L,2,b0,'+')))
  list(fit=fit,input=input)
}

test_that('A1 cells: posterior-mean probability, truth and intervals from synthetic draws', {
  need_archives()
  x <- a1_fixture()
  expect_no_warning(cells <- a1_cells(x$fit,x$input,sc))
  expect_identical(cells$phase,'A1');expect_identical(c(cells$n,cells$S,cells$n_original),c(100L,3L,100L))
  expect_equal(cells$estimate[1,],c(mean(probs[[1]]),.49,.88875),tolerance=1e-12)
  expect_equal(cells$truth[1,],c(.04,.5,.9),tolerance=1e-12)
  expect_equal(cells$lower[7,],vapply(probs,q7,numeric(1),p=.025),tolerance=1e-12)
  expect_equal(cells$upper[7,],vapply(probs,q7,numeric(1),p=.975),tolerance=1e-12)
  expect_identical(cells$covered[1,],c(FALSE,TRUE,TRUE))   # .04 is below the lower bound .05175
  expect_lt(cells$checks$estimate_vs_score_fit,1e-12)
  expect_equal(cells$b0$bias,vapply(probs,function(p) mean(qlogis(p)),numeric(1))-qlogis(c(.04,.5,.9)),tolerance=1e-12)
  g <- group_table(cells,sc)
  expect_equal(row_of(g,'primary','low')$signed_error,100*(mean(probs[[1]])-.04),tolerance=1e-10)
  expect_identical(row_of(g,'primary','rare_below_20pct')$cells,100L)
  expect_equal(row_of(g,'primary','all')$coverage,2/3)
  # 120 fitted sites: original sites 1 to 100 are primary; the rare group uses all 120.
  y <- a1_fixture(120L);expect_no_warning(cy <- a1_cells(y$fit,y$input,sc));gy <- group_table(cy,sc)
  expect_identical(row_of(gy,'primary','all')$cells,300L);expect_identical(row_of(gy,'allsites','all')$cells,360L)
  expect_identical(row_of(gy,'primary','rare_below_20pct')$cells,120L)
})

a2_fixture <- function() {
  n <- 4L;S <- 4L;ps <- 2L;ni <- 4L;nc <- 2L
  p <- list(c(.009,.010,.011,.012,.013,.014,.015,.016),c(.03,.035,.04,.045,.05,.055,.06,.065),
    c(.15,.20,.25,.30,.35,.22,.28,.32),c(.70,.72,.74,.76,.78,.80,.82,.84))
  B0 <- array(NA_real_,c(S,ni,nc));for(s in 1:S) B0[s,,] <- matrix(qlogis(p[[s]]),ni,nc)
  xs <- matrix(c(0,1,0,1,0,0,1,1),4,2);X <- matrix(c(-1,0,1,2),n,1)
  js <- list(B0_output=B0,B_output=array(0,c(1,S,ni,nc)),Bs_output=array(0,c(ps,S,ni,nc)),
    idx_ls_output=matrix(1L,ni,nc),sigmabs_output=matrix(1,ni,nc))
  species <- paste0('sp',1:S)
  fit <- list(results_output=list(jsdm_output=js),X_psi=X,Xs=xs,
    infos=list(model='binary',ps=ps,n_factors=0L,speciesNames=species,l_s_grid=c(.5,1),
      list_Xs=list(X_tilde=xs[c(1,4),],X_s=xs,Xs_index=1:4)))
  truth <- list(psi=unclass(hand_a2()$truth),B0=qlogis(c(.012,.05,.3,.8)),X=X,Xs=xs,target_prevalence=c(.01,.05,.25,.75))
  input <- list(truth=truth,settings=list(n=n,S=S),data=list(binary=list(OTU=matrix(0L,n,S,dimnames=list(NULL,species)))))
  list(fit=fit,input=input,probs=p)
}

test_that('A2 cells: probability draws through the archived spatial reconstruction, on a synthetic fit', {
  need_archives()
  x <- a2_fixture()
  expect_no_warning(cells <- a2_cells(x$fit,x$input,sc,knots=2L))
  expect_identical(cells$phase,'A2');expect_identical(c(cells$n,cells$S,cells$n_original),c(4L,4L,4L))
  expect_equal(cells$estimate[2,],c(.0125,.0475,.25875,.77),tolerance=1e-12)
  expect_equal(cells$lower[3,],vapply(x$probs,q7,numeric(1),p=.025),tolerance=1e-12)
  expect_identical(cells$covered,unclass(hand_a2()$covered))
  expect_equal(cells$b0$bias,vapply(x$probs,function(p) mean(qlogis(p)),numeric(1))-qlogis(c(.012,.05,.3,.8)),tolerance=1e-12)
  expect_identical(cells$b0$covered,vapply(seq_along(x$probs),function(s) {
    v <- qlogis(x$probs[[s]]);t <- qlogis(c(.012,.05,.3,.8))[s];q7(v,.025)<=t && t<=q7(v,.975)},logical(1)))
  g <- group_table(cells,sc);h <- group_table(hand_a2(),sc)
  expect_equal(g[c('signed_error','mean_abs_cell_error','coverage')],h[c('signed_error','mean_abs_cell_error','coverage')],tolerance=1e-10)
  expect_error(a2_cells(x$fit,x$input,sc,knots=100L))
})

# ---------------------------------------------------------------------------
# The gate on synthetic community scores
# ---------------------------------------------------------------------------

# Per-fit primary rows for one arm: every community has the bands, the gated
# rare group and 'all' (whose mean absolute cell error is the overall MAE).
arm_rows <- function(phase,sd,values=list(),communities=NULL) {
  keys <- phase_keys(phase)
  base <- list(low=10,middle=5,high=5,rare_1_5pct=8,rare_below_20pct=20,all=0)
  cov <- list(low=.8,middle=.8,high=.8,rare_1_5pct=.8,rare_below_20pct=.8,all=.8)
  # A1 and B: bands and the rare group of replicate 07; A2: bands and the 1% plus 5% group.
  groups <- if(phase=='A2') c('low','middle','high','all','rare_1_5pct') else c('low','middle','high','all','rare_below_20pct')
  rows <- list()
  for(k in keys) for(g in groups) {
    if(g=='rare_below_20pct' && !grepl('-07$',k)) next
    rows[[length(rows)+1L]] <- data.frame(phase=phase,sd=sd,key=k,stratum=phase_stratum(phase,k),community=phase_community(phase,k),
      scope='primary',group=g,cells=10L,signed_error=base[[g]],abs_signed_error=base[[g]],mean_abs_cell_error=if(g=='all') 8 else base[[g]],
      coverage=cov[[g]],stringsAsFactors=FALSE)
  }
  x <- do.call(rbind,rows)
  for(v in values) {
    i <- x$group==v$group & (is.null(v$stratum) | x$stratum %in% (v$stratum %||% x$stratum)) &
      (is.null(v$community) | x$community %in% (v$community %||% x$community))
    x[i,v$column %||% 'abs_signed_error'] <- v$value
    if(is.null(v$column)) x$signed_error[i] <- v$value
  }
  x
}
conv_rows <- function(phase,sd,flagged=0L) {
  strata <- switch(phase,A1=c('n100','n300'),B=c('qnear','qfar'),'all');n <- if(phase=='A2') 9L else 10L
  data.frame(phase=phase,stratum=strata,sd=sd,role=if(sd==1) 'control' else 'new',selected_fits=n,
    selected_flagged=rep_len(flagged,length(strata)),stringsAsFactors=FALSE)
}
gate <- function(phase,new,control=arm_rows(phase,1),conv=rbind(conv_rows(phase,1),conv_rows(phase,2)),missing=NULL)
  gate_subphase(phase,rbind(control,new),conv,sd=2,missing=missing)
detail <- function(g,criterion,component=NULL,stratum=NULL) {
  d <- g$detail[g$detail$criterion==criterion,,drop=FALSE]
  if(!is.null(component)) d <- d[d$component==component,,drop=FALSE]
  if(!is.null(stratum)) d <- d[d$stratum==stratum,,drop=FALSE]
  d
}
v <- function(group,value,stratum=NULL,community=NULL,column=NULL) list(group=group,value=value,stratum=stratum,community=community,column=column)
reps <- sprintf('%02d',1:10)

test_that('A1 criterion 1 pools the n100 and n300 fits of each generating community (R21)', {
  # Communities 01-07: n100 6, n300 13, mean 9.5 < 10 improved, although n300 alone is worse.
  # Communities 08-10: n100 9, n300 11, mean exactly 10: a tie, not improved.
  new <- arm_rows('A1',2,list(v('low',6,'n100',reps[1:7]),v('low',13,'n300',reps[1:7]),
    v('low',9,'n100',reps[8:10]),v('low',11,'n300',reps[8:10])))
  g <- gate('A1',new)
  c1 <- detail(g,'c1','low')
  expect_identical(nrow(c1),1L);expect_identical(c1$stratum,'pooled')
  expect_identical(c1$communities,10L);expect_identical(c1$needed,7L);expect_identical(c1$improved,7L)
  expect_equal(c1$value_new,(7*9.5+3*10)/10);expect_equal(c1$value_control,10);expect_true(c1$pass)
  # six of ten is not enough
  new6 <- arm_rows('A1',2,list(v('low',6,'n100',reps[1:6]),v('low',13,'n300',reps[1:6]),
    v('low',9,'n100',reps[7:10]),v('low',11,'n300',reps[7:10])))
  c16 <- detail(gate('A1',new6),'c1','low')
  expect_identical(c16$improved,6L);expect_false(c16$pass)
  # the rare group is not gated in A1, however it moves
  worse_rare <- arm_rows('A1',2,list(v('low',9),v('rare_below_20pct',90)))
  expect_false('rare_below_20pct' %in% gate('A1',worse_rare)$detail$component)
  expect_true(detail(gate('A1',worse_rare),'c1','low')$pass)
})

test_that('criterion 1 also needs a strictly lower across-community mean', {
  # nine of ten improved by 1, one community 12 points worse: the mean rises from 10 to 10.3
  new <- arm_rows('A1',2,list(v('low',9),v('low',22,community='10')))
  c1 <- detail(gate('A1',new),'c1','low')
  expect_identical(c1$improved,9L);expect_equal(c1$value_new,10.3);expect_false(c1$pass)
  # equal means are not strictly lower
  same <- detail(gate('A1',arm_rows('A1',2)),'c1','low')
  expect_identical(same$improved,0L);expect_false(same$pass)
})

test_that('A2 criterion 1 gates the low band and the 1% and 5% rare group, 6 of 9 each', {
  # one point better than the control (low 10, rare group 8) in the first k communities
  improved <- function(group,k) v(group,if(group=='low') 9 else 7,community=phase_keys('A2')[seq_len(k)])
  g <- gate('A2',arm_rows('A2',2,list(improved('low',6),improved('rare_1_5pct',6))))
  expect_setequal(detail(g,'c1')$component,c('low','rare_1_5pct'))
  expect_identical(detail(g,'c1','low')$needed,6L);expect_true(all(detail(g,'c1')$pass))
  expect_true(g$row$c1_pass);expect_identical(g$row$result,'PASS')
  g5 <- gate('A2',arm_rows('A2',2,list(improved('low',6),improved('rare_1_5pct',5))))
  expect_false(detail(g5,'c1','rare_1_5pct')$pass);expect_true(detail(g5,'c1','low')$pass)
  expect_identical(g5$row$result,'FAIL')
})

test_that('criterion 2: no band worse by more than 1 point, overall MAE by more than 0.5, at each A1 size', {
  ok <- list(v('low',9))
  # middle band exactly 1 point worse at n300 passes; slightly more fails
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('middle',6,'n300')))))
  d <- detail(g,'c2','middle','n300');expect_equal(d$value_new-d$value_control,1);expect_true(d$pass)
  expect_false(d$tolerance_decisive);expect_true(g$row$c2_pass)   # exactly 1 in doubles: no allowance needed
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('middle',6.001,'n300')))))
  expect_false(detail(g,'c2','middle','n300')$pass);expect_true(detail(g,'c2','middle','n100')$pass)
  expect_false(g$row$c2_pass);expect_identical(g$row$result,'FAIL')
  # overall MAE: +0.5 passes, +0.51 fails, judged at n100 and n300 separately
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('all',8.5,'n100',column='mean_abs_cell_error')))))
  expect_true(detail(g,'c2','mae','n100')$pass);expect_true(g$row$c2_pass)
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('all',8.51,'n100',column='mean_abs_cell_error')))))
  expect_false(detail(g,'c2','mae','n100')$pass);expect_false(g$row$c2_pass)
  # a better band never fails
  expect_true(detail(gate('A1',arm_rows('A1',2,c(ok,list(v('high',1))))),'c2','high','n100')$pass)
})

test_that('criterion 3: band coverage may fall by at most 0.03, at each A1 size', {
  ok <- list(v('low',9))
  base <- arm_rows('A1',1);base$coverage[base$group=='middle'] <- .53
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('middle',.50,column='coverage')))),control=base)
  d <- detail(g,'c3','middle','n100')
  expect_equal(d$value_control-d$value_new,.03);expect_true(d$pass);expect_true(g$row$c3_pass)
  # .53 - .50 is 0.030000000000000027 in doubles: GATE_TOL decides it, and gate-detail says so
  expect_gt(d$statistic,.03);expect_true(d$tolerance_decisive)
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('middle',.4999,column='coverage')))),control=base)
  expect_false(detail(g,'c3','middle','n100')$pass);expect_false(detail(g,'c3','middle','n100')$tolerance_decisive)
  expect_false(g$row$c3_pass);expect_identical(g$row$result,'FAIL')
  g <- gate('A2',arm_rows('A2',2,list(v('low',9),v('rare_1_5pct',7),v('high',.7,column='coverage'))))
  expect_false(detail(g,'c3','high','all')$pass);expect_identical(g$row$result,'FAIL')
})

test_that('criterion 4 is strict (R20): any flagged selected new-arm fit fails when no control is flagged', {
  ok <- arm_rows('A1',2,list(v('low',9)))
  g <- gate('A1',ok);expect_identical(g$row$result,'PASS');expect_true(all(detail(g,'c4')$pass))
  flagged <- rbind(conv_rows('A1',1),transform(conv_rows('A1',2),selected_flagged=c(0L,1L)))
  g <- gate('A1',ok,conv=flagged)
  expect_true(detail(g,'c4',stratum='n100')$pass);expect_false(detail(g,'c4',stratum='n300')$pass)
  expect_identical(detail(g,'c4',stratum='n300')$value_new,1);expect_identical(g$row$result,'FAIL')
  # as many flags as the control is allowed
  both <- rbind(transform(conv_rows('A2',1),selected_flagged=1L),transform(conv_rows('A2',2),selected_flagged=1L))
  expect_true(detail(gate('A2',arm_rows('A2',2,list(v('low',9),v('rare_1_5pct',7))),conv=both),'c4')$pass)
})

test_that('a community without a valid selected fit fails criterion 4; criteria 1 to 3 use the rest (R27)', {
  ok <- arm_rows('A1',2,list(v('low',9)))
  gone <- ok[ok$key!='jsdm-n0300-04',]
  conv <- rbind(conv_rows('A1',1),transform(conv_rows('A1',2),selected_fits=c(10L,9L)))
  g <- gate('A1',gone,conv=conv)
  expect_false(detail(g,'c4',stratum='n300')$pass);expect_identical(detail(g,'c4',stratum='n300')$missing_new,1L)
  expect_true(detail(g,'c4',stratum='n100')$pass);expect_identical(g$row$result,'FAIL')
  # pooled criterion 1 drops the whole generating community 04: 9 communities, 6 needed
  c1 <- detail(g,'c1','low');expect_identical(c1$communities,9L);expect_identical(c1$needed,6L);expect_true(c1$pass)
  expect_identical(detail(g,'c2','low','n300')$communities,9L);expect_identical(detail(g,'c2','low','n100')$communities,10L)
  # a manual record of a missing fit (AMENDMENT-2) for a community that has no score has the same effect
  g <- gate('A1',ok[ok$key!='jsdm-n0100-09',],missing=data.frame(sd=2,key='jsdm-n0100-09',reason='quarantined'))
  expect_false(detail(g,'c4',stratum='n100')$pass);expect_identical(detail(g,'c1','low')$communities,9L)
  expect_identical(g$row$result,'FAIL')
  # a community both scored and recorded as missing is refused, as analysis.R and supplement.R refuse it
  expect_error(gate('A1',ok,missing=data.frame(sd=2,key='jsdm-n0100-09',reason='quarantined')),'both scored and recorded as missing')
  # the control arm must be complete
  expect_error(gate('A1',ok,control=arm_rows('A1',1)[arm_rows('A1',1)$key!='jsdm-n0100-01',]),'control')
})

test_that('a criterion over no community is NA, never PASS; the result is then FAIL or UNDETERMINED', {
  # R27 as select.R can leave it: the SD has no selected fit at all, and no convergence count
  conv <- rbind(conv_rows('A1',1),data.frame(phase='A1',stratum=c('n100','n300'),sd=2,role='new',selected_fits=0L,selected_flagged=NA_integer_))
  g <- gate('A1',arm_rows('A1',2)[0,],conv=conv)
  expect_true(all(is.na(g$detail$pass[g$detail$criterion %in% c('c1','c2','c3')])))
  expect_identical(detail(g,'c1','low')$communities,0L);expect_identical(detail(g,'c4',stratum='n100')$missing_new,10L)
  expect_false(any(detail(g,'c4')$pass));expect_true(is.na(g$row$c1_pass) && is.na(g$row$c2_pass) && is.na(g$row$c3_pass))
  expect_false(g$row$c4_pass);expect_identical(g$row$result,'FAIL')
  # a band with no cells in any community, all else passing: UNDETERMINED, not PASS
  nohigh <- function(x) x[x$group!='high',]
  g <- gate('A1',nohigh(arm_rows('A1',2,list(v('low',9)))),control=nohigh(arm_rows('A1',1)))
  expect_true(is.na(detail(g,'c2','high','n100')$pass));expect_true(is.na(g$row$c2_pass));expect_identical(g$row$result,'UNDETERMINED')
})

test_that('a band empty in a community is left out for that community', {
  new <- arm_rows('A2',2,list(v('low',9),v('rare_1_5pct',7)))
  ctrl <- arm_rows('A2',1)
  drop <- function(x) x[!(x$group=='high' & x$key=='range8-rep03-binary-k100'),]
  g <- gate('A2',drop(new),control=drop(ctrl))
  expect_identical(detail(g,'c2','high','all')$communities,8L);expect_identical(detail(g,'c3','high','all')$communities,8L)
  expect_identical(detail(g,'c2','low','all')$communities,9L);expect_true(detail(g,'c4')$pass)
})

test_that('phase A passes only for an SD that passes A1 and A2; the smallest passing SD is marked', {
  a1 <- data.frame(sd=c(2,3,5),subphase='A1',result=c('PASS','PASS','FAIL'),stringsAsFactors=FALSE)
  a2 <- data.frame(sd=c(2,3,5),subphase='A2',result=c('FAIL','PASS','PASS'),stringsAsFactors=FALSE)
  x <- combine_phase_a(a1,a2);a <- x[x$subphase=='A',]
  expect_identical(a$sd,c(2,3,5));expect_identical(a$result,c('FAIL','PASS','FAIL'))
  expect_identical(a$smallest_passing,c(FALSE,TRUE,FALSE));expect_identical(nrow(x),9L)
  expect_error(combine_phase_a(a1,a2[a2$sd!=5,]),'both')
})

# ---------------------------------------------------------------------------
# Real fits as mechanics fixtures
# ---------------------------------------------------------------------------

test_that('an archived control is scored exactly as its archived result (A1 jsdm-n0100-01)', {
  need_archives()
  items <- control_items('A1',study,archives,repo)
  item <- items[items$key=='jsdm-n0100-01',,drop=FALSE]
  out <- tempfile('scores')
  path <- score_item(item,sc,archives,inputs_root,out)
  r <- readRDS(path);old <- readRDS(file.path(archives,'pr11-current-20260927/initial/jsdm-n0100-01-result.rds'))
  expect_identical(r$fit_md5,item$expected_md5);expect_identical(r$role,'control');expect_identical(r$sd,1)
  for(b in c('all','low','middle','high')) {
    mine <- row_of(r$groups,'primary',b);theirs <- old$groups[old$groups$scope=='original100' & old$groups$group==b,]
    expect_identical(mine$cells,as.integer(theirs$cells))
    expect_lt(abs(mine$signed_error-100*theirs$bias),1e-10);expect_lt(abs(mine$mean_abs_cell_error-100*theirs$mae),1e-10)
  }
  expect_lt(max(abs(r$cells$estimate-as.vector(old$estimate))),1e-12)
  # resuming leaves the record as it is
  mtime <- file.mtime(path);expect_identical(score_item(item,sc,archives,inputs_root,out),path)
  expect_identical(file.mtime(path),mtime)
})

test_that('the new-arm loading path checks the fit it scores (sd 3 pilot fit, mechanics only)', {
  need_archives()
  pilot <- file.path(study,'fits/A1/sd3/pilot/jsdm-n0100-01-fit.rds');if(!file.exists(pilot)) skip('pilot fit not available')
  item <- data.frame(role='new',phase='A1',sd=3,key='jsdm-n0100-01',schedule='pilot',fit=pilot,
    fit_label='fits/A1/sd3/pilot/jsdm-n0100-01-fit.rds',expected_md5=unname(tools::md5sum(pilot)),kind='selected',stringsAsFactors=FALSE)
  out <- tempfile('scores')
  expect_no_warning(p1 <- score_item(item,sc,archives,inputs_root,out));r <- readRDS(p1)
  # Structure only: no outcome of the pilot fit is compared or printed.
  expect_identical(nrow(r$cells),1000L);expect_true(all(c('primary') %in% r$groups$scope))
  expect_true(all(is.finite(r$groups$signed_error)));expect_identical(r$schedule,'pilot')
  bad <- item;bad$sd <- 2
  expect_error(score_item(bad,sc,archives,inputs_root,tempfile('scores')),'not the expected fit')
  bad <- item;bad$expected_md5 <- strrep('0',32)
  expect_error(score_item(bad,sc,archives,inputs_root,tempfile('scores')),'md5')
  # The A2 path, on the A2 sd 3 pilot fit (200 draws per chain), again structure only.
  pilot2 <- file.path(study,'fits/A2/sd3/pilot/range6-rep01-binary-k100-fit.rds');if(!file.exists(pilot2)) skip('A2 pilot fit not available')
  item2 <- data.frame(role='new',phase='A2',sd=3,key='range6-rep01-binary-k100',schedule='pilot',fit=pilot2,
    fit_label='fits/A2/sd3/pilot/range6-rep01-binary-k100-fit.rds',expected_md5=unname(tools::md5sum(pilot2)),kind='selected',stringsAsFactors=FALSE)
  expect_no_warning(p2 <- score_item(item2,sc,archives,inputs_root,out));r2 <- readRDS(p2)
  expect_identical(nrow(r2$cells),800L);expect_true(all(c('rare_1_5pct','prevalence_1pct') %in% r2$groups$group))
  expect_identical(r2$b0$species,8L);expect_true(all(is.finite(r2$groups$coverage)))
})

# ---------------------------------------------------------------------------
# Command lines and tables
# ---------------------------------------------------------------------------

cli_args <- function(study_dir,...) c(paste0('--repo=',repo),paste0('--study=',study_dir),paste0('--inputs-root=',inputs_root),
  paste0('--archives=',archives),...)

test_that('analysis.R refuses new arms before the selections exist, and invalid arms', {
  need_archives()
  empty <- tempfile('study');dir.create(empty)
  r <- run_cli('analysis.R',cli_args(empty,'--phase=A2','--arms=2',paste0('--out=',tempfile())))
  expect_identical(r$status,1L);expect_match(r$output,'selected-fits.csv')
  # A1 selected but A2 not: AMENDMENT-1 forbids new-arm errors before both selections exist
  dir.create(file.path(empty,'selection/A1'),recursive=TRUE);writeLines('x',file.path(empty,'selection/A1/selected-fits.csv'))
  r <- run_cli('analysis.R',cli_args(empty,'--phase=A1','--arms=2',paste0('--out=',tempfile())))
  expect_identical(r$status,1L);expect_match(r$output,'AMENDMENT-1')
  r <- run_cli('analysis.R',cli_args(empty,'--phase=A1','--arms=4'))
  expect_identical(r$status,1L);expect_match(r$output,'--arms')
  r <- run_cli('analysis.R',cli_args(empty,'--phase=B','--arms=1'))
  expect_identical(r$status,1L);expect_match(r$output,'phase')
  expect_identical(list.files(empty,recursive=TRUE),'selection/A1/selected-fits.csv')
})

# A temporary study whose score records are synthetic, built with the scorer's
# own record constructor, for summarise.R end to end. new_sds: SDs recorded by
# select.R (the first with flag_n300 flagged 300-site fits); absent_sd: an SD
# select.R could not record, with one manual-missing community and nothing
# else; supplement_sd: the same with the manual record of its other 19 fits,
# processed by supplement.R's rule (AMENDMENT-2, R27).
synthetic_study <- function(new_sds=NULL,flag_n300=0L,absent_sd=NULL,supplement_sd=NULL,gone='jsdm-n0100-03') {
  st <- tempfile('study');dir.create(st);out <- file.path(st,'summary')
  items <- control_items('A1',st,archives,repo)
  hashes <- current_score_hashes(repo,archives)
  write_record <- function(item,cells) {
    rec <- score_record(item,cells,sc,hashes,input_md5='synthetic');p <- score_path(out,'A1',item$sd,item$key,item$schedule)
    dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);saveRDS(rec,p)
  }
  base <- hand_a1()
  for(i in seq_len(nrow(items))) write_record(items[i,,drop=FALSE],base)
  if(length(c(new_sds,absent_sd,supplement_sd))) {
    better <- base;better$estimate <- better$truth+(base$estimate-base$truth)/2
    mk <- function(sd) {x <- items;x$role <- 'new';x$sd <- sd;x$schedule <- 'initial'
      x$fit_label <- sprintf('fits/A1/sd%d/initial/%s-fit.rds',sd,x$key);x$fit <- file.path(st,x$fit_label)
      x$expected_md5 <- sprintf('%02d%030d',sd,seq_len(nrow(x)));x}
    d <- file.path(st,'selection/A1');dir.create(d,recursive=TRUE)
    cf <- data.frame(role='control',phase='A1',sd=1,key=items$key,schedule=items$schedule,fit=items$fit_label,fit_md5=items$expected_md5,
      rule='pr11',warnings=0L,max_group_rhat=1,max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA,spatial_field_flags=NA,
      flagged=FALSE,reasons='',stringsAsFactors=FALSE)
    utils::write.csv(cf,file.path(d,'control-flags.csv'),row.names=FALSE)
    sf <- NULL
    for(sd in new_sds) {
      sel <- mk(sd);for(i in seq_len(nrow(sel))) write_record(sel[i,,drop=FALSE],better)
      flag <- sd==new_sds[1] & sel$key %in% sprintf('jsdm-n0300-%02d',seq_len(flag_n300))
      sf <- rbind(sf,data.frame(phase='A1',sd=sd,key=sel$key,role='new',first_schedule='initial',selected_schedule='initial',long_repeat=FALSE,
        fit=sel$fit_label,fit_md5=sel$expected_md5,rule='pr11',first_flagged=flag,flagged=flag,reasons=ifelse(flag,'stub',''),stringsAsFactors=FALSE))
    }
    sf <- rbind(sf,data.frame(phase='A1',sd=1,key=items$key,role='control',first_schedule=NA,selected_schedule=items$schedule,long_repeat=NA,
      fit=items$fit_label,fit_md5=items$expected_md5,rule='pr11',first_flagged=NA,flagged=FALSE,reasons='',stringsAsFactors=FALSE))
    utils::write.csv(sf,file.path(d,'selected-fits.csv'),row.names=FALSE)
    utils::write.csv(convergence_counts(sf),file.path(d,'convergence.csv'),row.names=FALSE)
    manual <- c(absent_sd,supplement_sd)
    if(length(manual)) utils::write.csv(data.frame(phase='A1',sd=manual,key=gone,reason='quarantined, see its log'),
      file.path(d,'manual-missing.csv'),row.names=FALSE)
    if(!is.null(supplement_sd)) {
      sel <- mk(supplement_sd);sel <- sel[sel$key!=gone,]
      ms <- data.frame(phase='A1',sd=supplement_sd,key=sel$key,schedule='initial',fit=sel$fit_label,fit_md5=sel$expected_md5,stringsAsFactors=FALSE)
      utils::write.csv(ms,file.path(d,'manual-selected.csv'),row.names=FALSE)
      flags <- data.frame(role='new',phase='A1',sd=supplement_sd,key=sel$key,schedule='initial',fit=sel$fit_label,fit_md5=sel$expected_md5,
        rule='pr11',warnings=0L,max_group_rhat=1,max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA,spatial_field_flags=NA,
        flagged=FALSE,reasons='',stringsAsFactors=FALSE)
      ms2 <- read_manual_selected(d,'A1');check_supplement('A1',ms2,read_manual_missing(d,'A1'),unique(sf$sd),control_schedule_of('A1',archives))
      tab <- supplement_table('A1',ms2,flags,control_schedule_of('A1',archives))
      dir.create(file.path(out,'A1'),recursive=TRUE,showWarnings=FALSE)
      write_frozen_table(flags,file.path(out,'A1','manual-flags.csv'));write_frozen_table(tab,file.path(out,'A1','manual-selection.csv'))
      for(i in seq_len(nrow(sel))) write_record(sel[i,,drop=FALSE],better)
    }
  }
  list(study=st,out=out)
}

test_that('summarise.R with --arms=1 writes the control tables and no gate', {
  need_archives()
  s <- synthetic_study()
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1',paste0('--out=',s$out)))
  expect_identical(r$status,0L)
  d <- file.path(s$out,'A1')
  for(f in c('band-error.csv','coverage.csv','mae.csv','b0.csv','convergence.csv','summary-means.csv','b0-means.csv','provenance.csv'))
    expect_true(file.exists(file.path(d,f)),info=f)
  expect_false(file.exists(file.path(d,'gate.csv')))
  be <- utils::read.csv(file.path(d,'band-error.csv'),stringsAsFactors=FALSE)
  expect_setequal(unique(be$sd),1);expect_identical(length(unique(be$key)),20L)
  low <- be[be$key=='jsdm-n0100-01' & be$scope=='primary' & be$group=='low',]
  expect_equal(low$signed_error,5,tolerance=1e-12)
  cv <- utils::read.csv(file.path(d,'convergence.csv'),stringsAsFactors=FALSE)
  expect_identical(cv$selected_fits,c(10L,10L));expect_identical(cv$selected_flagged,c(0L,0L))
  mae <- utils::read.csv(file.path(d,'mae.csv'),stringsAsFactors=FALSE)
  expect_equal(unique(mae$mae[mae$scope=='primary']),100*.40/6,tolerance=1e-12)
})

test_that('summarise.R gates a new arm end to end, including R20, and prints no verdict', {
  need_archives()
  s <- synthetic_study(new_sds=2)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))
  expect_identical(r$status,0L);expect_false(grepl('PASS|FAIL',r$output))   # the verdict is read only after verify.R
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE)
  # the new arm halves every cell error: all bands better, coverage unchanged, no flags
  expect_identical(nrow(g),1L);expect_identical(g$sd,2L);expect_identical(g$result,'PASS')
  expect_identical(g$c1_low_improved,10L);expect_identical(g$c1_low_needed,7L)
  gd <- utils::read.csv(file.path(s$out,'A1/gate-detail.csv'),stringsAsFactors=FALSE)
  expect_true('tolerance_decisive' %in% names(gd));expect_false(any(gd$tolerance_decisive,na.rm=TRUE))
  # one flagged selected fit at 300 sites fails the gate (R20)
  s <- synthetic_study(new_sds=2,flag_n300=1L)
  expect_identical(run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))$status,0L)
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE)
  expect_identical(g$result,'FAIL');expect_identical(g$c4_n300_flagged_new,1L);expect_true(g$c1_pass)
  # a score record that does not match the selected fit is refused
  s <- synthetic_study(new_sds=2)
  p <- score_path(s$out,'A1',2,'jsdm-n0100-05','initial');x <- readRDS(p);x$fit_md5 <- strrep('f',32);saveRDS(x,p)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))
  expect_identical(r$status,1L);expect_match(r$output,'jsdm-n0100-05')
})

test_that('the reachable R27 state: an SD that select.R could not record is reported by the manual record (AMENDMENT-2)', {
  need_archives()
  # SD 2 is absent from selected-fits.csv and has only a manual-missing row: criteria 1-3 NA, criterion 4 FAIL
  s <- synthetic_study(new_sds=3,absent_sd=2)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2,3',paste0('--out=',s$out)))
  expect_identical(r$status,0L)
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE);g2 <- g[g$sd==2,];g3 <- g[g$sd==3,]
  expect_identical(g2$result,'FAIL');expect_true(is.na(g2$c1_pass) && is.na(g2$c2_pass) && is.na(g2$c3_pass))
  expect_false(g2$c4_pass);expect_identical(c(g2$c4_n100_missing_new,g2$c4_n300_missing_new),c(10L,10L));expect_identical(g3$result,'PASS')
  cv <- utils::read.csv(file.path(s$out,'A1/convergence.csv'),stringsAsFactors=FALSE)
  expect_identical(unique(cv$source[cv$sd==2]),'absent (R27)');expect_identical(cv$selected_fits[cv$sd==2],c(0L,0L))
  # with the manual record of its other 19 fits: criteria 1-3 reported over the valid communities, criterion 4 FAIL
  s <- synthetic_study(new_sds=3,supplement_sd=2)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2,3',paste0('--out=',s$out)))
  expect_identical(r$status,0L)
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE);g2 <- g[g$sd==2,]
  expect_identical(g2$result,'FAIL');expect_identical(g2$c1_low_communities,9L);expect_true(g2$c1_pass)
  expect_true(g2$c2_pass);expect_true(g2$c3_pass);expect_false(g2$c4_n100_pass);expect_true(g2$c4_n300_pass)
  expect_identical(g2$c4_n100_missing_new,1L)
  cv <- utils::read.csv(file.path(s$out,'A1/convergence.csv'),stringsAsFactors=FALSE)
  expect_identical(unique(cv$source[cv$sd==2]),'manual supplement (R27)');expect_identical(cv$selected_fits[cv$sd==2],c(9L,10L))
  be <- utils::read.csv(file.path(s$out,'A1/band-error.csv'),stringsAsFactors=FALSE)
  expect_identical(length(unique(be$key[be$sd==2])),19L)
  # the supplement must list exactly the controller's manual record
  f <- file.path(s$study,'selection/A1/manual-selected.csv');m <- utils::read.csv(f,stringsAsFactors=FALSE,colClasses=c(fit_md5='character'));utils::write.csv(m[-1,],f,row.names=FALSE)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2,3',paste0('--out=',s$out)))
  expect_identical(r$status,1L);expect_match(r$output,'manual-selected.csv')
})

test_that('supplement.R enforces the manual-record rules and the selection rule with the frozen flags', {
  need_archives()
  cs <- control_schedule_of('A1',archives);keys <- phase_keys('A1')
  mm <- data.frame(phase='A1',sd=2,key='jsdm-n0100-03',reason='quarantined',stringsAsFactors=FALSE)
  ms <- data.frame(phase='A1',sd=2,key=setdiff(keys,'jsdm-n0100-03'),schedule='initial',stringsAsFactors=FALSE)
  ms$fit <- sprintf('fits/A1/sd2/initial/%s-fit.rds',ms$key);ms$fit_md5 <- sprintf('%032d',seq_len(nrow(ms)))
  expect_silent(check_supplement('A1',ms,mm,c(1,3),cs))
  expect_error(check_supplement('A1',ms,mm,c(1,2,3),cs),'in selected-fits.csv')
  expect_error(check_supplement('A1',ms,mm[0,],c(1,3),cs),'no manual-missing')
  expect_error(check_supplement('A1',ms[-1,],mm,c(1,3),cs),'exactly once')
  expect_error(check_supplement('A1',rbind(ms,transform(ms[1,],key='jsdm-n0100-03')),mm,c(1,3),cs),'both selected')
  a2 <- control_schedule_of('A2',archives);k2 <- phase_keys('A2')
  ms2 <- data.frame(phase='A2',sd=2,key=k2[-1],schedule='initial',fit=sprintf('fits/A2/sd2/initial/%s-fit.rds',k2[-1]),
    fit_md5=sprintf('%032d',1:8),stringsAsFactors=FALSE)
  expect_error(check_supplement('A2',ms2,data.frame(sd=2,key=k2[1]),numeric(),a2),'longer schedule only')
  flag <- function(x,flagged) data.frame(role='new',phase='A1',sd=x$sd,key=x$key,schedule=x$schedule,fit=x$fit,fit_md5=x$fit_md5,rule='pr11',
    warnings=0L,max_group_rhat=1,max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA,spatial_field_flags=NA,flagged=flagged,
    reasons=ifelse(flagged,'stub',''),stringsAsFactors=FALSE)
  t <- supplement_table('A1',ms,flag(ms,FALSE),cs)
  expect_identical(nrow(t),19L);expect_false(any(t$long_repeat));expect_identical(t$first_schedule,rep('initial',19L))
  expect_error(supplement_table('A1',ms,flag(ms,ms$key=='jsdm-n0100-05'),cs),'single longer repeat')
  # a longer repeat is selected only for a flagged first fit, and the first fit is recorded
  rep <- ms;rep$schedule[1] <- 'long';rep$fit[1] <- sub('/initial/','/long/',rep$fit[1])
  items <- supplement_flag_items('A1',rep,'/study',cs)
  expect_identical(sum(items$kind=='first'),1L);expect_identical(items$schedule[items$kind=='first'],'initial')
  first <- flag(transform(rep[1,],schedule='initial',fit=sub('/long/','/initial/',rep$fit[1]),fit_md5=strrep('a',32)),TRUE)
  t <- supplement_table('A1',rep,rbind(flag(rep,FALSE),first),cs)
  expect_true(t$long_repeat[t$key==rep$key[1]]);expect_true(t$first_flagged[t$key==rep$key[1]])
  expect_identical(t$first_fit_md5[t$key==rep$key[1]],strrep('a',32))
  first$flagged <- FALSE
  expect_error(supplement_table('A1',rep,rbind(flag(rep,FALSE),first),cs),'no flagged first fit')
  # like analysis.R, it refuses to score before the other sub-phase's selection exists, and without a record
  empty <- tempfile('study');dir.create(file.path(empty,'selection/A1'),recursive=TRUE)
  r <- run_cli('supplement.R',cli_args(empty,'--phase=A1'))
  expect_identical(r$status,1L);expect_match(r$output,'AMENDMENT-1')
  dir.create(file.path(empty,'selection/A2'));writeLines('x',file.path(empty,'selection/A2/selected-fits.csv'))
  r <- run_cli('supplement.R',cli_args(empty,'--phase=A1'))
  expect_identical(r$status,1L);expect_match(r$output,'No manual record')
  expect_false(dir.exists(file.path(empty,'summary')))
  # the manual record must name the protocol path of each fit
  d <- tempfile('sel');dir.create(d);bad <- ms;bad$fit[1] <- 'elsewhere.rds';utils::write.csv(bad,file.path(d,'manual-selected.csv'),row.names=FALSE)
  expect_error(read_manual_selected(d,'A1'),'Malformed')
})

# ---------------------------------------------------------------------------
# The independent audit
# ---------------------------------------------------------------------------

test_that('verify.R reimplements the metrics independently and agrees on the synthetic fits', {
  need_archives()
  expect_false(any(grepl('analysis.R',readLines(file.path(here,'verify.R')),fixed=TRUE) &
    grepl('source',readLines(file.path(here,'verify.R')),fixed=TRUE)))
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  x <- a1_fixture(120L)
  mine <- ve$audit_a1_cells(x$fit,x$input)
  expect_no_warning(cx <- a1_cells(x$fit,x$input,sc));theirs <- group_table(cx,sc)
  m <- merge(ve$audit_groups(mine),theirs,by=c('scope','group'))
  expect_identical(nrow(m),nrow(theirs))
  expect_lt(max(abs(m$signed_error.x-m$signed_error.y)),1e-10);expect_lt(max(abs(m$coverage.x-m$coverage.y)),1e-12)
  y <- a2_fixture()
  mine <- ve$audit_groups(ve$audit_cells_from_draws('A2',ve$audit_a2_draws(y$fit),y$input$truth$psi,y$input$truth$target_prevalence))
  theirs <- group_table(hand_a2(),sc)
  m <- merge(mine,theirs,by=c('scope','group'))
  expect_identical(nrow(m),nrow(theirs));expect_lt(max(abs(m$signed_error.x-m$signed_error.y)),1e-10)
  expect_lt(max(abs(m$coverage.x-m$coverage.y)),1e-12)
})

# summarise.R's tables of a synthetic sub-phase and this audit's own values for
# the same communities, for the table and gate audit.
gate_fixture <- function(phase,new,control=arm_rows(phase,1),flagged_new=character(),missing=NULL) {
  x <- rbind(control,new);x$arm <- ifelse(x$sd==1,'control','new');x$kind <- 'selected'
  k <- unique(x[c('sd','key','stratum')])
  flags <- data.frame(sd=k$sd,key=k$key,stratum=k$stratum,flagged=k$sd!=1 & k$key %in% flagged_new,stringsAsFactors=FALSE)
  strata <- unique(phase_stratum(phase,phase_keys(phase)))
  conv <- do.call(rbind,lapply(c(1,2),function(a) do.call(rbind,lapply(strata,function(h) {f <- flags[flags$sd==a & flags$stratum==h,];n <- nrow(f)
    data.frame(phase=phase,stratum=h,sd=a,role=if(a==1) 'control' else 'new',selected_fits=n,selected_flagged=if(n) sum(f$flagged) else NA,
      first_flagged=if(a==1 || !n) NA else sum(f$flagged),long_repeats=if(a==1 || !n) NA else 0L,
      missing_fits=sum(phase_stratum(phase,phase_keys(phase))==h)-n,stringsAsFactors=FALSE)}))))
  g <- gate_subphase(phase,x,conv,2,missing)
  m <- gate_subphase(phase,x,NULL,2,missing)$detail;m$first_fits_used <- 0L
  b0 <- unique(x[c('phase','sd','arm','kind','key','stratum','community')])
  b0$species <- 10L;b0$b0_bias <- ifelse(b0$sd==1,.3,.1);b0$b0_abs_bias <- .4;b0$b0_coverage <- .9
  list(own=cbind(x,flagged=flags$flagged[match(paste(x$sd,x$key),paste(flags$sd,flags$key))]),own_b0=b0,flags=flags,
    tabs=list(means=summary_means(x),b0means=b0_means(b0),convergence=conv,detail=g$detail,gate=g$row,matched=m),gate=g)
}
failing_tables <- function(ag) unique(ag$checks$table[!ag$checks$pass])

test_that('verify.R recomputes the across-community tables and every gate criterion, and catches planted errors', {
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  audit <- function(f,tabs=f$tabs,phase='A1') ve$audit_gate_tables(phase,c(1,2),f$own,f$own_b0,f$flags,NULL,tabs)
  # the R21 pooling case: 7 of 10 improved when pooled; ties; a middle-band coverage drop decided by GATE_TOL
  base <- arm_rows('A1',1);base$coverage[base$group=='middle'] <- .53
  new <- arm_rows('A1',2,list(v('low',6,'n100',reps[1:7]),v('low',13,'n300',reps[1:7]),v('low',9,'n100',reps[8:10]),
    v('low',11,'n300',reps[8:10]),v('middle',.50,column='coverage')))
  f <- gate_fixture('A1',new,control=base)
  ag <- audit(f);expect_true(all(ag$checks$pass));expect_identical(unique(ag$gate$result),f$gate$row$result)
  expect_true(all(c('summary-means.csv','b0-means.csv','convergence.csv','gate-detail.csv','gate.csv','schedule-matched.csv') %in% ag$checks$table))
  expect_true(ag$gate$tolerance_decisive[ag$gate$criterion=='c3' & ag$gate$component=='middle' & ag$gate$stratum=='n100'])
  plant <- function(part,fn) {t <- f$tabs;t[[part]] <- fn(t[[part]]);failing_tables(audit(f,t))}
  expect_identical(plant('gate',function(x) transform(x,c1_low_improved=c1_low_improved-1L)),'gate.csv')
  expect_identical(plant('gate',function(x) transform(x,result='PASS',c2_pass=TRUE)),'gate.csv')
  expect_identical(plant('detail',function(x) {i <- x$criterion=='c2' & x$component=='middle' & x$stratum=='n300';x$value_new[i] <- x$value_new[i]+1e-6;x}),'gate-detail.csv')
  expect_identical(plant('detail',function(x) {i <- x$criterion=='c1';x$pass[i] <- !x$pass[i];x}),'gate-detail.csv')
  expect_identical(plant('detail',function(x) {i <- x$criterion=='c3' & x$component=='middle';x$tolerance_decisive[i] <- FALSE;x}),'gate-detail.csv')
  expect_identical(plant('means',function(x) {x$mean_coverage[1] <- x$mean_coverage[1]+1e-6;x}),'summary-means.csv')
  expect_identical(plant('b0means',function(x) {x$mean_b0_bias[2] <- 0;x}),'b0-means.csv')
  expect_identical(plant('convergence',function(x) {x$selected_flagged[x$sd==2][1] <- 1L;x}),'convergence.csv')
  expect_identical(plant('matched',function(x) x[-1,]),'schedule-matched.csv')
  # criterion 4 from the audit's own flags: a flagged 300-site fit fails, and a table that says otherwise is caught
  g <- gate_fixture('A1',arm_rows('A1',2,list(v('low',9))),flagged_new='jsdm-n0300-02')
  ag <- audit(g);expect_true(all(ag$checks$pass));expect_identical(unique(ag$gate$result),'FAIL')
  forged <- g$tabs;forged$convergence$selected_flagged[forged$convergence$sd==2] <- 0L
  expect_true('convergence.csv' %in% failing_tables(audit(g,forged)))
  # A2: both gated groups
  a2 <- gate_fixture('A2',arm_rows('A2',2,list(v('low',9),v('rare_1_5pct',7,community=phase_keys('A2')[1:5]))))
  ag <- audit(a2,phase='A2');expect_true(all(ag$checks$pass));expect_identical(unique(ag$gate$result),'FAIL')
})

test_that('verify.R audits the reachable R27 states: an absent SD, and an SD with its manual record', {
  need_archives()
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  audit <- function(f) ve$audit_gate_tables('A1',c(1,2),f$own,f$own_b0,f$flags,NULL,f$tabs)
  # no selected fit at all: criteria 1-3 NA in both implementations, criterion 4 FAIL
  f <- gate_fixture('A1',arm_rows('A1',2)[0,])
  ag <- audit(f);expect_true(all(ag$checks$pass));expect_identical(unique(ag$gate$result),'FAIL')
  expect_true(all(is.na(ag$gate$pass[ag$gate$criterion!='c4'])))
  # one community missing under the manual record
  ok <- arm_rows('A1',2,list(v('low',9)))
  f <- gate_fixture('A1',ok[ok$key!='jsdm-n0100-03',],missing=data.frame(sd=2,key='jsdm-n0100-03'))
  ag <- audit(f);expect_true(all(ag$checks$pass));expect_identical(unique(ag$gate$result),'FAIL')
  # the fits the audit expects, read from the selection and the manual record
  s <- synthetic_study(new_sds=3,supplement_sd=2)
  sel <- ve$audit_selection('A1',c(1,2,3),repo,s$study,s$out)
  expect_identical(sum(sel$selected$sd==2),19L);expect_identical(sum(sel$selected$sd==3),20L)
  expect_identical(sel$missing$key[sel$missing$sd==2],'jsdm-n0100-03');expect_true(any(sel$flag_rows$sd==2))
  s <- synthetic_study(new_sds=3,absent_sd=2)
  sel <- ve$audit_selection('A1',c(1,2,3),repo,s$study,s$out)
  expect_identical(sum(sel$selected$sd==2),0L);expect_identical(sum(sel$missing$sd==2),20L)
})

# ---------------------------------------------------------------------------
# Phase B (two-stage): estimand, R24 pooling, secondary outcomes, gate and the
# final decision (README "Occupancy probability, truth and scored cells",
# "Gate", "Secondary outcomes"; AMENDMENT-2.md R24 and R27; R22)
# ---------------------------------------------------------------------------

test_that('phase B keys: qnear and qfar of one replicate are one generating community at two contamination levels (R24)', {
  k <- phase_keys('B')
  expect_identical(phase_stratum('B',k),rep(c('qnear','qfar'),each=10))
  expect_identical(phase_community('B',k),rep(reps,2))
  # the phase A strata and communities are those of analysis.R
  for(p in c('A1','A2')) {
    expect_identical(phase_stratum(p,phase_keys(p)),stratum_of(p,phase_keys(p)))
    expect_identical(phase_community(p,phase_keys(p)),community_of(p,phase_keys(p)))
  }
  expect_error(b_stratum('design-qnear_K6-q20-01'),'phase B key')
})

# The A1 hand example as a phase B community: 4 fitted sites, of which sites 1
# and 2 are original; sp3 is rare (all-site mean .125). The scored cells are the
# original sites of every species and every fitted site of the rare species, so
# the other cells (sites 3 and 4 of sp1 and sp2) carry no interval.
hand_b <- function(estimate=NULL) {
  a <- hand_a1();cov <- a$covered;cov[3:4,1:2] <- NA
  b0 <- data.frame(species=1:3,truth=c(-1,0,2),estimate=c(-.5,.25,1),lower=c(-2,-1,.5),upper=c(0,1,1.8))
  b0$bias <- b0$estimate-b0$truth;b0$covered <- b0$lower<=b0$truth & b0$truth<=b0$upper
  # collection intercept: species 1 and 3 covered, species 2 not (its truth -.2 lies below the lower bound -.1)
  th <- data.frame(species=1:3,truth=c(.5,-.2,1),estimate=c(.7,.1,.9),lower=c(0,-.1,.5),upper=c(1,.4,1.2))
  th$bias <- th$estimate-th$truth;th$covered <- th$lower<=th$truth & th$truth<=th$upper;th$correlation <- c(-.6,-.2,.2)
  b_make_cells(truth=a$truth,estimate=estimate %||% a$estimate,covered=cov,n_original=2L,b0=b0,theta=th)
}

test_that('phase B signed error, MAE and coverage: bands on the original sites, the rare group on all fitted sites, no all-site scope', {
  need_archives()
  expect_identical(b_scored_cells(hand_a1()$truth,2L),c(1L,2L,5L,6L,9L,10L,11L,12L))
  g <- b_group_table(hand_b(),bsc)
  expect_identical(unique(g$scope),'primary');expect_setequal(g$group,c('all','low','middle','high','rare_below_20pct'))
  low <- row_of(g,'primary','low')
  expect_identical(low$cells,2L);expect_equal(low$signed_error,5,tolerance=1e-12);expect_equal(low$coverage,1)
  mid <- row_of(g,'primary','middle')
  expect_identical(mid$cells,3L);expect_equal(mid$signed_error,100*(-.10+.10+.05)/3,tolerance=1e-12)
  expect_equal(mid$mean_abs_cell_error,100*(.10+.10+.05)/3,tolerance=1e-12);expect_equal(mid$coverage,0)
  high <- row_of(g,'primary','high')
  expect_equal(high$signed_error,-5,tolerance=1e-12);expect_equal(high$abs_signed_error,5,tolerance=1e-12)
  all <- row_of(g,'primary','all')
  expect_identical(all$cells,6L);expect_equal(all$mean_abs_cell_error,100*.40/6,tolerance=1e-12);expect_equal(all$coverage,3/6)
  rare <- row_of(g,'primary','rare_below_20pct')
  expect_identical(rare$cells,4L);expect_equal(rare$signed_error,100*(.05+.05+.05-.03)/4,tolerance=1e-12)
  expect_equal(rare$mean_abs_cell_error,100*.18/4,tolerance=1e-12);expect_equal(rare$coverage,.75)
  # a scored cell must have an interval; an unscored one need not
  a <- hand_a1();cov <- a$covered;cov[4,3] <- NA
  expect_error(b_make_cells(a$truth,a$estimate,cov,2L),'scored cell')
  # no rare species: the group is absent and only the original sites are scored
  a$truth[,3] <- c(.25,.15,.30,.30);cov <- a$covered;cov[3:4,] <- NA
  cells <- b_make_cells(a$truth,a$truth,cov,2L)
  expect_false('rare_below_20pct' %in% b_group_table(cells,bsc)$group);expect_identical(b_scored_cells(a$truth,2L),c(1L,2L,5L,6L,9L,10L))
})

test_that('phase B beta_theta intercept bias and coverage, and the B0 correlation, average over species', {
  s <- b_theta_summary(hand_b()$theta)
  expect_identical(s$species,3L);expect_equal(s$bt_bias,mean(c(.2,.3,-.1)),tolerance=1e-12)
  expect_equal(s$bt_abs_bias,mean(c(.2,.3,.1)),tolerance=1e-12);expect_equal(s$bt_coverage,2/3)
  expect_equal(s$b0_bt_correlation,mean(c(-.6,-.2,.2)),tolerance=1e-12)
  expect_identical(names(s),c('species','bt_bias','bt_abs_bias','bt_coverage','b0_bt_correlation'))
})

# A synthetic two-stage fit with the real fits' array layout: 3 species at 4
# fitted sites (2 original), 2 covariates, 2 site factors, 2 samples per site
# and 2 primers, 6 retained draws in each of 2 chains. The package outputs are
# the posterior means the archived design scorer checks. The collection
# intercept draws of species 1 and 2 are exact linear functions of their B0
# draws, so their correlations are 1 and -1.
b_fixture <- function() {
  set.seed(20260929)
  n <- 4L;n0 <- 2L;S <- 3L;P <- 2L;d <- 2L;M <- 2L;np <- 2L;ni <- 6L;nc <- 2L;species <- paste0('sp',1:S)
  xraw <- cbind(c(-1,0,1,2),c(3,1,2,0));X <- unname(scale(xraw))
  traw <- c(.5,-.5,1,0,-1,2,.3,-.3);Xt <- unname(cbind(1,scale(traw)))
  site <- rep(seq_len(n),each=M*np);sample <- rep(seq_len(n*M),each=np);primer <- rep(seq_len(np),times=n*M)
  info <- data.frame(Site=site,Sample=sample,Primer=primer,X_psi.EnvCov.1=xraw[site,1],X_psi.EnvCov.2=xraw[site,2],X_theta=traw[sample])
  B0 <- qlogis(c(.08,.5,.93));B <- rbind(c(.2,-.3,.1),c(-.1,.2,.05))
  U <- cbind(c(.1,-.1,.2,0),c(0,.1,-.2,.1));L <- rbind(c(.5,.2,-.4),c(-.3,.1,.2))
  jp <- list(B0=B0,B=B,U=U,L=L,eta=sweep(X%*%B+U%*%L,2,B0,'+'))
  bt <- rbind(c(.5,-.2,1),c(1,-1,0))
  tp <- list(jsdmParams_true=jp,beta_theta_true=bt,p_true=matrix(c(.9,.8,.85,.95,.7,.75),np,S),q_true=matrix(.05,np,S),
    w_true=matrix(rbinom(n*M*S,1,.5),n*M,S))
  OTU <- matrix(rpois(nrow(info)*S,3),nrow(info),S,dimnames=list(NULL,species))
  input <- list(scenario=list(S=S,n=n,model='two_stage',ncov_psi=P,P=np),truth=list(params=list(theta0=c(.02,.03,.04))),
    sim=list(true_params=tp,data_list=list(info=info,OTU=OTU)),design=list(original_sites=n0))
  B0d <- array(rnorm(S*ni*nc,B0,.3),c(S,ni,nc));Bd <- array(rnorm(P*S*ni*nc,as.vector(B),.1),c(P,S,ni,nc))
  Ud <- array(rnorm(n*d*ni*nc,as.vector(U),.1),c(n,d,ni,nc));Ld <- array(rnorm(d*S*ni*nc,as.vector(L),.1),c(d,S,ni,nc))
  btd <- array(rnorm(2*S*ni*nc,as.vector(bt),.2),c(2,S,ni,nc));btd[1,1,,] <- 2*B0d[1,,]+1;btd[1,2,,] <- .5-B0d[2,,]
  # every draw's probability, written out cell by cell
  p <- array(NA_real_,c(n,S,ni,nc));th <- array(NA_real_,c(n*M,S,ni,nc))
  for(ch in 1:nc) for(it in 1:ni) for(s in 1:S) {
    for(i in 1:n) p[i,s,it,ch] <- plogis(B0d[s,it,ch]+sum(X[i,]*Bd[,s,it,ch])+sum(Ud[i,,it,ch]*Ld[,s,it,ch]))
    for(j in 1:(n*M)) th[j,s,it,ch] <- plogis(sum(Xt[j,]*btd[,s,it,ch]))
  }
  ro <- list(jsdm_output=list(B0_output=B0d,B_output=Bd,U_output=Ud,L_output=Ld,sigmah_output=matrix(runif(ni*nc,.5,1.5),ni,nc)),
    beta_theta_output=btd,p_output=array(runif(np*S*ni*nc,.6,.95),c(np,S,ni,nc)),q_output=array(runif(np*S*ni*nc,.01,.1),c(np,S,ni,nc)),
    theta0_output=array(runif(S*ni*nc,.01,.05),c(S,ni,nc)),psi_output=apply(p,c(1,2),mean),theta_output=apply(th,c(1,2),mean))
  fit <- list(results_output=ro,infos=list(ps=0,model='two_stage',speciesNames=species),X_psi=X,X_theta=Xt)
  list(fit=fit,input=input,job=list(family='design',arm='sites300',priors=list()),p=p,truth=plogis(jp$eta),
    bt_truth=bt[1,]+mean(traw)*bt[2,],B0d=B0d,btd=btd)
}
pearson <- function(x,y) sum((x-mean(x))*(y-mean(y)))/sqrt(sum((x-mean(x))^2)*sum((y-mean(y))^2))

test_that('phase B cells from a synthetic two-stage fit, through the archived design scorer', {
  need_archives()
  x <- b_fixture()
  cells <- b_cells(x$fit,x$input,x$job,bsc)
  expect_identical(cells$phase,'B');expect_identical(c(cells$n,cells$S,cells$n_original),c(4L,3L,2L))
  expect_equal(cells$truth,x$truth,tolerance=1e-12)
  expect_equal(cells$estimate,apply(x$p,c(1,2),mean),tolerance=1e-12)
  draws <- function(i,s) as.vector(x$p[i,s,,])
  # sp1 (truth mean about .08) is rare: all its sites are scored; sp2 and sp3 on the original sites only
  scored <- rbind(c(1,1),c(2,1),c(3,1),c(4,1),c(1,2),c(2,2),c(1,3),c(2,3))
  for(k in seq_len(nrow(scored))) {i <- scored[k,1];s <- scored[k,2]
    expect_equal(cells$lower[i,s],q7(draws(i,s),.025),tolerance=1e-12);expect_equal(cells$upper[i,s],q7(draws(i,s),.975),tolerance=1e-12)
    expect_identical(cells$covered[i,s],q7(draws(i,s),.025)<=x$truth[i,s] && x$truth[i,s]<=q7(draws(i,s),.975))}
  expect_true(all(is.na(cells$covered[3:4,2:3])))
  expect_lt(cells$checks$estimate_vs_score_fit,1e-12);expect_lt(cells$checks$groups_vs_score_fit,1e-10)
  # B0: the generating B0 on the logit scale, bias and type-7 containment per species
  B0 <- x$input$sim$true_params$jsdmParams_true$B0
  expect_equal(cells$b0$bias,vapply(1:3,function(s) mean(x$B0d[s,,]),numeric(1))-B0,tolerance=1e-12)
  expect_identical(cells$b0$covered,vapply(1:3,function(s) q7(x$B0d[s,,],.025)<=B0[s] && B0[s]<=q7(x$B0d[s,,],.975),logical(1)))
  # beta_theta intercept: truth on the standardised collection-covariate scale (alpha + mean(raw) beta)
  expect_equal(cells$theta$truth,x$bt_truth,tolerance=1e-12)
  expect_equal(cells$theta$bias,vapply(1:3,function(s) mean(x$btd[1,s,,]),numeric(1))-x$bt_truth,tolerance=1e-12)
  expect_identical(cells$theta$covered,vapply(1:3,function(s) {v <- x$btd[1,s,,];q7(v,.025)<=x$bt_truth[s] && x$bt_truth[s]<=q7(v,.975)},logical(1)))
  # B0 and beta_theta intercept correlation over the 12 pooled draws of each species
  expect_equal(cells$theta$correlation,c(1,-1,pearson(as.vector(x$B0d[3,,]),as.vector(x$btd[1,3,,]))),tolerance=1e-12)
  g <- b_group_table(cells,bsc)
  expect_identical(row_of(g,'primary','rare_below_20pct')$cells,4L);expect_identical(row_of(g,'primary','all')$cells,6L)
  e <- as.vector(cells$estimate-x$truth)[c(1,2,5,6,9,10)]
  expect_equal(row_of(g,'primary','all')$signed_error,100*mean(e),tolerance=1e-10)
  expect_equal(b_theta_summary(cells$theta)$b0_bt_correlation,mean(cells$theta$correlation),tolerance=1e-12)
})

test_that('phase B criterion 1 pools qnear and qfar of each replicate (R24); the no-harm criteria hold at each level', {
  # replicates 01-07: qnear 6, qfar 13, mean 9.5 < 10, improved although qfar alone is worse;
  # 08-10: qnear 9, qfar 11, mean exactly 10: a tie, not improved
  new <- arm_rows('B',2,list(v('low',6,'qnear',reps[1:7]),v('low',13,'qfar',reps[1:7]),v('low',9,'qnear',reps[8:10]),v('low',11,'qfar',reps[8:10])))
  g <- gate('B',new)
  c1 <- detail(g,'c1','low')
  expect_identical(nrow(c1),1L);expect_identical(c1$stratum,'pooled');expect_identical(c1$communities,10L)
  expect_identical(c1$needed,7L);expect_identical(c1$improved,7L);expect_equal(c1$value_new,(7*9.5+3*10)/10);expect_true(c1$pass)
  # the gated group is the low band only; the replicate 07 rare group is descriptive (R22)
  expect_identical(unique(detail(g,'c1')$component),'low')
  # criteria 2 to 4 are evaluated at qnear and at qfar separately: here the qfar low band is 2.4 points worse
  expect_setequal(unique(g$detail$stratum[g$detail$criterion %in% c('c2','c3','c4')]),c('qnear','qfar'))
  d <- detail(g,'c2','low','qfar');expect_equal(d$value_new-d$value_control,(7*13+3*11)/10-10);expect_false(d$pass)
  expect_true(detail(g,'c2','low','qnear')$pass)
  # so phase B fails although pooled criterion 1 passes
  expect_true(g$row$c1_pass);expect_false(g$row$c2_pass);expect_identical(g$row$result,'FAIL')
  # improved by one point at both levels: every criterion passes
  ok <- gate('B',arm_rows('B',2,list(v('low',9))))
  expect_identical(ok$row$result,'PASS');expect_identical(detail(ok,'c1','low')$improved,10L)
  # six of ten pooled communities improved fails
  six <- gate('B',arm_rows('B',2,list(v('low',9),v('low',11,community=reps[7:10]))))
  expect_identical(detail(six,'c1','low')$improved,6L);expect_false(detail(six,'c1','low')$pass)
  # however the replicate 07 rare group moves, it is not gated
  expect_identical(gate('B',arm_rows('B',2,list(v('low',9),v('rare_below_20pct',90))))$row$result,'PASS')
})

test_that('phase B criterion 4 compares flagged selected fits with the control at each contamination level', {
  ok <- arm_rows('B',2,list(v('low',9)))
  # the pr11 control selection flags 1 qnear and 3 qfar fits
  ctl <- transform(conv_rows('B',1),selected_flagged=c(1L,3L))
  pass <- gate('B',ok,conv=rbind(ctl,transform(conv_rows('B',2),selected_flagged=c(1L,3L))))
  expect_true(all(detail(pass,'c4')$pass));expect_identical(pass$row$result,'PASS')
  expect_identical(detail(pass,'c4',stratum='qfar')$value_control,3)
  near <- gate('B',ok,conv=rbind(ctl,transform(conv_rows('B',2),selected_flagged=c(2L,0L))))
  expect_false(detail(near,'c4',stratum='qnear')$pass);expect_true(detail(near,'c4',stratum='qfar')$pass)
  expect_identical(near$row$result,'FAIL')
  far <- gate('B',ok,conv=rbind(ctl,transform(conv_rows('B',2),selected_flagged=c(0L,4L))))
  expect_false(detail(far,'c4',stratum='qfar')$pass);expect_identical(far$row$result,'FAIL')
})

test_that('a phase B replicate valid at only one contamination level leaves pooled criterion 1 and fails criterion 4 at that level (R24, R27)', {
  ok <- arm_rows('B',2,list(v('low',9)))
  gone <- ok[ok$key!='design-qfar_K6-sites300-04',]
  conv <- rbind(conv_rows('B',1),transform(conv_rows('B',2),selected_fits=c(10L,9L)))
  g <- gate('B',gone,conv=conv)
  c1 <- detail(g,'c1','low');expect_identical(c1$communities,9L);expect_identical(c1$needed,6L);expect_true(c1$pass)
  expect_false(detail(g,'c4',stratum='qfar')$pass);expect_identical(detail(g,'c4',stratum='qfar')$missing_new,1L)
  expect_true(detail(g,'c4',stratum='qnear')$pass)
  expect_identical(detail(g,'c2','low','qfar')$communities,9L);expect_identical(detail(g,'c2','low','qnear')$communities,10L)
  expect_identical(detail(g,'c3','high','qnear')$communities,10L);expect_identical(g$row$result,'FAIL')
  # the same through the manual record of a community without a valid selected fit
  g <- gate('B',ok[ok$key!='design-qnear_K6-sites300-09',],missing=data.frame(sd=2,key='design-qnear_K6-sites300-09',reason='quarantined'))
  expect_false(detail(g,'c4',stratum='qnear')$pass);expect_identical(detail(g,'c1','low')$communities,9L);expect_identical(g$row$result,'FAIL')
})

# A phase B selection as select.R writes it: selected-fits.csv, and
# convergence.csv with the single stratum 'all' that its size_stratum gives
# every phase B key; the controls flagged as in the pr11 selection.
B_CONTROL_FLAGGED <- sprintf('design-%s_K6-sites300-%s',c('qnear','qfar','qfar','qfar'),c('02','05','07','09'))
b_selection <- function(dir,new_flagged=character(),repeated=character()) {
  keys <- phase_keys('B')
  sf <- rbind(data.frame(phase='B',sd=2,key=keys,role='new',first_schedule='initial',selected_schedule=ifelse(keys %in% repeated,'long','initial'),
      long_repeat=keys %in% repeated,fit=sprintf('fits/B/sd2/%s/%s-fit.rds',ifelse(keys %in% repeated,'long','initial'),keys),
      fit_md5=sprintf('%032d',seq_along(keys)),rule='pr11',first_flagged=keys %in% repeated,flagged=keys %in% new_flagged,reasons='',stringsAsFactors=FALSE),
    data.frame(phase='B',sd=1,key=keys,role='control',first_schedule=NA,selected_schedule='initial',long_repeat=NA,
      fit=sprintf('pr11-current-20260927/initial/%s-fit.rds',keys),fit_md5=sprintf('%032d',100+seq_along(keys)),rule='pr11',
      first_flagged=NA,flagged=keys %in% B_CONTROL_FLAGGED,reasons='',stringsAsFactors=FALSE))
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(sf,file.path(dir,'selected-fits.csv'),row.names=FALSE)
  utils::write.csv(convergence_counts(sf),file.path(dir,'convergence.csv'),row.names=FALSE)
  cf <- sf[sf$sd==1,];cf <- data.frame(role='control',phase='B',sd=1,key=cf$key,schedule='initial',fit=cf$fit,fit_md5=cf$fit_md5,rule='pr11',
    warnings=0L,max_group_rhat=1,max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA,spatial_field_flags=NA,flagged=cf$flagged,reasons='',
    stringsAsFactors=FALSE)
  utils::write.csv(cf,file.path(dir,'control-flags.csv'),row.names=FALSE)
  sf
}

test_that('phase B convergence counts are per contamination level, checked against select.R convergence.csv (R24)', {
  need_archives()
  d <- tempfile('selB')
  b_selection(d,new_flagged=c('design-qnear_K6-sites300-03','design-qfar_K6-sites300-06'),repeated=c('design-qnear_K6-sites300-03','design-qnear_K6-sites300-04'))
  expect_identical(unique(utils::read.csv(file.path(d,'convergence.csv'))$stratum),'all')
  cv <- convergence_table('B',c(1,2),d,repo,read_manual_missing(tempfile(),'B'))
  row <- function(sd,h) cv[cv$sd==sd & cv$stratum==h,,drop=FALSE]
  expect_identical(nrow(cv),4L);expect_setequal(cv$stratum,c('qnear','qfar'))
  expect_identical(row(1,'qnear')$selected_flagged,1L);expect_identical(row(1,'qfar')$selected_flagged,3L)
  expect_identical(row(2,'qnear')$selected_flagged,1L);expect_identical(row(2,'qfar')$selected_flagged,1L)
  expect_identical(row(2,'qnear')$long_repeats,2L);expect_identical(row(2,'qnear')$first_flagged,2L);expect_identical(row(2,'qfar')$long_repeats,0L)
  expect_identical(row(2,'qfar')$selected_fits,10L);expect_identical(row(2,'qfar')$expected_fits,10L)
  # a select.R total that the per-level counts do not add up to is refused
  x <- utils::read.csv(file.path(d,'convergence.csv'));x$selected_flagged[x$sd==2] <- 5L;utils::write.csv(x,file.path(d,'convergence.csv'),row.names=FALSE)
  expect_error(convergence_table('B',c(1,2),d,repo,read_manual_missing(tempfile(),'B')),'convergence.csv')
  # before the phase B selection exists, the controls' counts per level come from the archived flag cross-check
  cv <- convergence_table('B',1,tempfile('none'),repo,read_manual_missing(tempfile(),'B'))
  expect_identical(cv$selected_flagged[cv$stratum=='qnear'],1L);expect_identical(cv$selected_flagged[cv$stratum=='qfar'],3L)
})

test_that('the final decision: the smallest SD passing phases A and B is recommended; PASS then FAIL is binary-only; UNDETERMINED is not adopted', {
  a <- function(r) data.frame(sd=c(2,3,5),subphase='A',result=r,stringsAsFactors=FALSE)
  b <- function(r,sds=c(2,3,5)) data.frame(sd=sds,subphase='B',result=r,stringsAsFactors=FALSE)
  d <- final_decision(a(c('PASS','PASS','PASS')),b(c('FAIL','PASS','PASS')))
  expect_identical(names(d),c('sd','phase_a_result','phase_b_result','passes_both','recommended','consequence','study_decision'))
  expect_identical(d$sd,c(2,3,5));expect_identical(d$phase_b_result,c('FAIL','PASS','PASS'))
  expect_identical(d$passes_both,c(FALSE,TRUE,TRUE));expect_identical(d$recommended,c(FALSE,TRUE,FALSE))
  expect_match(d$consequence[1],'binary-only improvement');expect_match(d$consequence[1],'does not change the default')
  expect_match(d$consequence[2],'passes phases A and B');expect_match(d$study_decision[1],'sigma_b0 = 3')
  expect_match(d$study_decision[1],'separate reviewed pull request');expect_identical(length(unique(d$study_decision)),1L)
  # every SD fails phase B or is undetermined: no default change
  d <- final_decision(a(c('PASS','PASS','PASS')),b(c('UNDETERMINED','FAIL','FAIL')))
  expect_false(any(d$recommended));expect_match(d$consequence[1],'undetermined');expect_match(d$consequence[1],'does not change the default')
  expect_match(d$consequence[2],'binary-only');expect_match(d$study_decision[1],'default stays at sigma_b0 = 1')
  # phase B runs only for SDs that pass phase A
  d <- final_decision(a(c('FAIL','PASS','PASS')),b(c('PASS','PASS'),sds=c(3,5)))
  expect_identical(d$phase_b_result,c('NOT RUN','PASS','PASS'));expect_match(d$consequence[1],'fails phase A')
  expect_identical(d$recommended,c(FALSE,TRUE,FALSE))
  expect_error(final_decision(a(c('PASS','PASS','PASS')),b(c('PASS','PASS'),sds=c(3,5))),'phase B result')
  expect_error(final_decision(a(c('FAIL','PASS','PASS')),b(c('PASS','PASS','PASS'))),'did not pass phase A')
})

# A temporary study with synthetic phase B score records built by the phase B
# record constructor: the 20 controls, and 20 records of SD 2 whose cell errors
# are half the control's. The synthetic selection names the real control fits.
synthetic_study_b <- function(new_flagged=character()) {
  st <- tempfile('study');dir.create(st);out <- file.path(st,'summary')
  items <- control_items('B',st,archives,repo);hashes <- b_score_hashes(repo,archives)
  write_record <- function(item,cells) {
    rec <- b_score_record(item,cells,bsc,hashes,input_md5='synthetic');p <- score_path(out,'B',item$sd,item$key,item$schedule)
    dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);saveRDS(rec,p)
  }
  for(i in seq_len(nrow(items))) write_record(items[i,,drop=FALSE],hand_b())
  sel <- file.path(st,'selection/B');sf <- b_selection(sel,new_flagged=new_flagged)
  i <- match(sf$key[sf$sd==1],items$key)
  sf$fit[sf$sd==1] <- items$fit_label[i];sf$fit_md5[sf$sd==1] <- items$expected_md5[i];sf$selected_schedule[sf$sd==1] <- items$schedule[i]
  utils::write.csv(sf,file.path(sel,'selected-fits.csv'),row.names=FALSE)
  utils::write.csv(convergence_counts(sf),file.path(sel,'convergence.csv'),row.names=FALSE)
  cf <- utils::read.csv(file.path(sel,'control-flags.csv'),stringsAsFactors=FALSE);j <- match(cf$key,items$key)
  cf$fit <- items$fit_label[j];cf$fit_md5 <- items$expected_md5[j];cf$schedule <- items$schedule[j]
  utils::write.csv(cf,file.path(sel,'control-flags.csv'),row.names=FALSE)
  better <- hand_b(estimate=hand_a1()$truth+(hand_a1()$estimate-hand_a1()$truth)/2);new <- sf[sf$sd==2,]
  for(k in seq_len(nrow(new))) write_record(data.frame(role='new',phase='B',sd=2,key=new$key[k],schedule='initial',fit=file.path(st,new$fit[k]),
    fit_label=new$fit[k],expected_md5=new$fit_md5[k],kind='selected',stringsAsFactors=FALSE),better)
  list(study=st,out=out)
}

test_that('summarise.R scores phase B end to end: per-level tables, beta_theta tables, the gate in summary/B, and the final decision', {
  need_archives()
  s <- synthetic_study_b()
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=B','--arms=1,2',paste0('--out=',s$out)))
  expect_identical(r$status,0L);expect_false(grepl('PASS|FAIL',r$output))
  d <- file.path(s$out,'B')
  for(f in c('band-error.csv','coverage.csv','mae.csv','b0.csv','beta-theta.csv','convergence.csv','summary-means.csv','b0-means.csv',
    'beta-theta-means.csv','gate.csv','gate-detail.csv','schedule-matched.csv','provenance.csv')) expect_true(file.exists(file.path(d,f)),info=f)
  g <- utils::read.csv(file.path(d,'gate.csv'),stringsAsFactors=FALSE)
  expect_identical(g$subphase,'B');expect_identical(g$result,'PASS');expect_identical(g$c1_low_communities,10L)
  expect_identical(c(g$c4_qnear_flagged_control,g$c4_qfar_flagged_control),c(1L,3L))
  expect_true(all(c('c2_low_qnear_diff','c2_mae_qfar_diff','c3_high_qfar_drop') %in% names(g)))
  cv <- utils::read.csv(file.path(d,'convergence.csv'),stringsAsFactors=FALSE)
  expect_setequal(cv$stratum,c('qnear','qfar'));expect_identical(sum(cv$selected_flagged[cv$sd==1]),4L)
  bt <- utils::read.csv(file.path(d,'beta-theta.csv'),stringsAsFactors=FALSE)
  expect_identical(nrow(bt),40L);expect_equal(unique(bt$bt_bias),mean(c(.2,.3,-.1)),tolerance=1e-12)
  bm <- utils::read.csv(file.path(d,'beta-theta-means.csv'),stringsAsFactors=FALSE)
  expect_identical(nrow(bm),4L);expect_equal(unique(bm$mean_b0_bt_correlation),mean(c(-.6,-.2,.2)),tolerance=1e-12)
  be <- utils::read.csv(file.path(d,'band-error.csv'),stringsAsFactors=FALSE)
  expect_identical(unique(be$scope),'primary');expect_setequal(unique(be$stratum),c('qnear','qfar'))
  # two flagged selected qnear fits, one more than the control's one, fail criterion 4 at qnear
  s2 <- synthetic_study_b(new_flagged=c('design-qnear_K6-sites300-01','design-qnear_K6-sites300-05'))
  expect_identical(run_cli('summarise.R',cli_args(s2$study,'--phase=B','--arms=1,2',paste0('--out=',s2$out)))$status,0L)
  g2 <- utils::read.csv(file.path(s2$out,'B/gate.csv'),stringsAsFactors=FALSE)
  expect_identical(g2$result,'FAIL');expect_false(g2$c4_qnear_pass);expect_true(g2$c4_qfar_pass)
  # the final decision reads the phase A and phase B gates
  dir.create(file.path(s$out,'A'));utils::write.csv(data.frame(sd=2,subphase='A',result='PASS'),file.path(s$out,'A/gate.csv'),row.names=FALSE)
  r <- run_cli('summarise.R',c(paste0('--repo=',repo),paste0('--study=',s$study),'--phase=final',paste0('--out=',s$out)))
  expect_identical(r$status,0L);expect_false(grepl('PASS|FAIL',r$output))
  fd <- utils::read.csv(file.path(s$out,'final/decision.csv'),stringsAsFactors=FALSE)
  expect_identical(fd$sd,2L);expect_true(fd$recommended);expect_match(fd$study_decision,'sigma_b0 = 2')
})

test_that('analysis-b.R refuses new arms before the phase B selection exists, invalid arms and other phases', {
  need_archives()
  empty <- tempfile('study');dir.create(empty)
  r <- run_cli('analysis-b.R',cli_args(empty,'--phase=B','--arms=1,2',paste0('--out=',tempfile())))
  expect_identical(r$status,1L);expect_match(r$output,'selection/B/selected-fits.csv')
  r <- run_cli('analysis-b.R',cli_args(empty,'--phase=B','--arms=4'))
  expect_identical(r$status,1L);expect_match(r$output,'--arms')
  r <- run_cli('analysis-b.R',cli_args(empty,'--phase=A1','--arms=1'))
  expect_identical(r$status,1L);expect_match(r$output,'phase B')
  expect_identical(list.files(empty,recursive=TRUE),character())
})

test_that('an archived phase B control is scored exactly as its archived result (design-qnear_K6-sites300-01)', {
  need_archives()
  items <- control_items('B',study,archives,repo);item <- items[items$key=='design-qnear_K6-sites300-01',,drop=FALSE]
  expect_identical(item$schedule,'initial')
  out <- tempfile('scores')
  expect_no_warning(path <- b_score_item(item,bsc,archives,inputs_root,out))
  r <- readRDS(path);old <- readRDS(file.path(archives,'pr11-current-20260927/initial/design-qnear_K6-sites300-01-result.rds'))
  expect_identical(r$version,B_SCORE_VERSION);expect_identical(r$role,'control');expect_identical(r$sd,1)
  expect_identical(r$stratum,'qnear');expect_identical(r$community,'01')
  expect_identical(r$cells$estimate,as.vector(old$estimate));expect_identical(r$cells$truth,as.vector(old$truth))
  expect_identical(sum(!is.na(r$cells$covered)),1000L)
  for(b in c('all','low','middle','high')) {
    mine <- row_of(r$groups,'primary',b);theirs <- old$groups[old$groups$metric=='occupancy_original_sites' & old$groups$group==if(b=='middle') 'medium' else b,]
    expect_identical(mine$cells,as.integer(theirs$n_elements))
    expect_lt(abs(mine$signed_error-100*theirs$bias),1e-10);expect_lt(abs(mine$mean_abs_cell_error-100*theirs$mae),1e-10)
  }
  el <- function(m) old$elements[old$elements$metric==m,]
  expect_lt(max(abs(r$b0_species$bias-el('B0')$bias)),1e-12);expect_lt(max(abs(r$theta_species$bias-el('collection_intercept')$bias)),1e-12)
  expect_identical(r$theta_species$truth,el('collection_intercept')$truth)
  expect_true(all(abs(r$theta_species$correlation)<=1))
  mtime <- file.mtime(path);expect_identical(b_score_item(item,bsc,archives,inputs_root,out),path);expect_identical(file.mtime(path),mtime)
})

test_that('the phase B new-arm loading path checks the fit it scores (a synthetic new-arm copy of a control fit; no new-arm fit is read)', {
  need_archives()
  key <- 'design-qfar_K6-sites300-10';items <- control_items('B',study,archives,repo);ctl <- items[items$key==key,,drop=FALSE]
  expect_identical(ctl$schedule,'initial')
  st <- tempfile('study');saved <- readRDS(ctl$fit);spec <- phase_jobs('B',archives,inputs_root,key)[[1]]
  saved$fit$infos$intercept_prior <- list(mean=0,sd=2)
  dest <- fit_path(st,'B',2,'initial',key);dir.create(dirname(dest),recursive=TRUE)
  saveRDS(list(fit=saved$fit,warnings=saved$warnings,phase='B',key=key,sd=2,schedule='initial',mcmc=schedule_mcmc('B','initial',archives),
    input_md5=spec$input_md5),dest,compress=FALSE)
  item <- data.frame(role='new',phase='B',sd=2,key=key,schedule='initial',fit=dest,fit_label=path_within(dest,st),
    expected_md5=unname(tools::md5sum(dest)),kind='selected',stringsAsFactors=FALSE)
  out <- tempfile('scores')
  expect_no_warning(p <- b_score_item(item,bsc,archives,inputs_root,out));r <- readRDS(p)
  expect_identical(r$role,'new');expect_identical(r$sd,2);expect_identical(nrow(r$cells),3000L);expect_identical(r$stratum,'qfar')
  expect_true(all(is.finite(r$groups$signed_error)));expect_identical(r$theta$species,10L)
  bad <- item;bad$sd <- 3;expect_error(b_score_item(bad,bsc,archives,inputs_root,tempfile('scores')),'not the expected fit')
  bad <- item;bad$expected_md5 <- strrep('0',32);expect_error(b_score_item(bad,bsc,archives,inputs_root,tempfile('scores')),'md5')
  unlink(st,recursive=TRUE)
})

test_that('verify.R recomputes the phase B cells, groups, B0, beta_theta intercept and correlation independently', {
  need_archives()
  src <- readLines(file.path(here,'verify.R'))
  expect_false(any(grepl('source',src,fixed=TRUE) & grepl('analysis',src,fixed=TRUE)))
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  x <- b_fixture()
  mine <- ve$audit_b_cells(x$fit,x$input)
  cx <- b_cells(x$fit,x$input,x$job,bsc);theirs <- b_group_table(cx,bsc)
  own <- ve$audit_groups(mine);m <- merge(own,theirs,by=c('scope','group'))
  expect_identical(nrow(m),nrow(theirs));expect_identical(nrow(own),nrow(theirs))
  expect_lt(max(abs(m$signed_error.x-m$signed_error.y)),1e-10);expect_lt(max(abs(m$mean_abs_cell_error.x-m$mean_abs_cell_error.y)),1e-10)
  expect_lt(max(abs(m$coverage.x-m$coverage.y)),1e-12);expect_identical(m$cells.x,m$cells.y)
  b0 <- ve$audit_b0(x$fit$results_output$jsdm_output$B0_output,x$input$sim$true_params$jsdmParams_true$B0);s0 <- b0_summary(cx$b0)
  expect_lt(max(abs(b0[c('b0_bias','b0_abs_bias','b0_coverage')]-unlist(s0[c('b0_bias','b0_abs_bias','b0_coverage')]))),1e-12)
  th <- ve$audit_theta(x$fit,x$input);st <- b_theta_summary(cx$theta)
  expect_lt(max(abs(th[c('bt_bias','bt_abs_bias','bt_coverage','b0_bt_correlation')]-unlist(st[c('bt_bias','bt_abs_bias','bt_coverage','b0_bt_correlation')]))),1e-12)
  expect_identical(unname(th[['species']]),3)
})

test_that("verify.R's phase B flag diagnostics reproduce flags.R's frozen pr11 rule", {
  need_archives()
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  x <- b_fixture();fl_sc <- load_flag_scorers(repo,archives)
  mine <- ve$audit_b_cells(x$fit,x$input)
  for(w in list(character(),'a fitting warning')) {
    fl <- fit_flags('B',x$fit,w,x$input,fl_sc);d <- ve$audit_diag_b(x$fit,x$input,mine,length(w))$summary
    expect_lt(abs(d$max_group_rhat-fl$max_group_rhat),1e-12);expect_lt(abs(d$max_element_rhat-fl$max_element_rhat),1e-12)
    expect_identical(as.integer(d$unresolved_rhat),as.integer(fl$unresolved_rhat));expect_identical(d$flagged,fl$flagged)
  }
  expect_true(ve$audit_diag_b(x$fit,x$input,mine,1L)$summary$flagged)
})

test_that('verify.R recomputes the phase B tables and gate (R24) and catches planted errors', {
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  new <- arm_rows('B',2,list(v('low',6,'qnear',reps[1:7]),v('low',9.5,'qfar',reps[1:7]),v('low',9,'qnear',reps[8:10]),v('low',11,'qfar',reps[8:10])))
  f <- gate_fixture('B',new,flagged_new='design-qfar_K6-sites300-03')
  th <- f$own_b0[c('phase','sd','arm','kind','key','stratum','community')];th$species <- 10L
  th$bt_bias <- ifelse(th$sd==1,.2,.1);th$bt_abs_bias <- .3;th$bt_coverage <- .9;th$b0_bt_correlation <- -.4
  tabs <- f$tabs;tabs$btmeans <- beta_theta_means(th)
  audit <- function(t=tabs) ve$audit_gate_tables('B',c(1,2),f$own,f$own_b0,f$flags,NULL,t,own_theta=th)
  ag <- audit()
  expect_true(all(ag$checks$pass));expect_identical(unique(ag$gate$result),f$gate$row$result)
  expect_true('beta-theta-means.csv' %in% ag$checks$table)
  expect_identical(ag$gate$stratum[ag$gate$criterion=='c1'],'pooled');expect_setequal(ag$gate$stratum[ag$gate$criterion=='c4'],c('qnear','qfar'))
  # no control fit is flagged in this fixture, so one flagged qfar fit fails criterion 4 at qfar only
  expect_false(ag$gate$pass[ag$gate$criterion=='c4' & ag$gate$stratum=='qfar']);expect_true(ag$gate$pass[ag$gate$criterion=='c4' & ag$gate$stratum=='qnear'])
  plant <- function(part,fn) {t <- tabs;t[[part]] <- fn(t[[part]]);failing_tables(audit(t))}
  expect_identical(plant('btmeans',function(x) {x$mean_b0_bt_correlation[1] <- 0;x}),'beta-theta-means.csv')
  expect_identical(plant('gate',function(x) transform(x,c1_low_improved=c1_low_improved-1L)),'gate.csv')
  expect_identical(plant('detail',function(x) {i <- x$criterion=='c2' & x$stratum=='qfar' & x$component=='low';x$pass[i] <- !x$pass[i];x}),'gate-detail.csv')
  expect_identical(plant('convergence',function(x) {x$selected_flagged[x$sd==2 & x$stratum=='qfar'] <- 0L;x}),'convergence.csv')
})

test_that('verify.R audits the final decision table from its own phase A and phase B results', {
  ve <- new.env();sys.source(file.path(here,'verify.R'),envir=ve)
  a <- data.frame(sd=c(2,3,5),result=c('PASS','PASS','PASS'),stringsAsFactors=FALSE)
  cols <- c('phase_a_result','phase_b_result','passes_both','recommended','consequence','study_decision')
  for(rb in list(c('FAIL','PASS','PASS'),c('UNDETERMINED','FAIL','FAIL'),c('PASS','PASS','PASS'))) {
    b <- data.frame(sd=c(2,3,5),result=rb,stringsAsFactors=FALSE)
    mine <- ve$audit_decision(a,b);theirs <- final_decision(transform(a,subphase='A'),transform(b,subphase='B'))
    expect_true(all(ve$compare_rows('final/decision.csv',mine,theirs,'sd',cols)$pass))
  }
  # every SD passes both phases: the smallest, SD 2, is recommended, so marking SD 3 is an error
  theirs$recommended <- c(FALSE,TRUE,FALSE)
  expect_false(all(ve$compare_rows('final/decision.csv',ve$audit_decision(a,b),theirs,'sd',cols)$pass))
  # phase A failed for SD 2: phase B was not run for it
  a2 <- transform(a,result=c('FAIL','PASS','PASS'));b2 <- data.frame(sd=c(3,5),result=c('FAIL','PASS'),stringsAsFactors=FALSE)
  mine <- ve$audit_decision(a2,b2);expect_identical(mine$phase_b_result,c('NOT RUN','FAIL','PASS'));expect_identical(mine$recommended,c(FALSE,FALSE,TRUE))
})
