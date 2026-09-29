#!/usr/bin/env Rscript
# Research tests for the phase A scorer: analysis.R (per-community outcomes),
# summarise.R (tables and the gate) and verify.R (independent audit). Run from
# this directory with `Rscript test-analysis.R`; the exit status is 1 if any
# expectation fails or any block errors.
#
# Blind development (ruling R28): every metric and gate criterion is checked on
# hand-computed synthetic examples. Real fits enter only as mechanics fixtures:
# one archived control fit (whose archived result the scorer must reproduce)
# and the sd 3 pilot fit (a pilot-schedule fit that is never scored; only its
# loading path is exercised and no outcome of it is printed or compared).
# No new-arm fit is read. Tests that need the archives skip when absent.
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
for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R','summarise.R')) source(f)
here <- normalizePath('.')
repo <- normalizePath('../../..')
archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
inputs_root <- file.path(archives,'intercept-prior-inputs')
study <- file.path(archives,'intercept-prior-20260929')
have_archives <- all(dir.exists(file.path(archives,c(unname(ARCHIVES),'intercept-prior-inputs'))))
need_archives <- function() if(!have_archives) skip('saved study archives not available')
rscript <- file.path(R.home('bin'),'Rscript')
run_cli <- function(script,args) {
  out <- suppressWarnings(system2(rscript,c(shQuote(file.path(here,script)),shQuote(args)),stdout=TRUE,stderr=TRUE))
  list(status=attr(out,'status') %||% 0L,output=paste(out,collapse='\n'))
}
# The archived scorer definitions, md5-checked and loaded on first use.
delayedAssign('sc',if(have_archives) load_metric_scorers(repo,archives) else NULL)
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
  cells <- suppressWarnings(a1_cells(x$fit,x$input,sc))
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
  y <- a1_fixture(120L);cy <- suppressWarnings(a1_cells(y$fit,y$input,sc));gy <- group_table(cy,sc)
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
  cells <- suppressWarnings(a2_cells(x$fit,x$input,sc,knots=2L))
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
  keys <- if(phase=='A1') phase_keys('A1') else phase_keys('A2')
  base <- list(low=10,middle=5,high=5,rare_1_5pct=8,rare_below_20pct=20,all=0)
  cov <- list(low=.8,middle=.8,high=.8,rare_1_5pct=.8,rare_below_20pct=.8,all=.8)
  groups <- if(phase=='A1') c('low','middle','high','all','rare_below_20pct') else c('low','middle','high','all','rare_1_5pct')
  rows <- list()
  for(k in keys) for(g in groups) {
    if(g=='rare_below_20pct' && !grepl('-07$',k)) next
    rows[[length(rows)+1L]] <- data.frame(phase=phase,sd=sd,key=k,stratum=stratum_of(phase,k),community=community_of(phase,k),
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
  strata <- if(phase=='A1') c('n100','n300') else 'all';n <- if(phase=='A1') 10L else 9L
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
  expect_true(g$row$c2_pass)
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
  g <- gate('A1',arm_rows('A1',2,c(ok,list(v('middle',.4999,column='coverage')))),control=base)
  expect_false(detail(g,'c3','middle','n100')$pass);expect_false(g$row$c3_pass);expect_identical(g$row$result,'FAIL')
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
  # a manual record of a missing fit (AMENDMENT-2) has the same effect even if a score exists
  g <- gate('A1',ok,missing=data.frame(sd=2,key='jsdm-n0100-09',reason='quarantined'))
  expect_false(detail(g,'c4',stratum='n100')$pass);expect_identical(detail(g,'c1','low')$communities,9L)
  expect_identical(g$row$result,'FAIL')
  # the control arm must be complete
  expect_error(gate('A1',ok,control=arm_rows('A1',1)[arm_rows('A1',1)$key!='jsdm-n0100-01',]),'control')
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
  r <- readRDS(suppressWarnings(score_item(item,sc,archives,inputs_root,out)))
  # Structure only: no outcome of the pilot fit is compared or printed.
  expect_identical(nrow(r$cells),1000L);expect_true(all(c('primary') %in% r$groups$scope))
  expect_true(all(is.finite(r$groups$signed_error)));expect_identical(r$schedule,'pilot')
  bad <- item;bad$sd <- 2
  expect_error(score_item(bad,sc,archives,inputs_root,tempfile('scores')),'not the expected fit')
  bad <- item;bad$expected_md5 <- strrep('0',32)
  expect_error(score_item(bad,sc,archives,inputs_root,tempfile('scores')),'md5')
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
# own record constructor, for summarise.R end to end.
synthetic_study <- function(new_sd=NULL,missing=NULL,flag_n300=0L) {
  st <- tempfile('study');dir.create(st);out <- file.path(st,'summary')
  items <- control_items('A1',st,archives,repo)
  hashes <- current_score_hashes(repo,archives)
  write_record <- function(item,cells) {
    rec <- score_record(item,cells,sc,hashes,input_md5='synthetic');p <- score_path(out,'A1',item$sd,item$key,item$schedule)
    dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);saveRDS(rec,p)
  }
  base <- hand_a1()
  for(i in seq_len(nrow(items))) write_record(items[i,,drop=FALSE],base)
  if(!is.null(new_sd)) {
    better <- base;better$estimate <- better$truth+(base$estimate-base$truth)/2
    sel <- items;sel$role <- 'new';sel$sd <- new_sd;sel$schedule <- 'initial'
    sel$fit_label <- sprintf('fits/A1/sd%d/initial/%s-fit.rds',new_sd,sel$key);sel$fit <- file.path(st,sel$fit_label)
    sel$expected_md5 <- sprintf('%032d',seq_len(nrow(sel)))
    keep <- !sel$key %in% (missing %||% character())
    for(i in which(keep)) write_record(sel[i,,drop=FALSE],better)
    d <- file.path(st,'selection/A1');dir.create(d,recursive=TRUE)
    cf <- data.frame(role='control',phase='A1',sd=1,key=items$key,schedule=items$schedule,fit=items$fit_label,fit_md5=items$expected_md5,
      rule='pr11',warnings=0L,max_group_rhat=1,max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA,spatial_field_flags=NA,
      flagged=FALSE,reasons='',stringsAsFactors=FALSE)
    utils::write.csv(cf,file.path(d,'control-flags.csv'),row.names=FALSE)
    flag <- sel$key %in% sprintf('jsdm-n0300-%02d',seq_len(flag_n300))
    sf <- rbind(data.frame(phase='A1',sd=new_sd,key=sel$key,role='new',first_schedule='initial',selected_schedule='initial',long_repeat=FALSE,
        fit=sel$fit_label,fit_md5=sel$expected_md5,rule='pr11',first_flagged=flag,flagged=flag,reasons=ifelse(flag,'stub',''),stringsAsFactors=FALSE)[keep,],
      data.frame(phase='A1',sd=1,key=items$key,role='control',first_schedule=NA,selected_schedule=items$schedule,long_repeat=NA,
        fit=items$fit_label,fit_md5=items$expected_md5,rule='pr11',first_flagged=NA,flagged=FALSE,reasons='',stringsAsFactors=FALSE))
    utils::write.csv(sf,file.path(d,'selected-fits.csv'),row.names=FALSE)
    utils::write.csv(convergence_counts(sf),file.path(d,'convergence.csv'),row.names=FALSE)
    if(length(missing)) utils::write.csv(data.frame(phase='A1',sd=new_sd,key=missing,reason='quarantined, see log'),
      file.path(d,'manual-missing.csv'),row.names=FALSE)
  }
  list(study=st,out=out)
}

test_that('summarise.R with --arms=1 writes the control tables and no gate', {
  need_archives()
  s <- synthetic_study()
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1',paste0('--out=',s$out)))
  expect_identical(r$status,0L)
  d <- file.path(s$out,'A1')
  for(f in c('band-error.csv','coverage.csv','mae.csv','b0.csv','convergence.csv','summary-means.csv','provenance.csv'))
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

test_that('summarise.R gates a new arm end to end, including R20 and R27', {
  need_archives()
  s <- synthetic_study(new_sd=2)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))
  expect_identical(r$status,0L)
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE)
  # the new arm halves every cell error: all bands better, coverage unchanged, no flags
  expect_identical(nrow(g),1L);expect_identical(g$sd,2L);expect_identical(g$result,'PASS')
  expect_identical(g$c1_low_improved,10L);expect_identical(g$c1_low_needed,7L)
  expect_true(file.exists(file.path(s$out,'A1/gate-detail.csv')))
  # one flagged selected fit at 300 sites fails the gate (R20)
  s <- synthetic_study(new_sd=2,flag_n300=1L)
  expect_identical(run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))$status,0L)
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE)
  expect_identical(g$result,'FAIL');expect_identical(g$c4_n300_flagged_new,1L);expect_true(g$c1_pass)
  # a community recorded by hand as lacking a valid fit fails criterion 4 (R27)
  s <- synthetic_study(new_sd=2,missing='jsdm-n0100-03')
  expect_identical(run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))$status,0L)
  g <- utils::read.csv(file.path(s$out,'A1/gate.csv'),stringsAsFactors=FALSE)
  expect_identical(g$result,'FAIL');expect_identical(g$c4_n100_missing_new,1L);expect_identical(g$c1_low_communities,9L)
  # a score record that does not match the selected fit is refused
  s <- synthetic_study(new_sd=2)
  p <- score_path(s$out,'A1',2,'jsdm-n0100-05','initial');x <- readRDS(p);x$fit_md5 <- strrep('f',32);saveRDS(x,p)
  r <- run_cli('summarise.R',cli_args(s$study,'--phase=A1','--arms=1,2',paste0('--out=',s$out)))
  expect_identical(r$status,1L);expect_match(r$output,'jsdm-n0100-05')
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
  theirs <- group_table(suppressWarnings(a1_cells(x$fit,x$input,sc)),sc)
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
