# Standalone research tests for verdicts.R (Task 3c of PLAN.md in this
# directory). Run with `Rscript test-verdicts.R` from this directory.
# Hand-built tables only: no fit or archive is read.
library(testthat)
source('modes.R')
source('verdicts.R')

# anatomy_regions() rows from hand-built labels: each chain's labels are all
# near-truth, all mirror, or alternate ('both'); `far` is each chain's share of
# draws far from both anchored components.
regions_of <- function(regions,far=rep(.005,length(regions)),n=1000L) {
  labels <- vapply(regions,function(r) switch(r,'near-truth'=rep(1L,n),mirror=rep(2L,n),
    both=rep(1:2,length.out=n),stop('unknown test region ',r)),integer(n))
  anchored_regions(named_result('anchored',unname(labels),
    extra=list(atypical=data.frame(chain=seq_along(regions),share=far))))
}

# Species-6 anatomy rows: every species-6 quantity for each chain.
anatomy_of <- function(n_chains,split=1.001,rhat=1.001,quantities=SPECIES6_QUANTITIES) {
  g <- expand.grid(chain=seq_len(n_chains),quantity=quantities,stringsAsFactors=FALSE)
  data.frame(key='synthetic',species=6L,chain=g$chain,quantity=g$quantity,
    split_rhat_within_chain=split,rhat_all_chains=rhat,stringsAsFactors=FALSE)
}
set_rows <- function(a,quantity,column,value,chain=NULL) {
  hit <- a$quantity==quantity & (if(is.null(chain)) TRUE else a$chain %in% chain)
  a[[column]][hit] <- value;a
}

# chain_separation()$separation rows for species 6.
separation_of <- function(label='agrees',rhat=1.002,fixed=FALSE,quantities=SPECIES6_QUANTITIES)
  data.frame(key='synthetic',species=6L,quantity=quantities,label=label,separation=.4,
    max_split_rhat_within_chain=1.001,rhat=rhat,fixed=fixed,stringsAsFactors=FALSE)
set_quantity <- function(s,quantity,column,value) {s[[column]][s$quantity==quantity] <- value;s}

# ---- Extended run -------------------------------------------------------------------

test_that('eight chains in each region, every split-Rhat at most 1.05, are Separated modes', {
  reg <- regions_of(rep(c('near-truth','mirror'),8))
  a <- anatomy_of(16,split=1.04,rhat=1.4)
  v <- extended_verdict(reg,a)
  expect_identical(v$verdict,'Separated modes')
  expect_identical(v$pattern,'each chain in one region')
  expect_identical(c(v$n_near_truth,v$n_mirror,v$n_both,v$n_unknown),c(8L,8L,0L,0L))
  expect_equal(v$max_split_rhat_within_chain,1.04)
  expect_true(v$split_rhat_ok)
  expect_false(v$rhat_ok)
  # One chain's own split-Rhat above 1.05 for one quantity makes it Mixed.
  v2 <- extended_verdict(reg,set_rows(a,'B0','split_rhat_within_chain',1.06,chain=5L))
  expect_identical(v2$verdict,'Mixed')
  expect_identical(v2$max_split_chain,5L)
  expect_identical(v2$max_split_quantity,'B0')
  # Uneven splits count too: one chain in the mirror region.
  expect_identical(extended_verdict(regions_of(c(rep('near-truth',15),'mirror')),anatomy_of(16))$verdict,'Separated modes')
})

test_that('every chain visiting both with all-chain Rhat at most 1.01 is Slow mixing', {
  reg <- regions_of(rep('both',16))
  a <- anatomy_of(16,split=1.2,rhat=1.008)
  v <- extended_verdict(reg,a)
  expect_identical(v$verdict,'Slow mixing')
  expect_identical(v$pattern,'every chain visits both')
  expect_true(v$rhat_ok)
  expect_equal(v$max_rhat_all_chains,1.008)
  v2 <- extended_verdict(reg,set_rows(a,'theta0','rhat_all_chains',1.02))
  expect_identical(v2$verdict,'Mixed')
  expect_identical(v2$max_rhat_quantity,'theta0')
  # A chain that stays in one region while the others visit both is Mixed.
  expect_identical(extended_verdict(regions_of(c(rep('both',15),'mirror')),anatomy_of(16))$verdict,'Mixed')
})

test_that('one mode only names the mode, whatever the Rhats', {
  v <- extended_verdict(regions_of(rep('near-truth',16)),anatomy_of(16,split=1.3,rhat=1.3))
  expect_identical(v$verdict,'one mode only: near-truth')
  expect_identical(v$n_near_truth,16L)
  expect_identical(extended_verdict(regions_of(rep('mirror',16)),anatomy_of(16))$verdict,'one mode only: mirror')
})

test_that('a chain in an unknown region makes the pattern other and the verdict Mixed', {
  far <- c(rep(.005,15),.2)
  v <- extended_verdict(regions_of(rep('near-truth',16),far=far),anatomy_of(16))
  expect_identical(v$pattern,'other')
  expect_identical(v$verdict,'Mixed')
  expect_identical(v$n_unknown,1L)
  expect_match(v$note,'chain 16 in an unknown region')
  # At 5% exactly the chain is not unknown (the rule is more than 5%).
  expect_identical(extended_verdict(regions_of(rep('near-truth',16),far=c(rep(.005,15),.05)),anatomy_of(16))$verdict,
    'one mode only: near-truth')
})

test_that('the extended verdict refuses incomplete anatomy', {
  reg <- regions_of(rep(c('near-truth','mirror'),8))
  a <- anatomy_of(16)
  expect_error(extended_verdict(reg,a[a$quantity!='B0',]),'lacks species-6 quantities: B0')
  expect_error(extended_verdict(reg,a[a$chain!=16L,]),'chains')
  expect_error(extended_verdict(reg,rbind(a,a[1,])),'one row per chain')
  # Loadings are ignored if present.
  extra <- rbind(a,transform(anatomy_of(16,split=2,rhat=2,quantities='L1')))
  expect_identical(extended_verdict(reg,extra)$verdict,'Separated modes')
})

# ---- Variants (a) and (b) -------------------------------------------------------------

test_that('a variant explains the split when every species-6 quantity agrees', {
  e <- explains_verdict(separation_of())
  expect_identical(e$verdict$verdict,'explains')
  expect_true(e$verdict$explains)
  expect_identical(e$verdict$n_considered,11L)
  expect_identical(e$verdict$n_agree,11L)
  expect_false(e$verdict$uninformative)
  expect_identical(nrow(e$quantities),11L)
})

test_that('agrees needs the label agrees and all-chain Rhat at most 1.05 (R14)', {
  # Species 6 B0 in the pr11 fit: labelled agrees with Rhat 1.071.
  e <- explains_verdict(set_quantity(separation_of(),'B0','rhat',1.071))
  expect_identical(e$verdict$verdict,'does not explain')
  expect_identical(e$verdict$not_agreeing,'B0')
  expect_false(e$quantities$agrees[e$quantities$quantity=='B0'])
  expect_identical(explains_verdict(set_quantity(separation_of(),'B0','rhat',1.05))$verdict$verdict,'explains')
  expect_identical(explains_verdict(set_quantity(separation_of(),'B_slope1','label','separated'))$verdict$verdict,'does not explain')
  expect_identical(explains_verdict(set_quantity(separation_of(),'theta0','label','drifting'))$verdict$verdict,'does not explain')
  # A missing Rhat does not count as agreeing unless the quantity is fixed.
  expect_identical(explains_verdict(set_quantity(separation_of(),'q_primer1','rhat',NA_real_))$verdict$verdict,'does not explain')
})

test_that('a quantity fixed by design is left out (theta0 in variant b)', {
  s <- set_quantity(set_quantity(separation_of(),'theta0','fixed',TRUE),'theta0','rhat',NA_real_)
  e <- explains_verdict(s)
  expect_identical(e$verdict$verdict,'explains')
  expect_identical(e$verdict$n_considered,10L)
  expect_identical(e$verdict$excluded_fixed,'theta0')
  expect_false(e$quantities$considered[e$quantities$quantity=='theta0'])
})

test_that('variant verdicts are flagged uninformative when the extended run shows one mode', {
  e <- explains_verdict(separation_of(),extended='one mode only: mirror')
  expect_identical(e$verdict$verdict,'explains')
  expect_true(e$verdict$uninformative)
  expect_false(explains_verdict(separation_of(),extended='Separated modes')$verdict$uninformative)
  expect_error(explains_verdict(separation_of(quantities=SPECIES6_QUANTITIES[-1])),'lacks species-6 quantities: B0')
})

# ---- Variant (c) ------------------------------------------------------------------------

changes_of <- function(change,mcse,names=paste0('OTU_',c(1:5,7:10)))
  data.frame(species_name=names,change=change,combined_mcse=mcse,stringsAsFactors=FALSE)

test_that('a change below 0.01 isolates whatever its standard error', {
  x <- isolation_checks(changes_of(.008,1e-4,'OTU_1'))
  expect_true(x$below_change_limit);expect_false(x$within_2_mcse);expect_true(x$pass)
  expect_true(isolation_checks(changes_of(-.0099,1e-4,'OTU_1'))$pass)
  expect_false(isolation_checks(changes_of(.01,1e-4,'OTU_1'))$below_change_limit)
})

test_that('a change within 2 combined MCSEs isolates although above 0.01', {
  x <- isolation_checks(changes_of(.03,.02,'OTU_7'))
  expect_false(x$below_change_limit);expect_true(x$within_2_mcse);expect_true(x$pass)
  y <- isolation_checks(changes_of(-.03,.01,'OTU_7'))
  expect_false(y$within_2_mcse);expect_false(y$pass)
})

test_that('variant (c) isolates only if all nine species pass', {
  ok <- changes_of(c(.001,-.004,.03,.002,0,.009,-.002,.006,.005),c(rep(.001,2),.02,rep(.001,6)))
  v <- isolates_verdict(isolation_checks(ok))
  expect_identical(v$verdict,'isolates')
  expect_identical(c(v$n_species,v$n_pass,v$n_by_change,v$n_by_mcse_only),c(9L,9L,8L,1L))
  expect_equal(v$max_abs_change,.03)
  bad <- ok;bad$combined_mcse[3] <- .01
  v2 <- isolates_verdict(isolation_checks(bad),extended='Separated modes')
  expect_identical(v2$verdict,'does not isolate')
  expect_identical(v2$failing,'OTU_3')
  expect_false(v2$uninformative)
  expect_true(isolates_verdict(isolation_checks(ok),extended='one mode only: near-truth')$uninformative)
})

test_that('occupancy_change pools chains, matches species by name and combines MCSEs in quadrature', {
  set.seed(501)
  mk <- function(m,K) matrix(rnorm(2000L*K,m,.05),2000L,K)
  ext <- list(OTU_1=mk(.4,4),OTU_6=mk(.5,4),OTU_7=mk(.6,4))
  var <- list(OTU_7=mk(.63,2),OTU_1=mk(.4,2))
  x <- occupancy_change(ext,var)
  expect_identical(x$species_name,c('OTU_7','OTU_1'))
  expect_equal(x$extended_mean,c(mean(ext$OTU_7),mean(ext$OTU_1)))
  expect_equal(x$change,c(mean(var$OTU_7)-mean(ext$OTU_7),mean(var$OTU_1)-mean(ext$OTU_1)))
  expect_equal(x$extended_mcse[1],posterior::mcse_mean(ext$OTU_7))
  expect_equal(x$combined_mcse,sqrt(x$extended_mcse^2+x$variant_mcse^2))
  expect_true(all(c('below_change_limit','within_2_mcse','pass') %in% names(x)))
  expect_false(x$pass[1]);expect_true(x$pass[2])
  expect_error(occupancy_change(ext,list(OTU_2=mk(.4,2))),'not in the extended run: OTU_2')
})
