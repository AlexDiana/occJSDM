#!/usr/bin/env Rscript
# Task 3c of PLAN.md in this directory: apply the frozen decision rules
# (README.md as amended by AMENDMENT-1.md) to the 40 diagnostic chains of
# design-qfar_K6-sites300-05 (community 5). It reads the saved fits and inputs
# and never fits a model; the archive is read, never written.
#
#   Rscript analyse.R [--study=ARCHIVE] [--archives=DIR] [--inputs-root=DIR]
#     [--out=DIR] [--chains=N]
#
# Defaults: ARCHIVE is dev/simstudy/results/convergence-diagnosis-20261001 of
# the main checkout, --archives its parent, --inputs-root
# ARCHIVES/intercept-prior-inputs, --out this directory's results/. --chains=N
# uses only the first N chains of each run: a plumbing check, written to an
# --out outside results/ and never read as a result.
#
# Steps, in the order the protocol gives them:
#   0. Refuse unless modes.R, anatomy.R and the anchored classifier have the
#      md5s AMENDMENT-1 records, and run.R the md5 the launcher froze (its
#      functions give the frozen seeds, schedules and variant definitions).
#      Check every fit's metadata against the frozen run, chain, seed and
#      schedule, and the inputs against their frozen md5s (variant (c): the
#      saved derived input, also identical to a fresh drop_species()).
#   1. chain_anatomy() and chain_separation() of the extended run (16 files)
#      and of each variant (8 files; variant (c) with its nine species).
#   2. Species 6 in the extended run and variants (a) and (b): the primary
#      anchored assignment (R17), the far-from-both share and unknown rule
#      (R18), and the secondary checks (theta0 cut where theta0 is free, the
#      refitted mixture), compared per chain.
#   3. The verdicts of verdicts.R; posterior mass per mode; species 6's means
#      per region against the generating values.
#   4. results/extended/ and results/variants/ (compact CSVs and figures).

`%||%` <- function(x,y) if(is.null(x)) y else x

HERE <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
args <- commandArgs(trailingOnly=TRUE)
opt <- function(name,default=NULL) {
  hit <- grep(paste0('^--',name,'='),args,value=TRUE)
  if(length(hit)) sub(paste0('^--',name,'='),'',hit[length(hit)]) else default
}
known <- c('study','archives','inputs-root','out','chains')
bad <- setdiff(sub('=.*$','',sub('^--','',args)),known)
if(length(bad)) stop('Unknown option: --',bad[1])

# ---- 0. Frozen files -------------------------------------------------------------

FROZEN_MD5 <- c('modes.R'='bab7c1c3934ba2e0c83b8f26232a0aaf','anatomy.R'='bd98ff3f40cad972f96bdeca5a07e206',
  'results/modes/anchor-classifier.csv'='e235c2fa641eb05c36232bec0d513bc8',
  'run.R'='f9e167d8cc8879e27b32f80a4135eb26')
found <- unname(tools::md5sum(file.path(HERE,names(FROZEN_MD5))))
if(!identical(found,unname(FROZEN_MD5))) stop('BLOCKED: frozen file md5 differs: ',
  paste(names(FROZEN_MD5)[found!=FROZEN_MD5],collapse=', '))
cat(format(Sys.time()),'frozen md5s match:',paste(names(FROZEN_MD5),FROZEN_MD5,collapse='; '),'\n')

RUN <- new.env(parent=globalenv());source(file.path(HERE,'run.R'),local=RUN)
source(file.path(HERE,'anatomy.R'))
source(file.path(HERE,'modes.R'))
source(file.path(HERE,'verdicts.R'))
invisible(gc(reset=TRUE))

STUDY <- normalizePath(opt('study','/Users/douglasyu/src/occJSDM/dev/simstudy/results/convergence-diagnosis-20261001'),mustWork=TRUE)
ARCHIVES_DIR <- normalizePath(opt('archives',dirname(STUDY)),mustWork=TRUE)
INPUTS_ROOT <- normalizePath(opt('inputs-root',file.path(ARCHIVES_DIR,'intercept-prior-inputs')),mustWork=TRUE)
OUT <- opt('out',file.path(HERE,'results'))
N_CHAINS <- opt('chains')
if(!is.null(N_CHAINS) && normalizePath(OUT,mustWork=FALSE)==normalizePath(file.path(HERE,'results')))
  stop('--chains is a plumbing check; give an --out outside results/')

# Frozen in README.md: the community-5 input, the derived input of variant (c)
# and the clone of variant (b).
SOURCE_INPUT_MD5 <- '574ed4df46b9c1bd227c79d2185a7fcb'
DERIVED_INPUT_MD5 <- 'b2c80630eedaa0c1692a8aa9535d0571'
CLONE_MD5 <- '55725289edeb10835b5b16278f89822c'
RUNS <- list(extended=1:16,a=1:8,b=1:8,c=1:8)
if(!is.null(N_CHAINS)) RUNS <- lapply(RUNS,function(ch) ch[seq_len(min(length(ch),as.integer(N_CHAINS)))])
TARGET <- RUN$TARGET_SPECIES
STRIP_QUANTITIES <- c('theta0','B0','beta_theta_intercept','B_slope1','B_slope2','mean_psi_original_sites')
RUN_LABELS <- c(extended='extended run',a='variant (a), Beta(1, 100) theta0 prior',
  b='variant (b), species 6 theta0 fixed at 0.038',c='variant (c), species 6 removed')

src <- RUN$pr11_job(ARCHIVES_DIR,INPUTS_ROOT)
if(!identical(src$input_md5,SOURCE_INPUT_MD5)) stop('The pr11 job records input md5 ',src$input_md5)
SOURCE_INPUT <- RUN$checked_input(src$input_file,SOURCE_INPUT_MD5)
DERIVED_INPUT <- RUN$checked_input(file.path(STUDY,RUN$DERIVED_INPUT),DERIVED_INPUT_MD5)
source_input <- readRDS(SOURCE_INPUT)
if(!identical(readRDS(DERIVED_INPUT),RUN$drop_species(source_input)))
  stop('The saved derived input differs from drop_species() of the source input')
FIXED_THETA0 <- RUN$fixed_theta0_value(source_input)
SPECIES_ALL <- colnames(source_input$sim$data_list$OTU)
SPECIES_C <- setdiff(SPECIES_ALL,RUN$TARGET_SPECIES_NAME)
rm(source_input)
REVISION <- readLines(file.path(STUDY,'source-revision.txt'))
LAUNCHER_HASHES <- c('dev/simstudy/convergence-flag-diagnosis/run.R'='f9e167d8cc8879e27b32f80a4135eb26',
  'dev/simstudy/convergence-flag-diagnosis/fixed-theta0.R'='79df8885f21c206524bc8da900f1d869',
  'dev/simstudy/convergence-flag-diagnosis/results/library-fingerprint.csv'='7ef28e22ddddc076bc9438022aa419ba',
  'dev/simstudy/convergence-flag-diagnosis/results/fixed-theta0/hashes.csv'='3134069c78e0f38a2dfe43db6f529d0b')
cat(format(Sys.time()),'inputs: source',SOURCE_INPUT_MD5,'; derived',DERIVED_INPUT_MD5,
  '(identical to a fresh derivation); fixed theta0',format(FIXED_THETA0,digits=12),'\n')

run_files <- function(run) file.path(STUDY,'fits',run,sprintf('chain-%02d-fit.rds',RUNS[[run]]))
run_input <- function(run) if(run=='c') DERIVED_INPUT else SOURCE_INPUT

# Every field that identifies the frozen fit, checked before the file is used.
expected_definition <- function(run) switch(run,
  extended=list(change='none'),
  a=list(change='listPriors added',listPriors_added=RUN$VARIANT_A_PRIORS),
  b=list(change='theta0 of one species held fixed by the fixed-theta0.R clone of sample_theta0',
    species=TARGET,species_name=RUN$TARGET_SPECIES_NAME,fixed_value=FIXED_THETA0,
    fixed_value_source='input$truth$params$theta0[6]',clone_md5=CLONE_MD5),
  c=list(change='species removed from the input',species=TARGET,species_name=RUN$TARGET_SPECIES_NAME,
    source_input_md5=SOURCE_INPUT_MD5,derived_input_md5=DERIVED_INPUT_MD5,derived_input=RUN$DERIVED_INPUT))

check_fit <- function(file,run,chain) {
  s <- readRDS(file);fail <- character()
  chk <- function(ok,what) if(!isTRUE(ok)) fail <<- c(fail,what)
  input_md5 <- if(run=='c') DERIVED_INPUT_MD5 else SOURCE_INPUT_MD5
  chk(identical(s$key,RUN$KEY),'key')
  chk(identical(s$run,run) && identical(s$variant,run),'run')
  chk(identical(s$chain,as.integer(chain)),'chain')
  chk(identical(s$seed,RUN$chain_seed(chain)),'seed')
  chk(identical(s$seed_rule,RUN$SEED_RULE) && identical(s$rng_kind,RUN$RNG_KIND),'seed rule')
  set.seed(RUN$chain_seed(chain),kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
  chk(identical(s$rng_initial,.Random.seed),'initial RNG state')
  chk(identical(s$schedule,'diagnostic') && identical(s$mcmc,RUN$SCHEDULES$diagnostic),'schedule')
  chk(identical(s$input_md5,input_md5) && identical(s$job$input_md5,input_md5),'input md5')
  chk(identical(s$input_md5_after,input_md5) && identical(s$source_input_md5_after,SOURCE_INPUT_MD5),'input unchanged')
  chk(identical(s$source_revision,REVISION),'library revision')
  chk(identical(s$script_hashes[names(LAUNCHER_HASHES)],LAUNCHER_HASHES),'script hashes')
  chk(identical(s$listPriors_used,c(src$job$priors,if(run=='a') RUN$VARIANT_A_PRIORS else list())),'listPriors')
  chk(identical(s$variant_definition,expected_definition(run)),'variant definition')
  chk(identical(s$iterations_run,50000),'iterations run')
  species <- s$fit$infos$speciesNames;th <- s$fit$results_output$theta0_output
  chk(identical(species,if(run=='c') SPECIES_C else SPECIES_ALL),'species')
  chk(identical(dim(th),c(length(species),10000L,1L)),'retained draws')
  varies <- apply(th[,,1,drop=FALSE],1,function(v) length(unique(v))>1L)
  if(run=='b') chk(all(th[TARGET,,1]==FIXED_THETA0) && all(varies[-TARGET]),'theta0 fixed for species 6 only') else
    chk(all(varies),'theta0 varies')
  row <- data.frame(run=run,chain=as.integer(chain),fit_file=relative_to(file,ARCHIVES_DIR),
    fit_md5=unname(tools::md5sum(file)),seed=s$seed,schedule=paste(names(s$mcmc),unlist(s$mcmc),collapse=' '),
    input_md5=s$input_md5,n_species=length(species),warnings=length(s$warnings),
    started=format(s$started),finished=format(s$finished),elapsed_seconds=round(s$elapsed_seconds,1),
    checks_failed=paste(fail,collapse=';'),stringsAsFactors=FALSE)
  rm(s);invisible(gc(FALSE))
  row
}

provenance <- list()
for(run in names(RUNS)) for(k in seq_along(RUNS[[run]]))
  provenance[[length(provenance)+1L]] <- check_fit(run_files(run)[k],run,RUNS[[run]][k])
provenance <- do.call(rbind,provenance)
if(any(nzchar(provenance$checks_failed))) {
  print(provenance[nzchar(provenance$checks_failed),c('run','chain','checks_failed')],row.names=FALSE)
  stop('Fit metadata differ from the frozen runs')
}
cat(format(Sys.time()),'metadata of',nrow(provenance),'fits match the frozen runs; package warnings',
  sum(provenance$warnings),'\n')

# ---- Helpers -------------------------------------------------------------------------

anchor <- read_anchor(file.path(HERE,ANCHOR_FILE))
dir.create(file.path(OUT,'extended'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUT,'variants'),recursive=TRUE,showWarnings=FALSE)

named_anatomy <- function(a,names) {
  a$species_name <- names[a$species]
  a[c('key','species','species_name',setdiff(names(a),c('key','species','species_name')))]
}
with_names <- function(x,names) {
  if(!nrow(x)) return(x)
  x$species_name <- names[x$species]
  x[c('key','species','species_name',setdiff(names(x),c('key','species','species_name')))]
}

# Anatomy of one run, keeping the species-6 draws of STRIP_QUANTITIES and every
# species' mean_psi_original_sites draws; the full draws are dropped.
run_anatomy_draws <- function(run) {
  started <- Sys.time()
  a <- chain_anatomy(run_files(run),run_input(run),keep_draws=TRUE)
  names <- if(run=='c') SPECIES_C else SPECIES_ALL
  draws <- attr(a,'draws')
  d6 <- if(run=='c') NULL else species_mode_draws(a,TARGET,STRIP_QUANTITIES)
  psi <- stats::setNames(lapply(seq_along(names),function(s) matrix(draws$mean_psi_original_sites[s,,],
    dim(draws[[1]])[2],dim(draws[[1]])[3])),names)
  attr(a,'draws') <- NULL;rm(draws);invisible(gc(FALSE))
  s <- chain_separation(a)
  cat(format(Sys.time()),run,sprintf('anatomy of %d chains in %.0f s; psi check %.2e; species labels: %s',
    length(RUNS[[run]]),as.numeric(difftime(Sys.time(),started,units='secs')),max(attr(a,'psi_check')$psi_max_abs_diff),
    paste0(names[s$species$species],' ',s$species$label,collapse=', ')),'\n')
  flush.console()
  list(run=run,anatomy=a,separation=s,names=names,d6=d6,psi=psi,psi_mean=attr(a,'psi_mean'))
}

write_anatomy <- function(r,dir) {
  write_compact(named_anatomy(r$anatomy,r$names),file.path(dir,'chain-summary.csv'))
  write_compact(with_names(r$separation$separation,r$names),file.path(dir,'separation.csv'))
  write_compact(with_names(r$separation$species,r$names),file.path(dir,'species-labels.csv'))
  write_compact(with_names(r$separation$chain_groups,r$names),file.path(dir,'chain-groups.csv'))
}

truth6 <- function(r) {
  x <- r$anatomy[r$anatomy$species==TARGET,,drop=FALSE]
  stats::setNames(x$truth[match(STRIP_QUANTITIES,x$quantity)],STRIP_QUANTITIES)
}

summary_rows <- function(d6,pick,truth,run,grouping,group,chains) do.call(rbind,lapply(STRIP_QUANTITIES,function(q) {
  x <- d6[[q]][pick];qs <- stats::quantile(x,c(.025,.975),names=FALSE)
  data.frame(run=run,species=TARGET,grouping=grouping,group=group,chains=chains,n_draws=length(x),quantity=q,
    mean=mean(x),sd=stats::sd(x),q025=qs[1],q975=qs[2],truth=truth[[q]],error=mean(x)-truth[[q]],
    abs_error=abs(mean(x)-truth[[q]]),truth_in_interval=truth[[q]]>=qs[1] & truth[[q]]<=qs[2],stringsAsFactors=FALSE)
}))

# Species 6's means per anchored region (chains grouped by their region) and
# per primary draw label (all chains), against the generating values.
region_means <- function(d6,primary,regions,truth,run) {
  ni <- nrow(d6[[1]]);K <- ncol(d6[[1]]);rows <- list()
  for(g in intersect(c('near-truth','mirror','both','unknown'),regions$region)) {
    ch <- regions$chain[regions$region==g]
    pick <- matrix(FALSE,ni,K);pick[,ch] <- TRUE
    rows[[length(rows)+1L]] <- summary_rows(d6,pick,truth,run,'chain region',g,paste(ch,collapse=';'))
  }
  for(m in 1:2) if(any(primary$labels==m))
    rows[[length(rows)+1L]] <- summary_rows(d6,primary$labels==m,truth,run,'draw label',ANCHOR_MODES[m],'all')
  do.call(rbind,rows)
}

# Posterior mass per mode over all draws, by assignment method.
overall_mass <- function(run,primary,cut,refit) {
  ni <- nrow(primary$labels)
  far <- data.frame(mode=NA_integer_,mode_name='far from both (R18)',draws=as.integer(round(sum(primary$atypical$share*ni))),
    share=mean(primary$atypical$share))
  one <- function(method,r,names,n_modes=NA_integer_) data.frame(run=run,method=method,n_modes=n_modes,
    mode=r$overall$mode,mode_name=names[r$overall$mode],draws=r$overall$draws,share=r$overall$share,stringsAsFactors=FALSE)
  out <- rbind(one('anchored',primary,ANCHOR_MODES),data.frame(run=run,method='anchored',n_modes=NA_integer_,far,stringsAsFactors=FALSE))
  if(!is.null(cut)) out <- rbind(out,one('theta0 cut',cut,ANCHOR_MODES))
  coded <- refit_in_anchor_names(refit,primary)
  out <- rbind(out,one('refitted mixture (named by majority overlap with the anchored labels)',
    named_result('refit',coded),ANCHOR_MODES,refit$n_modes))
  out
}

species6_modes <- function(r) {
  run <- r$run;d6 <- r$d6
  primary <- assign_anchored(d6,anchor)
  regions <- anchored_regions(primary)
  theta0_free <- stats::sd(as.vector(d6$theta0))>0
  cut <- if(theta0_free) theta0_cut_assignment(d6) else NULL
  refit <- assign_modes(d6[MODE_QUANTITIES])
  cmp <- compare_assignments(primary,cut,refit)
  agreement <- attr(cmp,'draw_agreement')
  cat(format(Sys.time()),run,'species 6: anchored regions',paste(regions$region,collapse=','),
    '; pattern',region_pattern(regions),'; max far share',sprintf('%.4f',max(regions$share_far)),
    '; refit modes',refit$n_modes,'(dropped:',paste(refit$dropped_constant,collapse=',') ,')',
    '; draw agreement',paste(names(agreement),sprintf('%.5f',agreement),collapse=' '),'\n')
  flush.console()
  list(primary=primary,regions=regions,cut=cut,refit=refit,comparison=cmp,agreement=agreement,
    mass=mode_mass_table(primary,run),overall=overall_mass(run,primary,cut,refit),
    means=region_means(d6,primary,regions,truth6(r),run))
}

write_modes <- function(r,m,dir) {
  write_compact(m$mass,file.path(dir,'mode-mass.csv'))
  write_compact(m$overall,file.path(dir,'mode-mass-overall.csv'))
  cmp <- m$comparison;cmp <- data.frame(run=r$run,cmp,stringsAsFactors=FALSE)
  write_compact(cmp,file.path(dir,'assignment-comparison.csv'))
  write_compact(data.frame(run=r$run,method=names(m$agreement),share_of_draws_labelled_as_primary=unname(m$agreement)),
    file.path(dir,'assignment-draw-agreement.csv'))
  write_compact(m$means,file.path(dir,'region-means-vs-truth.csv'))
  K <- ncol(r$d6[[1]]);ni <- nrow(r$d6[[1]])
  plot_chain_strips(r$d6,file.path(dir,'chain-strips.png'),quantities=STRIP_QUANTITIES,truth=truth6(r),modes=m$primary,
    title=sprintf('%s, species 6 (OTU_6), %s: %d chains of %d retained draws (thin line 90%%, thick 50%%, dot mean)',
      RUN$KEY,RUN_LABELS[[r$run]],K,ni))
}

source_hashes <- function() {
  files <- file.path(HERE,c('analyse.R','verdicts.R','test-verdicts.R','anatomy.R','modes.R',ANCHOR_FILE,'run.R'))
  data.frame(file=vapply(files,relative_to,'',root=REPO),md5=unname(tools::md5sum(files)),row.names=NULL)
}

# ---- 1-3. Extended run -----------------------------------------------------------------

ext <- run_anatomy_draws('extended')
ext_dir <- file.path(OUT,'extended')
write_anatomy(ext,ext_dir)
ext_modes <- species6_modes(ext)
write_modes(ext,ext_modes,ext_dir)
ext_verdict <- extended_verdict(ext_modes$regions,ext$anatomy[ext$anatomy$species==TARGET,,drop=FALSE])
share_of <- function(o,name) o$share[o$method=='anchored' & o$mode_name==name]
ext_verdict <- data.frame(run='extended',ext_verdict,
  share_near_truth=share_of(ext_modes$overall,'near-truth'),share_mirror=share_of(ext_modes$overall,'mirror'),
  share_far=share_of(ext_modes$overall,'far from both (R18)'),
  chains_near_truth=paste(ext_modes$regions$chain[ext_modes$regions$region=='near-truth'],collapse=';'),
  chains_mirror=paste(ext_modes$regions$chain[ext_modes$regions$region=='mirror'],collapse=';'),
  chains_both=paste(ext_modes$regions$chain[ext_modes$regions$region=='both'],collapse=';'),
  chains_unknown=paste(ext_modes$regions$chain[ext_modes$regions$region=='unknown'],collapse=';'),
  rule='AMENDMENT-1, extended-run verdict: Slow mixing if every chain visits both anchored regions and rhat_all_chains <= 1.01 for the eleven species-6 quantities; Separated modes if every chain stays in one region, each region holds a chain, and every split_rhat_within_chain <= 1.05; one mode only (R15) if all chains stay in the same region; Mixed otherwise, including any chain in an unknown region (R18)',
  stringsAsFactors=FALSE)
write_compact(ext_verdict,file.path(ext_dir,'verdict.csv'))
write_compact(provenance[provenance$run=='extended',,drop=FALSE],file.path(ext_dir,'provenance.csv'))
utils::write.csv(source_hashes(),file.path(ext_dir,'source-hashes.csv'),row.names=FALSE)
EXTENDED_VERDICT <- ext_verdict$verdict
cat(format(Sys.time()),'EXTENDED VERDICT:',EXTENDED_VERDICT,'; pattern',ext_verdict$pattern,
  sprintf('; max all-chain Rhat %.4f (%s); max chain split-Rhat %.4f (chain %d, %s)',ext_verdict$max_rhat_all_chains,
    ext_verdict$max_rhat_quantity,ext_verdict$max_split_rhat_within_chain,ext_verdict$max_split_chain,
    ext_verdict$max_split_quantity),'\n')
flush.console()
ext_psi <- ext$psi;ext_psi_mean <- ext$psi_mean
rm(ext);invisible(gc(FALSE))

# ---- Variants (a) and (b) ------------------------------------------------------------------

verdict_rows <- list()
for(run in c('a','b')) {
  r <- run_anatomy_draws(run)
  dir <- file.path(OUT,'variants',run);dir.create(dir,showWarnings=FALSE)
  write_anatomy(r,dir)
  m <- species6_modes(r)
  write_modes(r,m,dir)
  sep6 <- r$separation$separation[r$separation$separation$species==TARGET,,drop=FALSE]
  e <- explains_verdict(sep6,EXTENDED_VERDICT)
  write_compact(data.frame(run=run,e$quantities,stringsAsFactors=FALSE),file.path(dir,'species6-agreement.csv'))
  reg <- m$regions
  beside <- sprintf('beside the verdict (not deciding): anchored regions %s (pattern %s%s); refitted mixture %d mode(s)%s',
    paste(names(table(reg$region)),table(reg$region),sep=' x',collapse=', '),region_pattern(reg),
    if(nzchar(region_note(reg))) paste0('; ',region_note(reg)) else '',m$refit$n_modes,
    if(is.null(m$cut)) '' else sprintf('; theta0 cut: %d of %d chains agree with the primary region',sum(m$comparison$agree_cut),nrow(reg)))
  v <- e$verdict
  verdict_rows[[length(verdict_rows)+1L]] <- data.frame(variant=run,criterion='explains',item='all species-6 quantities',
    result=v$verdict,pass=v$explains,label=NA_character_,separation=v$max_separation,rhat=v$max_rhat,fixed=NA,
    change=NA_real_,combined_mcse=NA_real_,below_change_limit=NA,within_2_mcse=NA,uninformative=v$uninformative,
    extended_verdict=EXTENDED_VERDICT,
    detail=sprintf('%d of %d quantities agree (label agrees and all-chain Rhat <= 1.05)%s%s; %s',v$n_agree,v$n_considered,
      if(nzchar(v$not_agreeing)) paste0('; not agreeing: ',v$not_agreeing) else '',
      if(nzchar(v$excluded_fixed)) paste0('; left out as fixed by design: ',v$excluded_fixed) else '',beside),
    stringsAsFactors=FALSE)
  q <- e$quantities
  verdict_rows[[length(verdict_rows)+1L]] <- data.frame(variant=run,criterion='species-6 quantity agrees',item=q$quantity,
    result=ifelse(q$fixed,'left out (fixed by design)',ifelse(q$agrees,'agrees','does not agree')),
    pass=ifelse(q$fixed,NA,q$agrees),label=q$label,separation=q$separation,rhat=q$rhat,fixed=q$fixed,
    change=NA_real_,combined_mcse=NA_real_,below_change_limit=NA,within_2_mcse=NA,uninformative=v$uninformative,
    extended_verdict=EXTENDED_VERDICT,detail=sprintf('max chain split-Rhat %.4f',q$max_split_rhat_within_chain),
    stringsAsFactors=FALSE)
  write_compact(provenance[provenance$run==run,,drop=FALSE],file.path(dir,'provenance.csv'))
  cat(format(Sys.time()),'VARIANT',run,'VERDICT:',v$verdict,sprintf('(%d of %d agree%s)',v$n_agree,v$n_considered,
    if(nzchar(v$not_agreeing)) paste0('; not agreeing: ',v$not_agreeing) else ''),'\n')
  flush.console()
  rm(r,m);invisible(gc(FALSE))
}

# ---- Variant (c) -------------------------------------------------------------------------------

r <- run_anatomy_draws('c')
dir <- file.path(OUT,'variants','c');dir.create(dir,showWarnings=FALSE)
write_anatomy(r,dir)
changes <- occupancy_change(ext_psi,r$psi)
site <- vapply(changes$species_name,function(s) {
  d <- abs(r$psi_mean[,match(s,r$names)]-ext_psi_mean[,match(s,SPECIES_ALL)]);c(max(d),which.max(d))
},numeric(2))
truth_c <- r$anatomy[r$anatomy$quantity=='mean_psi_original_sites' & r$anatomy$chain==1L,,drop=FALSE]
changes <- data.frame(run='c',changes,truth=truth_c$truth[match(changes$species_name,r$names[truth_c$species])],
  max_abs_site_change=site[1,],site_of_max_change=as.integer(site[2,]),stringsAsFactors=FALSE)
write_compact(changes,file.path(dir,'occupancy-change.csv'))
iv <- isolates_verdict(changes,EXTENDED_VERDICT)
verdict_rows[[length(verdict_rows)+1L]] <- data.frame(variant='c',criterion='isolates',item='nine other species',
  result=iv$verdict,pass=iv$isolates,label=NA_character_,separation=NA_real_,rhat=NA_real_,fixed=NA,
  change=changes$change[which.max(abs(changes$change))],combined_mcse=changes$combined_mcse[which.max(abs(changes$change))],
  below_change_limit=NA,within_2_mcse=NA,uninformative=iv$uninformative,extended_verdict=EXTENDED_VERDICT,
  detail=sprintf('%d of %d species pass (%d by |change| < 0.01, %d only by 2 combined MCSEs)%s; largest |change| %.4f (%s); largest site change %.4f (%s, not deciding)',
    iv$n_pass,iv$n_species,iv$n_by_change,iv$n_by_mcse_only,if(nzchar(iv$failing)) paste0('; failing: ',iv$failing) else '',
    iv$max_abs_change,changes$species_name[which.max(abs(changes$change))],max(changes$max_abs_site_change),
    changes$species_name[which.max(changes$max_abs_site_change)]),stringsAsFactors=FALSE)
verdict_rows[[length(verdict_rows)+1L]] <- data.frame(variant='c',criterion='occupancy mean change within limits',
  item=changes$species_name,result=ifelse(changes$pass,'passes','fails'),pass=changes$pass,label=NA_character_,
  separation=NA_real_,rhat=NA_real_,fixed=NA,change=changes$change,combined_mcse=changes$combined_mcse,
  below_change_limit=changes$below_change_limit,within_2_mcse=changes$within_2_mcse,uninformative=iv$uninformative,
  extended_verdict=EXTENDED_VERDICT,detail=sprintf('extended mean %.4f, variant (c) mean %.4f, largest site change %.4f',
    changes$extended_mean,changes$variant_mean,changes$max_abs_site_change),stringsAsFactors=FALSE)
psi_truth <- stats::setNames(changes$truth,changes$species_name)
plot_chain_strips(r$psi[changes$species_name],file.path(dir,'chain-strips.png'),truth=psi_truth,
  title=sprintf('%s, %s: mean occupancy over the original sites, %d chains of %d retained draws',
    RUN$KEY,RUN_LABELS[['c']],ncol(r$psi[[1]]),nrow(r$psi[[1]])))
write_compact(provenance[provenance$run=='c',,drop=FALSE],file.path(dir,'provenance.csv'))
cat(format(Sys.time()),'VARIANT c VERDICT:',iv$verdict,sprintf('(%d of %d pass; largest |change| %.4f)',iv$n_pass,iv$n_species,iv$max_abs_change),'\n')
rm(r);invisible(gc(FALSE))

verdicts <- do.call(rbind,verdict_rows);rownames(verdicts) <- NULL
write_compact(verdicts,file.path(OUT,'variants','verdicts.csv'))
utils::write.csv(source_hashes(),file.path(OUT,'variants','source-hashes.csv'),row.names=FALSE)
print(verdicts[verdicts$criterion %in% c('explains','isolates'),c('variant','criterion','result','uninformative','detail')],row.names=FALSE)
mem <- gc()
cat(format(Sys.time()),sprintf('done; R max memory used %.0f Mb (Ncells %.0f Mb, Vcells %.0f Mb)',sum(mem[,ncol(mem)]),mem[1,ncol(mem)],mem[2,ncol(mem)]),'\n')
