#!/usr/bin/env Rscript
# Setup and verification for the Task 3 diagnostic fits (README.md in this
# directory). run.R neither sources nor hashes this file.
#
#   Rscript verify.R --mode=MODE --repo=REPO --study=ARCHIVE [--inputs-root=DIR]
#     [--archives=DIR] [--side=pr11|study]
#
# Modes:
#   fingerprint     write results/library-fingerprint.csv for ARCHIVE's build
#                   (revision, exported source and installed library md5s);
#   clone-archive   write results/fixed-theta0/ (the installed sample_theta0,
#                   the clone and their md5s);
#   equivalence     fit a short default fit of community 5 with the pr11
#                   library and with the study library, each in its own
#                   process (equivalence-fit --side=...), then compare
#                   (equivalence-compare); exit 1 unless every check passes;
#   pilot-checks    check the four pilot fits under ARCHIVE/fits/pilot/ and
#                   write results/pilot-checks.csv; exit 1 unless all pass;
#   mode-calibration  run assign_modes() on every species of the saved pr11
#                   fit of community 5 (Task 1 data) and write
#                   results/mode-calibration.csv;
#   anchor-calibration  calibrate the anchored classifier of AMENDMENT-1 on the
#                   same fit (species 6, chains 1 and 3 near-truth, 2 and 4
#                   mirror), write results/modes/anchor-classifier.csv, and
#                   record its validation (anchor-validation.csv) and the
#                   pseudo-chain imbalance checks of all three assignments
#                   (imbalance-check.csv) in results/modes/.
# Pilot fits are plumbing checks only and are never analysed.

args <- commandArgs(trailingOnly=TRUE)
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(here,'run.R'))
PR11$extract_functions(file.path(SCRIPT_REPO,'dev/simstudy/occupancy-intercept-prior/verify-helpers.R'),
  c('installed_library_files','max_abs_difference','count_numeric','write_library_fingerprint'),environment())

o <- parse_options(args,known=c('mode','repo','study','inputs-root','archives','side'),required=c('mode','repo','study'))
repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this verify.R')
archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
inputs_root <- o$`inputs-root` %||% file.path(archives,'intercept-prior-inputs')
results <- file.path(repo,STUDY_REL,'results')
TOLERANCE <- 1e-12
EQ_DIR <- file.path(study,'equivalence')

# ---- Equivalence -------------------------------------------------------------------

# The pr11 initial schedule shortened to 50 burn-in and 50 retained draws,
# keeping its storage types (2 chains, nthin 1).
equivalence_mcmc <- function() {
  m <- readRDS(file.path(archives,PR11_SETTINGS[1]))$mcmc
  stopifnot(identical(names(m),c('nchain','nburn','niter','nthin')),all(unlist(m)==c(2,3000,5000,1)))
  m$nburn <- as.vector(50,mode=typeof(m$nburn));m$niter <- as.vector(50,mode=typeof(m$niter))
  m
}

equivalence_fit <- function(side) {
  if(side=='study') {
    revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
    pkg <- load_study_library(study)
  } else if(side=='pr11') {
    lib <- file.path(archives,PR11_ARCHIVE,'library')
    Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    .libPaths(c(lib,.libPaths()))
    suppressPackageStartupMessages(library(occJSDM))
    pkg <- normalizePath(find.package('occJSDM'))
    stopifnot(identical(pkg,normalizePath(file.path(lib,'occJSDM'))))
    RcppParallel::setThreadOptions(numThreads=1)
    revision <- readLines(file.path(archives,PR11_ARCHIVE,'source-revision.txt'))
  } else stop('--side must be pr11 or study')
  runner <- pr11_runner(repo,archives);src <- pr11_job(archives,inputs_root)
  checked_input(src$input_file,src$input_md5);input <- readRDS(src$input_file)
  recorder <- new.env()
  if(side=='pr11') {expr <- runner$expr;bindings <- list()} else {
    # The extended run's own path: run.R's statement and fitter.
    expr <- with_intercept_prior(runner$expr)
    bindings <- list(.fitter=make_run_fitter('extended',input,recorder))
  }
  mcmc <- equivalence_mcmc()
  result <- evaluate_seeded(expr,input,src$job,mcmc,bindings=bindings,seed='input')
  saved <- c(list(side=side,library=pkg,revision=revision,library_hashes=tools::md5sum(installed_library_files(pkg)),
    runner=runner[c('file','md5')],expression=deparse(expr,width.cutoff=500L),mcmc=mcmc,input_file=src$input_file,
    input_md5_before=src$input_md5,input_md5_after=unname(tools::md5sum(src$input_file)),
    listPriors_used=recorder$listPriors,job_priors=src$job$priors),result,
    list(elapsed_seconds=as.numeric(difftime(result$finished,result$started,units='secs'))))
  dir.create(EQ_DIR,showWarnings=FALSE)
  atomic_save(saved,file.path(EQ_DIR,paste0(side,'.rds')))
  cat(side,'saved;',round(saved$elapsed_seconds,1),'seconds;',length(saved$warnings),'warnings\n')
}

equivalence_compare <- function() {
  a <- readRDS(file.path(EQ_DIR,'pr11.rds'));b <- readRDS(file.path(EQ_DIR,'study.rds'))
  rest <- setdiff(union(names(a$fit),names(b$fit)),c('results_output','infos'))
  fields <- union(names(a$fit$infos),names(b$fit$infos))
  differ <- fields[!vapply(fields,function(f) identical(a$fit$infos[[f]],b$fit$infos[[f]]),logical(1))]
  infos <- paste0(differ,ifelse(differ %in% names(a$fit$infos),'(changed)','(added)'),collapse=';')
  so <- function(x) unname(x$library_hashes[endsWith(names(x$library_hashes),'libs/occJSDM.so')])
  row <- data.frame(key=KEY,comparison='pr11 library vs study library, default short fit, saved input RNG state',
    pr11_revision=a$revision,pr11_so_md5=so(a),study_revision=b$revision,study_so_md5=so(b),
    mcmc=paste(unlist(b$mcmc),collapse='/'),numeric_values=count_numeric(b$fit$results_output),
    max_abs_diff=max_abs_difference(a$fit$results_output,b$fit$results_output),
    bitwise_identical=identical(a$fit$results_output,b$fit$results_output),
    other_components_identical=identical(a$fit[rest],b$fit[rest]),infos_differences=infos,
    warnings_pr11=length(a$warnings),warnings_study=length(b$warnings),warnings_identical=identical(a$warnings,b$warnings),
    rng_initial_identical=identical(a$rng_initial,b$rng_initial),rng_final_identical=identical(a$rng_final,b$rng_final),
    study_listPriors_as_archived=identical(b$listPriors_used,b$job_priors),
    inputs_unchanged=all(c(a$input_md5_after,b$input_md5_after)==b$input_md5_before),
    pr11_seconds=round(a$elapsed_seconds,2),study_seconds=round(b$elapsed_seconds,2),stringsAsFactors=FALSE)
  row$pass <- row$max_abs_diff<=TOLERANCE && row$other_components_identical &&
    (!nzchar(infos) || all(strsplit(infos,';',fixed=TRUE)[[1]] %in% 'intercept_prior(added)')) &&
    row$warnings_identical && row$rng_initial_identical && row$rng_final_identical &&
    row$study_listPriors_as_archived && row$inputs_unchanged
  utils::write.csv(row,file.path(EQ_DIR,'equivalence.csv'),row.names=FALSE)
  utils::write.csv(row,file.path(results,'equivalence.csv'),row.names=FALSE)
  print(t(row))
  row$pass
}

# ---- Pilot checks ---------------------------------------------------------------------

pilot_checks <- function() {
  pkg <- load_study_library(study)
  src <- pr11_job(archives,inputs_root);input <- readRDS(checked_input(src$input_file,src$input_md5))
  v6 <- fixed_theta0_value(input)
  anatomy <- new.env(parent=globalenv());sys.source(file.path(here,'anatomy.R'),envir=anatomy)
  rows <- lapply(VARIANTS,function(v) {
    spec <- run_spec('pilot',v,1L);f <- fit_file(study,spec)
    if(!file.exists(f)) return(data.frame(variant=v,check='fit exists',pass=FALSE,value='missing',stringsAsFactors=FALSE))
    s <- readRDS(f);ro <- s$fit$results_output;th <- ro$theta0_output
    S <- dim(th)[1];species <- s$fit$infos$speciesNames
    ck <- list()
    add <- function(check,pass,value) ck[[length(ck)+1L]] <<- data.frame(variant=v,check=check,pass=isTRUE(pass),
      value=paste(value,collapse=' '),stringsAsFactors=FALSE)
    add('schedule is the pilot schedule',identical(s$mcmc,SCHEDULES$pilot),paste(names(s$mcmc),unlist(s$mcmc),collapse=' '))
    add('one chain of niter retained draws',identical(dim(th)[2:3],c(200L,1L)),paste(dim(th),collapse='x'))
    add('seed is 20261001 + chain',identical(s$seed,chain_seed(1L)),s$seed)
    add('rng state is the seeded Mersenne-Twister state',{set.seed(s$seed,kind='Mersenne-Twister',normal.kind='Inversion',
      sample.kind='Rejection');identical(s$rng_initial,.Random.seed)},s$rng_initial[1])
    add('library fingerprint revision recorded',identical(s$source_revision,readLines(file.path(study,'source-revision.txt'))),s$source_revision)
    add('fit hashes match the installed library',identical(s$fit_hashes,hash_files(production_files(study,pkg),study)),length(s$fit_hashes))
    add('script hashes include run.R and fixed-theta0.R',all(c(RUN_SCRIPT,FIXED_THETA0_SCRIPT) %in% names(s$script_hashes)),length(s$script_hashes))
    expected_priors <- c(src$job$priors,if(v=='a') VARIANT_A_PRIORS else list())
    add('listPriors passed as the variant states',identical(s$listPriors_used,expected_priors),
      paste(names(s$listPriors_used),unlist(s$listPriors_used),sep='=',collapse=','))
    if(v=='a') add('tighter theta0 prior recorded in the variant definition',
      identical(s$variant_definition$listPriors_added,list(a_theta0=1,b_theta0=100)),'a_theta0=1,b_theta0=100')
    if(v=='b') {
      add('species 6 theta0 constant at the generating value',all(th[6,,1]==v6),
        sprintf('%.10g (unique draws %d)',v6,length(unique(th[6,,1]))))
      add('other species theta0 vary',all(apply(th[-6,,1,drop=FALSE],1,function(x) length(unique(x))>1L)),S-1L)
    } else add('every species theta0 varies',all(apply(th[,,1,drop=FALSE],1,function(x) length(unique(x))>1L)),S)
    if(v=='c') {
      add('species count is 9',S==9L && length(species)==9L,S)
      add('OTU_6 absent',!'OTU_6' %in% species,paste(species,collapse=','))
      add('fit input is the saved derived input',identical(s$input_md5,unname(tools::md5sum(file.path(study,DERIVED_INPUT)))),s$input_md5)
    } else add('species count is 10',S==10L,S)
    add('warnings captured as a character vector',is.character(s$warnings),length(s$warnings))
    # Plumbing for Task 3c: Task 1's chain_anatomy() reads the fit with the
    # input it was fitted to (variant c: the derived input). No value is used.
    a <- tryCatch(anatomy$chain_anatomy(f,s$input_file),error=function(e) conditionMessage(e))
    add('chain_anatomy() reads the fit and its fitted input',is.data.frame(a) && nrow(a)==S*13L &&
      max(attr(a,'psi_check')$psi_max_abs_diff)<1e-10,if(is.data.frame(a)) paste(nrow(a),'rows') else a)
    add('seconds per 1000 iterations',is.finite(s$seconds_per_1000_iterations),round(s$seconds_per_1000_iterations,2))
    do.call(rbind,ck)
  })
  table <- do.call(rbind,rows)
  utils::write.csv(table,file.path(results,'pilot-checks.csv'),row.names=FALSE)
  print(table,row.names=FALSE)
  all(table$pass)
}

# ---- Mode calibration on the saved pr11 fit -----------------------------------------

mode_calibration <- function() {
  anatomy <- new.env(parent=globalenv());sys.source(file.path(here,'anatomy.R'),envir=anatomy)
  modes <- new.env(parent=globalenv());sys.source(file.path(here,'modes.R'),envir=modes)
  a <- anatomy$chain_anatomy(anatomy$selected_fit(KEY),anatomy$selected_input(KEY),keep_draws=TRUE)
  S <- dim(attr(a,'draws')[[1]])[1];ni <- dim(attr(a,'draws')[[1]])[2]
  rows <- lapply(seq_len(S),function(s) {
    r <- modes$assign_modes(modes$species_mode_draws(a,s));cr <- modes$chain_regions(r)
    comp <- r$components
    data.frame(key=KEY,species=s,chain=cr$chain,draws_per_chain=ni,n_modes=r$n_modes,ridgeline_maxima=r$ridgeline_maxima,
      share_mode1=cr$share_mode1,share_mode2=cr$share_mode2,visits_both=cr$visits_both,region=cr$region,
      mode1_mean_theta0=if(r$n_modes==2L) comp$mean_theta0[1] else NA_real_,
      mode2_mean_theta0=if(r$n_modes==2L) comp$mean_theta0[2] else NA_real_,stringsAsFactors=FALSE)
  })
  table <- do.call(rbind,rows)
  for(nm in names(table)) if(is.double(table[[nm]])) table[[nm]] <- signif(table[[nm]],6)
  utils::write.csv(table,file.path(results,'mode-calibration.csv'),row.names=FALSE)
  print(table[table$species==6L | table$chain==1L,],row.names=FALSE)
  TRUE
}

# ---- Anchored classifier (AMENDMENT-1, ruling R17) ------------------------------------

anchor_calibration <- function() {
  anatomy <- new.env(parent=globalenv());sys.source(file.path(here,'anatomy.R'),envir=anatomy)
  modes <- new.env(parent=globalenv());sys.source(file.path(here,'modes.R'),envir=modes)
  fit <- anatomy$selected_fit(KEY)
  a <- anatomy$chain_anatomy(fit,anatomy$selected_input(KEY),keep_draws=TRUE)
  d <- modes$species_mode_draws(a,TARGET_SPECIES,unique(c(modes$MODE_QUANTITIES,modes$ANCHOR_QUANTITIES)))
  rm(a);invisible(gc(FALSE))
  anchor <- modes$calibrate_anchor(d,near_chains=c(1L,3L),mirror_chains=c(2L,4L))
  out <- file.path(results,'modes');dir.create(out,showWarnings=FALSE)
  file <- file.path(out,'anchor-classifier.csv')
  modes$write_anchor(anchor,file,source=c(key=KEY,species=as.character(TARGET_SPECIES),
    fit=file.path(PR11_ARCHIVE,basename(dirname(fit)),basename(fit)),fit_md5=unname(tools::md5sum(fit))))
  back <- modes$read_anchor(file,md5=NULL)
  stopifnot(identical(back$mean,anchor$mean),identical(back$cov,anchor$cov),identical(back$weight,anchor$weight))
  sub <- function(ch) lapply(d,function(m) m[,ch,drop=FALSE])
  half <- function(rows) lapply(d,function(m) m[rows,,drop=FALSE])
  ni <- nrow(d[[1]]);first <- seq_len(ni/2);second <- setdiff(seq_len(ni),first)
  share_mirror <- function(r) modes$anchored_regions(r)$share_mirror
  v <- list()
  addv <- function(check,chain,value) v[[length(v)+1L]] <<- data.frame(check=check,chain=chain,value=value,stringsAsFactors=FALSE)
  r <- modes$assign_anchored(d,anchor)
  addv('in-sample share of draws labelled mirror',1:4,share_mirror(r))
  addv('in-sample share of draws far from both components',1:4,r$atypical$share)
  h <- modes$calibrate_anchor(d,1L,2L);addv('held out: calibrated on chains 1 and 2, share mirror',3:4,share_mirror(modes$assign_anchored(sub(3:4),h)))
  h <- modes$calibrate_anchor(d,3L,4L);addv('held out: calibrated on chains 3 and 4, share mirror',1:2,share_mirror(modes$assign_anchored(sub(1:2),h)))
  h <- modes$calibrate_anchor(half(first),c(1L,3L),c(2L,4L))
  addv('held out: calibrated on first halves, share mirror in second halves',1:4,share_mirror(modes$assign_anchored(half(second),h)))
  addv('theta0 cut at 0.135: share above the cut',1:4,share_mirror(modes$theta0_cut_assignment(d)))
  # Far-from-both share of consecutive blocks of each chain (AMENDMENT-1, R18).
  far <- matrix(0,ni,4);for(ch in 1:4) {one <- modes$assign_anchored(sub(ch),anchor)
    Z <- sapply(anchor$quantities,function(q) d[[q]][,ch])
    d2 <- sapply(1:2,function(k) stats::mahalanobis(Z,anchor$mean[[k]],anchor$cov[[k]]))
    far[,ch] <- pmin(d2[,1],d2[,2])>stats::qchisq(modes$ATYPICAL_LEVEL,length(anchor$quantities))
    stopifnot(isTRUE(all.equal(mean(far[,ch]),one$atypical$share)))}
  for(b in c(1500L,3000L,6000L,12000L)) {
    s <- unlist(lapply(1:4,function(ch) colMeans(matrix(far[seq_len(ni%/%b*b),ch],b))))
    addv(sprintf('largest far-from-both share over consecutive blocks of %d draws (%d blocks)',b,length(s)),NA_integer_,max(s))
  }
  th <- as.vector(d$theta0);dd <- stats::density(th,n=2048L);w <- dd$x>.08 & dd$x<.2
  addv('valley of the pooled theta0 density between 0.08 and 0.2',NA_integer_,dd$x[w][which.min(dd$y[w])])
  addv('squared Mahalanobis distance between the component means (average covariance)',NA_integer_,
    stats::mahalanobis(anchor$mean[[1]],anchor$mean[[2]],(anchor$cov[[1]]+anchor$cov[[2]])/2))
  validation <- do.call(rbind,v);validation$value <- signif(validation$value,6)
  utils::write.csv(validation,file.path(out,'anchor-validation.csv'),row.names=FALSE)
  # Pseudo-chains: nl near-truth and nh mirror blocks of chains {1,3} and {2,4}.
  pseudo <- function(nl,nh,nb) {
    blocks <- function(ch) lapply(d,function(m) do.call(cbind,lapply(ch,function(j) matrix(m[,j],nrow(m)/nb,nb))))
    lo <- blocks(c(1,3));hi <- blocks(c(2,4))
    stats::setNames(lapply(names(d),function(q) cbind(lo[[q]][,seq_len(nl),drop=FALSE],hi[[q]][,seq_len(nh),drop=FALSE])),names(d))
  }
  rows <- list()
  for(cfg in list(c(15,1,8),c(14,2,8),c(8,8,8),c(2,14,8),c(1,15,8),c(7,1,4),c(1,7,4))) {
    x <- pseudo(cfg[1],cfg[2],cfg[3]);truth <- rep(c('near-truth','mirror'),cfg[1:2])
    prim <- modes$assign_anchored(x,anchor);cut <- modes$theta0_cut_assignment(x);refit <- modes$assign_modes(x[modes$MODE_QUANTITIES])
    cmp <- modes$compare_assignments(prim,cut,refit)
    for(m in c('primary','cut','refit')) {
      reg <- cmp[[paste0(m,'_region')]];sm <- cmp[[paste0(m,'_share_mirror')]]
      off <- ifelse(truth=='mirror',1-sm,sm)
      # chains_region_correct: the chain's region is its true one (a chain
      # visiting both, or unknown under R18, is not); chains_majority_correct:
      # its majority component is the true one; chains_unknown: R18 applied.
      rows[[length(rows)+1L]] <- data.frame(near_truth_chains=cfg[1],mirror_chains=cfg[2],draws_per_chain=nrow(x[[1]]),
        method=c(primary='anchored',cut='theta0 cut',refit='refitted mixture')[[m]],
        refit_modes=if(m=='refit') refit$n_modes else NA_integer_,chains=length(truth),
        chains_region_correct=sum(reg==truth),chains_majority_correct=sum(ifelse(sm>.5,'mirror','near-truth')==truth),
        chains_unknown=sum(reg=='unknown'),max_off_region=signif(max(off),4),
        max_share_far=if(m=='primary') signif(max(cmp$primary_share_far),4) else NA_real_,stringsAsFactors=FALSE)
    }
  }
  imbalance <- do.call(rbind,rows)
  utils::write.csv(imbalance,file.path(out,'imbalance-check.csv'),row.names=FALSE)
  print(validation,row.names=FALSE);print(imbalance,row.names=FALSE)
  cat('anchor-classifier.csv md5',unname(tools::md5sum(file)),'\n')
  TRUE
}

# ---- Main --------------------------------------------------------------------------------

status <- tryCatch({
  if(o$mode=='fingerprint') {
    revision <- readLines(file.path(study,'source-revision.txt'))
    write_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
    stopifnot(identical(check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE)),revision))
    cat('Wrote',FINGERPRINT_FILE,'for revision',revision,'\n');0L
  } else if(o$mode=='clone-archive') {
    revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
    load_study_library(study)
    write_fixed_theta0_archive(file.path(repo,CLONE_ARCHIVE_DIR),revision)
    print(check_fixed_theta0_archive(file.path(repo,CLONE_ARCHIVE_DIR)));0L
  } else if(o$mode=='equivalence-fit') {
    equivalence_fit(o$side %||% stop('--side is required'));0L
  } else if(o$mode=='equivalence-compare') {
    if(equivalence_compare()) 0L else 1L
  } else if(o$mode=='equivalence') {
    dir.create(EQ_DIR,showWarnings=FALSE)
    unlink(file.path(EQ_DIR,c('pr11.rds','study.rds')))
    for(side in c('pr11','study')) {
      code <- system2(file.path(R.home('bin'),'Rscript'),shQuote(c(file.path(here,'verify.R'),'--mode=equivalence-fit',
        paste0('--side=',side),paste0('--repo=',repo),paste0('--study=',study),paste0('--archives=',archives),
        paste0('--inputs-root=',inputs_root))),stdout=file.path(EQ_DIR,paste0(side,'.log')),stderr=file.path(EQ_DIR,paste0(side,'.log')))
      cat(side,'fit exit',code,'\n')
      if(code!=0L) stop('Equivalence fit ',side,' failed; see ',file.path(EQ_DIR,paste0(side,'.log')))
    }
    if(equivalence_compare()) {cat('Equivalence PASSES\n');0L} else {cat('Equivalence FAILS\n');1L}
  } else if(o$mode=='anchor-calibration') {
    anchor_calibration();0L
  } else if(o$mode=='mode-calibration') {
    mode_calibration();0L
  } else if(o$mode=='pilot-checks') {
    if(pilot_checks()) {cat('All pilot checks pass\n');0L} else {cat('Pilot checks FAIL\n');1L}
  } else stop('Unknown --mode: ',o$mode)
},error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
quit(save='no',status=status)
