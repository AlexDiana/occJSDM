args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1], mustWork=TRUE)
Sys.setenv(P2A_R_PROJECT=file.path(root,"r-runtime"))
source(file.path(root,"r-runtime","activate.R"))
library(occJSDM)
library(jsonlite)
setTimeLimit(elapsed=600, transient=FALSE)
out <- file.path(root,"notebooks","quickstart")
repo <- file.path(root,"repo","occJSDM")
source(file.path(repo,"tests","testthat","helper-fixtures.R"))
# Evaluate only the two native fixture definitions; never execute test_that blocks.
align_expr <- parse(file.path(repo,"tests","testthat","test-collection-alignment.R"))
for (expr in align_expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% c("alignment_samples","alignment_data")) eval(expr)
}
cases <- list(); probes <- list(); checks <- list(); fits <- list()
runtime <- list(R=R.version.string, Rscript="/usr/local/bin/Rscript",
  libPaths=.libPaths(), occJSDM_version=as.character(packageVersion("occJSDM")),
  occJSDM_path=find.package("occJSDM"), renv_version=as.character(packageVersion("renv")),
  session=capture.output(sessionInfo()),
  thread_env=as.list(Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"))))
capture_native <- function(expr) {
  warnings <- character(); messages <- character()
  ans <- tryCatch(withCallingHandlers(expr,
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")},
    message=function(m) {messages <<- c(messages,conditionMessage(m)); invokeRestart("muffleMessage")}),
    error=function(e) structure(list(error=conditionMessage(e),class=class(e)),class="native_failure"))
  if (inherits(ans,"native_failure")) list(status="error",error=ans$error,
    error_class=ans$class,warnings=warnings,messages=messages)
  else list(status="ok",value=ans,warnings=warnings,messages=messages)
}
check <- function(id, actual, expected, tol=0) {
  comp <- all.equal(actual,expected,tolerance=tol,check.attributes=TRUE)
  checks[[id]] <<- list(pass=isTRUE(comp),tolerance=tol,
    comparator="R all.equal(check.attributes=TRUE)",detail=if(isTRUE(comp)) NULL else comp)
}
native_prepare <- function(dat, occ=character(), coll=character()) {
  model <- occJSDM:::inferDataModel(dat)
  inf <- dat$info
  if(is.null(inf$Site)) inf$Site <- seq_len(nrow(inf))
  n <- length(unique(inf$Site))
  xp <- occJSDM:::process_covariates(inf,occ,"Site",n,remove_intercept=TRUE,spline_vars=FALSE)
  xt <- NULL
  if(model %in% c("occupancy","two_stage")) {
    if(is.null(inf$Sample)) inf$Sample <- seq_len(nrow(inf))
    sam <- dplyr::ungroup(dplyr::slice(dplyr::group_by(inf,Site,Sample),1))
    xt <- occJSDM:::process_covariates(sam,coll,NULL,nrow(sam),remove_intercept=FALSE,spline_vars=FALSE)
  }
  list(model=model,X_psi=xp,X_theta=xt)
}
defaults <- list(nchain=2L,nburn=20L,niter=20L,nthin=1L)
fit_args <- function(dat,occ,coll,mcmc=defaults) list(data=dat,
  listParams=list(n_factors=2L),threshold=1,
  occCovariates=occ,collCovariates=coll,spatCovariates=character(),
  MCMCparams=mcmc,summarisedLatentPresences=TRUE,listPriors=list())
post <- function(id,fit) {
  dr <- capture_native(returnConvergenceDiagnostics(fit))
  rr <- capture_native(returnOccupancyRates(fit))
  ans <- list(diagnostics=dr[setdiff(names(dr),"value")],rates=rr[setdiff(names(rr),"value")])
  if(dr$status=="ok") {
    write.csv(dr$value,file.path(out,"reference",paste0(id,"_diagnostics.csv")),row.names=FALSE,na="NA")
    saveRDS(dr$value,file.path(out,"reference",paste0(id,"_diagnostics.rds")))
    ans$diagnostics$rows <- nrow(dr$value); ans$diagnostics$columns <- names(dr$value)
    ans$diagnostics$blocks <- unique(dr$value$param)
    ans$diagnostics$nonfinite <- lapply(dr$value[vapply(dr$value,is.numeric,logical(1))],function(x) sum(!is.finite(x)))
  }
  if(rr$status=="ok") {
    rates <- rr$value
    saveRDS(rates,file.path(out,"reference",paste0(id,"_occupancy_draws.rds")))
    summary <- data.frame(species=colnames(rates),
      mean=apply(rates,2,mean),
      median=apply(rates,2,stats::quantile,probs=.5),
      q2.5=apply(rates,2,stats::quantile,probs=.025),
      q97.5=apply(rates,2,stats::quantile,probs=.975),row.names=NULL)
    write.csv(summary,file.path(out,"reference",paste0(id,"_occupancy_summary.csv")),row.names=FALSE)
    ans$rates$dim <- dim(rates); ans$rates$species <- colnames(rates)
    ans$rates$orientation <- "pooled iteration/chain draws x species"
  }
  ans
}
run_case <- function(id,dat,occ=character(),coll=character(),seed=1702L,mcmc=defaults,replay=FALSE) {
  cat("CASE",id,"\n")
  saveRDS(dat,file.path(out,"data",paste0(id,"_data.rds")))
  argv <- fit_args(dat,occ,coll,mcmc)
  saveRDS(list(seed=seed,args=argv),file.path(out,"data",paste0(id,"_call.rds")))
  prep <- capture_native(native_prepare(dat,occ,coll))
  if(prep$status=="ok") saveRDS(prep$value,file.path(out,"reference",paste0(id,"_preparation.rds")))
  set.seed(seed)
  r <- capture_native(do.call(runOccJSDM,argv))
  rec <- list(seed=seed,settings=mcmc,occCovariates=occ,collCovariates=coll,
    input_dim=dim(dat$OTU),input_species=colnames(dat$OTU),info_columns=names(dat$info),
    preparation=prep[setdiff(names(prep),"value")],
    fit=r[setdiff(names(r),"value")])
  if(r$status=="ok") {
    fit <- r$value; fits[[id]] <<- fit
    saveRDS(fit,file.path(out,"reference",paste0(id,"_fit.rds")))
    rec$inferred_model <- fit$infos$model; rec$species <- fit$infos$speciesNames
    rec$B0_dim <- dim(fit$results_output$jsdm_output$B0_output)
    rec$X_psi_dim <- dim(fit$X_psi); rec$X_psi_columns <- colnames(fit$X_psi)
    rec$X_theta_dim <- dim(fit$X_theta); rec$X_theta_columns <- colnames(fit$X_theta)
    rec$Tr_dim <- dim(fit$Tr); rec$Tr_columns <- colnames(fit$Tr)
    rec$site_names <- fit$infos$siteNames
    rec$postprocessing <- post(id,fit)
    check(paste0(id,"_kept_draws"),unname(rec$B0_dim),c(ncol(dat$OTU),mcmc$niter,mcmc$nchain))
    if(replay) {
      saved <- readRDS(file.path(out,"data",paste0(id,"_call.rds")))
      set.seed(saved$seed)
      again <- capture_native(do.call(runOccJSDM,saved$args))
      rec$replay <- again[setdiff(names(again),"value")]
      if(again$status=="ok") {
        saveRDS(again$value,file.path(out,"reference",paste0(id,"_replay_fit.rds")))
        for(key in c("jsdm_output","beta_theta_output","p_output","q_output","theta0_output"))
          check(paste0(id,"_replay_",key),again$value$results_output[[key]],fit$results_output[[key]])
        check(paste0(id,"_replay_X_psi"),again$value$X_psi,fit$X_psi)
        check(paste0(id,"_replay_X_theta"),again$value$X_theta,fit$X_theta)
        check(paste0(id,"_replay_species"),again$value$infos$speciesNames,fit$infos$speciesNames)
        dg2 <- capture_native(returnConvergenceDiagnostics(again$value))
        if(dg2$status=="ok" && rec$postprocessing$diagnostics$status=="ok")
          check(paste0(id,"_replay_diagnostics"),dg2$value,readRDS(file.path(out,"reference",paste0(id,"_diagnostics.rds"))))
        rate2 <- capture_native(returnOccupancyRates(again$value))
        if(rate2$status=="ok")
          check(paste0(id,"_replay_rates"),rate2$value,readRDS(file.path(out,"reference",paste0(id,"_occupancy_draws.rds"))))
      }
    }
  }
  cases[[id]] <<- rec
  invisible(r)
}
fixtures <- list()
for(model in c("binary","occupancy","two_stage")) {
  sim <- simulate_fixture(model=model,useSpatField=FALSE,seed=1701L)
  saveRDS(sim,file.path(out,"data",paste0(model,"_simulation.rds")))
  fixtures[[model]] <- sim$data_list
  occ <- fixture_occ_covariates()
  coll <- if(model=="binary") character() else grep("^X_theta",names(sim$data_list$info),value=TRUE)
  run_case(model,sim$data_list,occ,coll,replay=TRUE)
}
dat <- fixtures$two_stage; occ <- fixture_occ_covariates(); coll <- grep("^X_theta",names(dat$info),value=TRUE)
nt <- dat; nt$traits <- NULL
run_case("no_traits",nt,occ,coll)
ren <- dat; names(ren$info)[match(occ[1],names(ren$info))] <- "habitat"
names(ren$info)[match(coll[1],names(ren$info))] <- "effort"
run_case("renamed_single_covariates",ren,"habitat","effort")
per <- dat; per$traits <- per$traits[rev(seq_len(nrow(per$traits))),,drop=FALSE]
run_case("permuted_trait_rows",per,occ,coll)
if(!is.null(fits$permuted_trait_rows)) check("traits_reordered_by_species",fits$permuted_trait_rows$Tr,fits$two_stage$Tr)
catdat <- dat
set.seed(1701L)
catdat$traits <- data.frame(mass=rnorm(ncol(dat$OTU)),
  diet=factor(rep(c("carnivore","herbivore"),length.out=ncol(dat$OTU))),row.names=colnames(dat$OTU))
run_case("categorical_traits",catdat,occ,coll)
run_case("changed_seed_thinning2",dat,occ,coll,seed=1703L,
 mcmc=list(nchain=2L,nburn=20L,niter=10L,nthin=2L),replay=TRUE)
if(!is.null(fits$changed_seed_thinning2)) checks$changed_seed_changes_draws <- list(
 pass=!identical(fits$two_stage$results_output$jsdm_output$B0_output,
 fits$changed_seed_thinning2$results_output$jsdm_output$B0_output),
 comparator="different seed and thinning schedule; nonidentity only, no seed-only isolation",tolerance=NULL)
run_case("changed_seed_same_schedule",dat,occ,coll,seed=1703L)
if(!is.null(fits$changed_seed_same_schedule)) checks$seed_only_changes_draws <- list(
 pass=!identical(fits$two_stage$results_output$jsdm_output$B0_output,fits$changed_seed_same_schedule$results_output$jsdm_output$B0_output),
 comparator="R identical; seed changed, fixed input/schedule",tolerance=NULL)
for(ids in c("global","local","character")) {
  sam <- alignment_samples(ids)
  for(shuffled in c(FALSE,TRUE)) {
    id <- paste0("alignment_",ids,if(shuffled) "_shuffled" else "_ordered")
    run_case(id,alignment_data(sam,shuffled=shuffled),coll=c("effort","collection_type"))
    fit <- fits[[id]]
    if(!is.null(fit)) {
      expected_effort <- (sam$effort-mean(sam$effort))/sd(sam$effort)
      check(paste0(id,"_effort"),unname(fit$X_theta[,"effort"]),expected_effort,tol=1e-14)
      check(paste0(id,"_category"),unname(fit$X_theta[,"collection_typewet"]),as.numeric(sam$collection_type=="wet"))
      check(paste0(id,"_M"),unname(fit$infos$M),c(2,2,2))
      check(paste0(id,"_K"),unname(fit$infos$K),c(2,1,3,1,2,2,1,2,1))
      expected_sample <- rep(c(1,1,2,3,3,4,5,6,6),c(2,1,3,1,2,2,1,2,1))
      check(paste0(id,"_idx_w_k"),fit$infos$list_idx$idx_w_k,expected_sample)
      check(paste0(id,"_info_pair"),fit$infos$data_info$effort,sam$effort[expected_sample])
      check(paste0(id,"_otu_pair"),unname(fit$infos$OTU[,"sp1"]),c(5,0,8,2,0,3)[expected_sample])
    }
  }
}
run_case("alignment_identifier_covariates",alignment_data(alignment_samples("local"),shuffled=TRUE),coll=c("Site","Sample"))
run_case("alignment_occupancy_shuffled",alignment_data(alignment_samples(),"occupancy",TRUE),coll=c("effort","collection_type"))
# Native behavior probes: no wrapper validation, no repair and no substitutions.
probe_fit <- function(id,d,oc=occ,cc=coll,mm=defaults) {
  ans <- run_case(paste0("boundary_",id),d,oc,cc,mcmc=mm)
  probes[[id]] <<- cases[[paste0("boundary_",id)]]
}
binary <- fixtures$binary; occupancy <- fixtures$occupancy
u <- binary; u$info$Site <- rev(seq_len(nrow(u$info)))
probe_fit("binary_unsorted_uniqueSite",u,fixture_occ_covariates(),character())
if(!is.null(fits$boundary_binary_unsorted_uniqueSite)) {
  xp <- fits$boundary_binary_unsorted_uniqueSite$X_psi
  ordered <- order(u$info$Site)
  x <- u$info[[occ[1]]]
  expected <- (x-mean(x))/sd(x)
  probes$binary_unsorted_uniqueSite$X_psi_follows_sorted_site <- isTRUE(all.equal(unname(xp[,occ[1]]),expected[ordered],tolerance=1e-14))
  probes$binary_unsorted_uniqueSite$OTU_follows_input <- identical(fits$boundary_binary_unsorted_uniqueSite$infos$OTU,u$OTU)
  probes$binary_unsorted_uniqueSite$X_psi_follows_input_rows <- isTRUE(all.equal(unname(xp[,occ[1]]),expected,tolerance=1e-14))
}
probe_fit("occupancy_localSampleIDs",alignment_data(alignment_samples("local"),"occupancy"),character(),c("effort","collection_type"))
for(model in c("binary","occupancy","two_stage")) {
  d <- fixtures[[model]]; d$OTU[1,1] <- NA_real_
  probe_fit(paste0(model,"_NA_OTU"),d,occ,if(model=="binary") character() else coll)
}
probe_fit("niter1",dat,mm=list(nchain=2L,nburn=20L,niter=1L,nthin=1L))
for(val in c("NA","Inf","NaN")) {
  d <- dat; d$info[[occ[1]]][1] <- switch(val,"NA"=NA_real_,"Inf"=Inf,"NaN"=NaN)
  probe_fit(paste0(val,"_occupancy_covariate"),d)
}
d <- binary; d$OTU[1,1] <- 5; probe_fit("no_stage_count_response",d,occ,character())
d <- dat; d$OTU <- d$OTU[-1,,drop=FALSE]; probe_fit("wrong_dimensions",d)
probe_fit("missing_covariate",dat,"absent_habitat",coll)
d <- dat; rownames(d$traits) <- paste0("unknown",seq_len(nrow(d$traits))); probe_fit("missing_trait_species",d)
d <- dat; d$info$Primer <- NULL; probe_fit("missing_Primer",d)
d <- dat; d$info$Site[1] <- NA; probe_fit("NA_site_ID",d)
d <- dat; colnames(d$OTU)[2] <- colnames(d$OTU)[1]; probe_fit("duplicated_species",d)
d <- dat; d$info[[occ[1]]][1] <- d$info[[occ[1]]][1]+10
probe_fit("within_site_covariate_conflict",d)
d <- dat; d$info[[coll[1]]][1] <- d$info[[coll[1]]][1]+10
probe_fit("within_sample_covariate_conflict",d)
d <- dat; d$OTU[1,1] <- Inf
probe_fit("infinite_OTU",d)
data("sampledata",package="occJSDM")
data("sampleresults",package="occJSDM")
saveRDS(sampledata,file.path(out,"data","bundled_sampledata.rds"))
saveRDS(sampleresults,file.path(out,"reference","bundled_sampleresults_fit.rds"))
bundled_prepare <- capture_native(native_prepare(sampledata,
 grep("^X_psi",names(sampledata$info),value=TRUE),
 grep("^X_theta",names(sampledata$info),value=TRUE)))
if(bundled_prepare$status=="ok") saveRDS(bundled_prepare$value,file.path(out,"reference","bundled_sampledata_preparation.rds"))
bundled <- list(preparation=bundled_prepare[setdiff(names(bundled_prepare),"value")],
 input_dim=dim(sampledata$OTU), input_species=colnames(sampledata$OTU),
 fit_model=sampleresults$infos$model, B0_dim=dim(sampleresults$results_output$jsdm_output$B0_output),
 postprocessing=post("bundled_sampleresults",sampleresults),
 limitation="Bundled stored fit is spatial; fresh pilot fits are nonspatial. This is a direct postprocessing reference, not a fresh nonspatial fitting reference.")
result <- list(runtime=runtime,cases=cases,boundaries=probes,checks=checks,bundled=bundled,
 counts=list(cases=length(cases),fit_ok=sum(vapply(cases,function(x)x$fit$status=="ok",logical(1))),
  fit_errors=sum(vapply(cases,function(x)x$fit$status=="error",logical(1))),
  checks=length(checks),passed=sum(vapply(checks,function(x)isTRUE(x$pass),logical(1)))))
saveRDS(result,file.path(out,"reference","native_execution_results.rds"))
write_json(result,file.path(out,"reference","native_execution_results.json"),pretty=TRUE,auto_unbox=TRUE,na="null",null="null",digits=17)
cat("FINAL_COUNTS",toJSON(result$counts,auto_unbox=TRUE),"\n")
