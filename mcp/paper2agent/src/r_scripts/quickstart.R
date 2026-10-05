# File-based dispatcher for pinned occJSDM public APIs. Values arrive as JSON data.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop('Expected activation, request and result paths')
source(normalizePath(args[1],mustWork=TRUE))
suppressPackageStartupMessages(library(occJSDM))
suppressPackageStartupMessages(library(jsonlite))
project <- normalizePath(Sys.getenv('P2A_R_PROJECT'),mustWork=TRUE)
if(!startsWith(normalizePath(find.package('occJSDM')),paste0(project,'/'))) stop('occJSDM resolved outside isolated P2A_R_PROJECT')
if(!startsWith(normalizePath(find.package('jsonlite')),paste0(project,'/'))) stop('jsonlite resolved outside isolated P2A_R_PROJECT')
source_revision <- unique(as.character(stats::na.omit(read.dcf(file.path(find.package('occJSDM'),'DESCRIPTION'),fields='RemoteSha')[,'RemoteSha'])))
if(!identical(source_revision,'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6')) stop('Isolated occJSDM source revision differs from pinned source; restore the recorded runtime')
req <- read_json(args[2],simplifyVector=FALSE)
out <- req$output_dir
warns <- character()

runtime <- function() list(source_revision=source_revision,R=R.version.string,occJSDM_version=as.character(packageVersion('occJSDM')),
  occJSDM_path=find.package('occJSDM'),renv_version=as.character(packageVersion('renv')),
  libPaths=as.list(.libPaths()),thread_env=as.list(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS'))),
  Rscript=normalizePath(file.path(R.home('bin'),'Rscript'),mustWork=FALSE))

native_try <- function(expr) tryCatch(list(ok=TRUE,value=expr),error=function(e)list(ok=FALSE,error=conditionMessage(e)))

# Only metadata, counts and design labels leave this function, never OTU values.
preflight <- function(dat,oc,cc) {
  issues <- character(); notes <- character(); model <- NULL; design <- list()
  add <- function(x) issues <<- c(issues,x)
  done <- function() list(eligible=length(issues)==0L,model=model,issues=as.list(issues),notes=as.list(notes),
    dimensions=if(is.matrix(dat[['OTU',exact=TRUE]]))as.list(dim(dat[['OTU',exact=TRUE]])) else NULL,
    traits_present=!is.null(dat[['traits',exact=TRUE]]),design=design)
  if(!is.list(dat) || is.data.frame(dat)) return(list(eligible=FALSE,model=NULL,issues=list('RDS must contain a native list(info, OTU, traits)'),notes=list()))
  if(is.null(names(dat)) || anyDuplicated(names(dat))) add('Data components must have unique names')
  inf <- dat[['info',exact=TRUE]]; y <- dat[['OTU',exact=TRUE]]; tr <- dat[['traits',exact=TRUE]]
  if(is.null(tr) && any(startsWith(names(dat),'traits'))) add('Use the exact component name traits; partial trait aliases can bypass checks and are not accepted in the pilot')
  if(!is.data.frame(inf)) add('info must be a data.frame with one row per observation')
  if(!is.matrix(y) || !is.numeric(y)) add('OTU must be a numeric matrix with rows matching info and named species columns')
  if(length(issues)) return(done())
  if(nrow(y)!=nrow(inf)) add('OTU and info must have the same number of rows')
  if(nrow(y)<1L || ncol(y)<1L) add('OTU must contain observations and species')
  if(anyDuplicated(names(inf)) || any(!nzchar(names(inf)))) add('info columns must have unique nonempty names')
  sp <- colnames(y)
  if(is.null(sp) || anyNA(sp) || any(!nzchar(trimws(sp))) || anyDuplicated(sp)) add('OTU species column names must be present, unique and nonempty')
  if(!is.null(rownames(y)) && .row_names_info(inf,1L)>0L && !identical(rownames(inf),rownames(y))) add('Supplied row identifiers in info and OTU do not match in order; pair the observations before saving RDS')
  if(any(is.nan(y)) || any(is.infinite(y)) || any(y<0,na.rm=TRUE)) add('OTU contains nonfinite or negative reads; only NA observation values can be supported')
  for(id in intersect(c('Site','Sample','Primer'),names(inf))) {
    x <- inf[[id]]
    if(!is.atomic(x) || anyNA(x) || (is.numeric(x) && any(!is.finite(x))) || any(!nzchar(trimws(as.character(x))))) add(paste('Invalid or missing identifier in',id))
  }
  selected <- unique(c(oc,cc))
  missing <- setdiff(selected,names(inf))
  if(length(missing)) add(paste('Selected covariate names not in info:',paste(missing,collapse=', ')))
  for(nm in intersect(selected,names(inf))) {
    x <- inf[[nm]]
    if(anyNA(x)) add(paste('NA in selected covariate',nm))
    if(is.numeric(x) && any(!is.finite(x))) add(paste('Nonfinite values in selected covariate',nm))
  }
  if(!is.null(tr)) {
    if(!(is.matrix(tr) || is.data.frame(tr))) add('traits must be a matrix or data.frame with species row names')
    else {
      ids <- rownames(tr)
      if(is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) add('traits row names must be present, unique species names')
      else if(!all(sp %in% ids)) add('Species names in OTU not present in traits row names')
      if(anyNA(tr)) add('NA in traits covariates')
      num <- if(is.matrix(tr) && is.numeric(tr)) list(tr) else if(is.data.frame(tr)) tr[vapply(tr,is.numeric,logical(1))] else list()
      if(any(vapply(num,function(x)any(!is.finite(x)),logical(1)))) add('Nonfinite numeric traits')
    }
  } else notes <- c(notes,'Traits absent: native no-traits route retained.')
  if(length(issues)) return(done())
  inferred <- native_try(occJSDM:::inferDataModel(dat))
  if(!inferred$ok) {add(paste('Native model inference:',inferred$error));return(done())}
  model <- inferred$value
  if(!(model %in% c('binary','occupancy','two_stage'))) add(paste('Pilot supports non-spatial binary, occupancy and two_stage models; native inferred',model))
  if(anyNA(y) && model!='two_stage') add('NAs are allowed only in the native two-stage model')
  if(model=='two_stage' && !('Primer' %in% names(inf))) add('Native inference selected two_stage but Primer is missing. Reused local Sample IDs across sites can cause this inference; use globally unique IDs for one-stage occupancy or supply actual PCR/primer metadata for a two-stage survey.')
  if(model=='binary' && !is.null(inf$Site)) {
    ord <- dplyr::arrange(inf,Site)
    if(!identical(inf$Site,ord$Site)) add('Pinned binary preparation sorts covariates by unique Site but retains OTU input order. Explicit Site must already be ordered; reorder paired info and OTU yourself before saving, or omit Site when each row is a site.')
  }
  if(length(issues)) return(done())
  if(is.null(inf$Site)) inf$Site <- seq_len(nrow(inf))
  if(model %in% c('occupancy','two_stage') && is.null(inf$Sample)) inf$Sample <- seq_len(nrow(inf))
  n <- length(unique(inf$Site))
  invariant <- function(cols,groups,label) {
    grouped <- dplyr::group_by(inf,dplyr::across(dplyr::all_of(groups)))
    indices <- dplyr::group_rows(grouped)
    for(nm in cols) if(any(vapply(indices,function(i)length(unique(inf[[nm]][i]))>1L,logical(1)))) add(paste('Selected covariate',nm,'conflicts within',label,'; supply one consistent value per unit'))
  }
  invariant(oc,'Site','Site')
  if(model %in% c('occupancy','two_stage')) invariant(cc,c('Site','Sample'),'Site/Sample')
  else if(length(cc)) notes <- c(notes,'Binary fits skip collection detection; selected collection covariates are checked for existence but unused by the native sampler.')
  if(length(issues)) return(done())
  prepared <- native_try({
    xp <- occJSDM:::process_covariates(inf,oc,'Site',n,remove_intercept=TRUE,spline_vars=FALSE)$df
    if(nrow(xp)!=n || any(!is.finite(xp))) stop('Occupancy covariate preparation produced missing/nonfinite design rows; check constant covariates')
    xt <- NULL; ns <- NULL
    if(model %in% c('occupancy','two_stage')) {
      samples <- dplyr::ungroup(dplyr::slice(dplyr::group_by(inf,Site,Sample),1))
      ns <- nrow(samples)
      xt <- occJSDM:::process_covariates(samples,cc,NULL,ns,remove_intercept=FALSE,spline_vars=FALSE)$df
      if(nrow(xt)!=ns || any(!is.finite(xt))) stop('Collection covariate preparation produced missing/nonfinite design rows; check constant covariates')
    }
    traits <- NULL
    if(!is.null(tr)) {
      traits <- occJSDM:::create_covariates_matrix(tr[match(sp,rownames(tr)),,drop=FALSE],spline_vars=FALSE,remove_intercept=TRUE)$X
      if(nrow(traits)!=length(sp) || any(!is.finite(traits))) stop('Trait preparation produced missing/nonfinite rows; check constant traits')
    }
    list(occupancy_columns=as.list(colnames(xp)),collection_columns=as.list(colnames(xt)),trait_columns=as.list(colnames(traits)),
      n_sites=n,n_samples=ns,n_primers=if(model=='two_stage')length(unique(inf$Primer))else NULL,
      missing_observations=sum(is.na(y)),species_count=ncol(y),observations=nrow(y))
  })
  if(prepared$ok) design <- prepared$value else add(paste('Native covariate preparation:',prepared$error))
  done()
}

validate_fit <- function(fit,settings=NULL) {
  if(!is.list(fit) || !is.list(fit$infos) || !is.list(fit$results_output)) stop('RDS does not contain a native completed fit')
  b <- fit$results_output$jsdm_output$B0_output
  if(length(dim(b))!=3L || dim(b)[1]!=length(fit$infos$speciesNames) || any(!is.finite(b))) stop('Invalid occupancy intercept output in native fit')
  if(!is.null(settings) && !identical(as.integer(dim(b)[2:3]),as.integer(c(settings$niter,settings$nchain)))) stop('Native fit retained draw or chain dimensions differ from request')
  TRUE
}

screen <- function(tab) {
  rh <- tab$rhat; es <- tab$ess
  list(total_rows=nrow(tab),rhat_flagged=sum(is.finite(rh) & rh>1.01),ess_flagged=sum(is.finite(es) & es<400),
    flagged_rows=sum((is.finite(rh) & rh>1.01)|(is.finite(es)&es<400)),
    unavailable_rhat=sum(!is.finite(rh)),unavailable_ess=sum(!is.finite(es)),
    unavailable_rows=sum(!is.finite(rh)|!is.finite(es)),thresholds=list(rhat=1.01,ess=400),
    blocks=as.list(unique(tab$param)),scope='Selected occupancy/detection coefficient blocks; excludes latent variables and checks of model assumptions',
    unavailable_reason='Rhat is unavailable for constant chains or insufficient chains; unavailable or nonfinite statistics are inconclusive.')
}

run <- function() {
  if(req$operation=='validate_data') {
    d <- native_try(readRDS(req$data_path))
    if(!d$ok) return(list(status='ok',eligible=FALSE,model=NULL,issues=list('RDS could not be read as a native input object'),notes=list()))
    return(c(list(status='ok'),preflight(d$value,unlist(req$occCovariates,use.names=FALSE),unlist(req$collCovariates,use.names=FALSE))))
  }
  if(req$operation=='fit_model') {
    d <- readRDS(req$data_path);o <- req$options
    oc <- as.character(unlist(o$occCovariates,use.names=FALSE));cc <- as.character(unlist(o$collCovariates,use.names=FALSE))
    checked <- preflight(d,oc,cc)
    if(!checked$eligible) stop(paste(unlist(checked$issues),collapse='; '))
    set.seed(as.integer(req$seed))
    fit <- occJSDM::runOccJSDM(data=d,listParams=o$listParams,threshold=o$threshold,
      occCovariates=oc,collCovariates=cc,spatCovariates=character(),MCMCparams=o$MCMCparams,
      summarisedLatentPresences=TRUE)
    validate_fit(fit,o$MCMCparams)
    pending <- file.path(out,'fit.rds.pending')
    saveRDS(fit,pending)
    valid <- validate_fit(readRDS(pending),o$MCMCparams)
    notes <- list('Beta research software: faithful execution does not establish unbiased estimates or calibrated intervals.')
    if(o$MCMCparams$niter<400 || o$MCMCparams$nburn<400) notes <- c(notes,'Short fits are transport/teaching checks and are unsuitable for scientific interpretation; inspect convergence before extending sampling.')
    return(list(status='ok',model=fit$infos$model,readback_valid=valid,
      effective_n_factors=fit$infos$n_factors,kept_draws_per_chain=o$MCMCparams$niter,total_iterations_per_chain=o$MCMCparams$nburn+o$MCMCparams$niter*o$MCMCparams$nthin,
      dimensions=checked$dimensions,design=checked$design,runtime=runtime(),qualifications=notes))
  }
  if(!(req$operation %in% c('diagnostics','summarise_fit'))) stop('Unknown fixed dispatcher operation')
  fit <- readRDS(req$fit_path);validate_fit(fit)
  tab <- occJSDM::returnConvergenceDiagnostics(fit)
  write.csv(tab,file.path(out,'diagnostics.csv'),row.names=FALSE,na='NA')
  screening <- screen(tab)
  if(req$operation=='diagnostics') return(list(status='ok',screening=screening,
    rows=utils::head(tab,as.integer(req$max_rows)),total_rows=nrow(tab)))
  rates <- occJSDM::returnOccupancyRates(fit)
  if(!is.matrix(rates) || ncol(rates)!=length(fit$infos$speciesNames) || !identical(colnames(rates),fit$infos$speciesNames)) stop('Unexpected native occupancy draws orientation')
  summary <- data.frame(species=colnames(rates),mean=apply(rates,2,mean),
    median=apply(rates,2,stats::quantile,probs=.5),q2.5=apply(rates,2,stats::quantile,probs=.025),
    q97.5=apply(rates,2,stats::quantile,probs=.975),row.names=NULL)
  write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE,na='NA')
  list(status='ok',estimand='Baseline occupancy: inverse-logit species intercept, pooled across chains and retained draws',
    rows=utils::head(summary,as.integer(req$max_rows)),total_rows=nrow(summary),screening=screening,
    qualifications=list('Baseline occupancy is not occupancy at a particular site or observed detections.',
      'Beta research software: wrapper agreement does not establish unbiased estimates or calibrated intervals.',
      'Flagged or unavailable diagnostics require investigation before interpretation.'))
}

result <- tryCatch(withCallingHandlers(run(),warning=function(w){
  warns <<- c(warns,conditionMessage(w));invokeRestart('muffleWarning')
}),error=function(e)list(status='error',error=conditionMessage(e)))
result$warnings <- as.list(utils::head(unique(warns),20L))
result$warning_count <- length(warns)
# jsonlite serializes NA/NaN/Inf as null under na='null', including diagnostic previews.
write_json(result,args[3],auto_unbox=TRUE,pretty=TRUE,na='null',null='null',digits=17,dataframe='rows')
if(result$status=='error') quit(status=1L)
