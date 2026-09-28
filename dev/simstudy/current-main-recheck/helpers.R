# Research-only paired rerun. Archived files are read, never rewritten.
extract_functions <- function(path, names, env) {
  expressions <- as.list(parse(path))
  for (name in names) {
    hits <- Filter(function(e) is.call(e) && identical(e[[1]], as.name('<-')) &&
      identical(e[[2]], as.name(name)), expressions)
    stopifnot(length(hits) == 1L)
    eval(hits[[1]], envir = env)
  }
}

load_scoring <- function(repo) {
  design <- new.env(parent = globalenv())
  design$`%||%` <- function(x,y) if (is.null(x)) y else x
  extract_functions(file.path(repo,'dev/simstudy/nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R'),
    c('trace_summary','score_fit'), design)
  sys.source(file.path(repo,'dev/simstudy/nonspatial-design-recheck/design_helpers.R'), envir=design)
  # Four-chain checks must compare every chain, not just the first two.
  design$trace_summary <- function(x) {
    x <- as.matrix(x)
    c(rhat=posterior::rhat(x),ess_bulk=posterior::ess_bulk(x),
      ess_mean=posterior::ess_mean(x),mcse=posterior::mcse_mean(x),
      chain_gap=diff(range(colMeans(x))))
  }
  design$score <- design$make_design_scorer()
  original_body <- body(design$score)
  body(design$score) <- bquote({
    result <- .(original_body)
    diagnostic_block <- function(arr,metric) {
      dd <- dim(arr);np <- prod(head(dd,-2L));flat <- array(arr,c(np,ni,nc))
      do.call(rbind,lapply(seq_len(np),function(k)
        data.frame(metric=metric,element=k,as.list(trace_summary(flat[k,,])))))
    }
    result$additional_diagnostics <- rbind(
      diagnostic_block(jo$B_output,'environment_slope'),
      diagnostic_block(psi_draws[base_idx,,,drop=FALSE],'original_probability'),
      diagnostic_block(jo$sigmah_output,'sigma_h'))
    result$estimate <- reconstructed
    result$truth <- true_psi
    result
  })
  jsdm <- new.env(parent = globalenv())
  sys.source(file.path(repo,'dev/simstudy/jsdm-sample-size-recheck/helpers.R'), envir=jsdm)
  list(design=design, jsdm=jsdm)
}

make_jobs <- function(archives) {
  bias <- file.path(archives,'nonspatial-bias-recheck-20260914')
  design <- file.path(archives,'nonspatial-design-recheck-20260919')
  jsdm <- file.path(archives,'jsdm-sample-size-20260919')
  jobs <- list()
  add <- function(key,family,arm,scenario,replicate,input,old) {
    stopifnot(file.exists(input),file.exists(old))
    original <- readRDS(old)
    jobs[[key]] <<- list(key=key,family=family,arm=arm,scenario=scenario,
      replicate=replicate,input_file=normalizePath(input),
      input_md5=unname(tools::md5sum(input)),old_initial_file=normalizePath(old),
      priors=if(family=='design') original$job$priors else list())
  }
  for (r in 1:10) for (n in c(100L,300L,1000L)) {
    k <- sprintf('n%04d-%02d',n,r)
    add(paste0('jsdm-',k),'jsdm',paste0('sites',n),'perfect',r,
      file.path(jsdm,'inputs',paste0(k,'.rds')),file.path(jsdm,'fits',paste0(k,'-result.rds')))
  }
  for (sc in c('qnear_K6','qfar_K6')) for (r in 1:10) {
    k <- sprintf('%s-q20-%02d',sc,r)
    add(paste0('design-',k),'design','baseline',sc,r,
      file.path(bias,'operational-k6',sprintf('%s-data-%02d.rds',sc,r)),
      file.path(bias,'operational-k6',paste0(k,'-result.rds')))
    for (arm in c('field4','sites300','knownU')) {
      k <- sprintf('%s-%s-%02d',sc,arm,r)
      add(paste0('design-',k),'design',arm,sc,r,
        file.path(design,'inputs',paste0(k,'.rds')),file.path(design,'fits',paste0(k,'-result.rds')))
    }
  }
  stopifnot(length(jobs)==110L,!anyDuplicated(names(jobs)))
  jobs
}

check_current_fit <- function(fit,input,job) {
  dat <- if(job$family=='jsdm') input$data else input$sim$data_list
  raw <- dat$traits[match(fit$infos$speciesNames,rownames(dat$traits)),,drop=FALSE]
  meta <- fit$infos$list_Tr_mat
  stopifnot(fit$infos$ps==0L,identical(fit$infos$speciesNames,colnames(dat$OTU)),
    !is.null(meta),max(abs(meta$mean_df-colMeans(raw)))<1e-12,
    max(abs(meta$sd_df-apply(raw,2,sd)))<1e-12,
    max(abs(unname(fit$Tr)-unname(scale(raw))))<1e-12)
  invisible(TRUE)
}

score_current_fit <- function(fit,input,job,scoring) {
  if(job$family=='jsdm') {
    result <- scoring$jsdm$score_fit(fit,input)
    result$additional_diagnostics <- data.frame(metric='sigma_h',element=1L,
      as.list(scoring$jsdm$trace_stats(fit$results_output$jsdm_output$sigmah_output)))
    return(result)
  }
  if(job$arm=='baseline') input$design <- list(original_sites=100L,
    original_samples=seq_len(nrow(fit$X_theta)))
  scoring$design$score(fit,input,job$priors)
}

safe_max <- function(x) if(!length(x) || all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)
