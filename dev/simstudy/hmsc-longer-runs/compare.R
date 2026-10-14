#!/usr/bin/env Rscript
# Compare every refit with its saved longer attempt and recompute Lesson 6's
# Hmsc cells with the refits substituted (README.md in this directory).
#
#   Rscript compare.R --repo=REPO --root=EXTENSION_ARCHIVE --study=ARCHIVE [--out=DIR]
#
# This is the first script of the check that opens truth/. For each selected
# job it reads the saved attempt-2 draws, score and diagnostics from the
# extension archive and the attempt-3 draws, diagnostics and result from
# ARCHIVE, and computes with the extension's own code: the 300-test-site
# prediction (checked_posterior_mean, error_summary), the marginal-probability
# intervals at the five fixed environments (probability_intervals, truth from
# marginal_truth), the probit coefficient posterior means and 95% intervals on
# the fitted scale, and the Hmsc trait rows where traits were supplied. The
# saved side is checked against the calibration bundle before any comparison.
# Then the bundle's Hmsc rows and predictions are copied with the selected
# jobs' values replaced by the refits' and the lesson's grouped summaries
# recomputed. Writes, under results/ (or --out): comparison-by-job.csv,
# lesson6-probability-cells.csv, lesson6-prediction-cells.csv, summary.csv and
# decision.txt.
here <- local({
  a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  f <- NULL
  for(i in rev(seq_len(sys.nframe()))) {x <- sys.frame(i)$ofile;if(!is.null(x)) {f <- x;break}}
  dirname(normalizePath(f %||% sub('^--file=','',a[1])))
})
`%||%` <- function(x,y) if(is.null(x)) y else x
source(file.path(here,'common.R'))

# A Lesson 6 cell counts as changed when its coverage moves by at least this
# many percentage points, or its bias or RMSE by at least this many points of
# probability (README.md, "Decision rules").
FLAG_COVERAGE_POINTS <- 1
FLAG_ERROR_POINTS <- .5
GRID_LABELS <- c('mean environment','gradient 1 low','gradient 1 high','gradient 2 low','gradient 2 high')
CELL_KEYS <- c('scenario','n_sites','n_species','response','use_traits')

coefficient_summary <- function(beta) list(mean=apply(beta,1:2,mean),sd=apply(beta,1:2,stats::sd),
  lower=apply(beta,1:2,stats::quantile,probs=.025),upper=apply(beta,1:2,stats::quantile,probs=.975))

hmsc_trait_rows <- function(p,input) {
  out <- NULL
  for(t in 1:2) for(e in 1:2) {
    draws <- as.vector(p$trait[e+1,t+1,,])/input$spread[e]
    out <- rbind(out,c(estimate=mean(draws),stats::quantile(draws,c(.025,.975))))
  }
  colnames(out) <- c('estimate','lower','upper');out
}

compare_job <- function(j,repo,root,study,code,bundle,grid) {
  sdir <- file.path(root,'jobs',j$job);rdir <- study_job_dir(study,j$job)
  if(!file.exists(file.path(rdir,'result.rds'))) stop('Missing refit for ',j$job)
  refit <- readRDS(file.path(rdir,'result.rds'))
  stopifnot(isTRUE(refit$ok),!isTRUE(refit$truth_used),identical(refit$input_md5,j$input_md5),
    identical(refit$settings$warmup,SCHEDULE$warmup),identical(refit$settings$draws,SCHEDULE$draws),identical(refit$settings$thin,SCHEDULE$thin))
  p_new <- readRDS(file.path(rdir,sprintf('attempt-%d-parameters.rds',ATTEMPT)))
  saved <- readRDS(file.path(sdir,'result.rds'));p_old <- saved$parameters
  stopifnot(isTRUE(saved$ok),!isTRUE(saved$diagnostic_pass),saved$settings$warmup==ORIGINAL_LONGER$warmup,
    saved$settings$draws==ORIGINAL_LONGER$draws,identical(dim(p_old$beta)[1:2],dim(p_new$beta)[1:2]))
  score_old <- readRDS(file.path(sdir,'score.rds'))
  d_old <- utils::read.csv(file.path(sdir,'attempt-2-diagnostics.csv'),stringsAsFactors=FALSE)
  d_new <- utils::read.csv(file.path(rdir,sprintf('attempt-%d-diagnostics.csv',ATTEMPT)),stringsAsFactors=FALSE)
  input <- job_input(root,j);S <- as.integer(j$n_species)
  truth <- readRDS(file.path(root,'truth',paste0(j$dataset,'.rds')))
  # 300 test sites, as score-job.R scored the saved fit.
  train_grid <- input$x[c(1,21,41,61,81),,drop=FALSE]
  reference <- matrix(d_new$mean[d_new$block=='new_site_probability'],nrow=5)
  est_new <- code$checked_posterior_mean(p_new,input$test_x,train_grid,reference)$mean
  truth300 <- readRDS(file.path(root,'truth',paste0(j$dataset,'-marginal.rds')))[,seq_len(S),drop=FALSE]
  stopifnot(identical(dim(est_new),dim(score_old$estimate)),max(abs(score_old$truth-truth300))<1e-12)
  e_old <- score_old$estimate-truth300;e_new <- est_new-truth300
  stopifnot(abs(100*mean(e_old)-score_old$overall$bias_pp)<1e-8,abs(100*mean(abs(e_old))-score_old$overall$mae_pp)<1e-8)
  pred <- data.frame(bias300_old=100*mean(e_old),bias300_new=100*mean(e_new),mae300_old=100*mean(abs(e_old)),mae300_new=100*mean(abs(e_new)),
    rmse300_old=100*sqrt(mean(e_old^2)),rmse300_new=100*sqrt(mean(e_new^2)),max_abs_change300=100*max(abs(est_new-score_old$estimate)))
  # Five fixed environments, as export-calibration.R extracted the saved fit.
  raw <- code$basis(grid,j$response)
  x <- as.data.frame(sweep(sweep(as.matrix(raw),2,input$centre,'-'),2,input$spread,'/'))
  true_p <- as.vector(code$marginal_truth(truth,grid)[,seq_len(S),drop=FALSE])
  pi_old <- code$probability_intervals(p_old,x);pi_new <- code$probability_intervals(p_new,x)
  brow <- bundle$rows[bundle$rows$job==j$job & bundle$rows$target=='Marginal probability',,drop=FALSE]
  stopifnot(nrow(brow)==5L*S,max(abs(brow$estimate-pi_old$estimate))<1e-10,max(abs(brow$lower-pi_old$lower))<1e-10,
    max(abs(brow$upper-pi_old$upper))<1e-10,max(abs(brow$truth-true_p))<1e-10,all(!brow$diagnostic_pass))
  covered <- function(pi) mean(pi$lower<=true_p & true_p<=pi$upper)
  gridc <- data.frame(coverage5_old=covered(pi_old),coverage5_new=covered(pi_new),width5_old=mean(pi_old$upper-pi_old$lower),
    width5_new=mean(pi_new$upper-pi_new$lower),bias5_old=100*mean(pi_old$estimate-true_p),bias5_new=100*mean(pi_new$estimate-true_p),
    rmse5_old=100*sqrt(mean((pi_old$estimate-true_p)^2)),rmse5_new=100*sqrt(mean((pi_new$estimate-true_p)^2)),
    max_abs_change5=100*max(abs(pi_new$estimate-pi_old$estimate)))
  # Probit coefficients on the fitted (standardised) scale: no logit truth exists.
  c_old <- coefficient_summary(p_old$beta);c_new <- coefficient_summary(p_new$beta)
  overlap <- c_new$lower<=c_old$upper & c_old$lower<=c_new$upper
  coef <- data.frame(coef_elements=length(c_old$mean),coef_max_abs_mean_change=max(abs(c_new$mean-c_old$mean)),
    coef_rms_mean_change=sqrt(mean((c_new$mean-c_old$mean)^2)),coef_mean_change_in_saved_sd=mean(abs(c_new$mean-c_old$mean)/c_old$sd),
    coef_intervals_overlap=mean(overlap),coef_mean_width_old=mean(c_old$upper-c_old$lower),coef_mean_width_new=mean(c_new$upper-c_new$lower))
  traits <- data.frame(trait_rows=0L,trait_max_abs_change=NA_real_,trait_intervals_overlap=NA_real_)
  if(isTRUE(j$use_traits) && !is.null(input$traits)) {
    t_old <- hmsc_trait_rows(p_old,input);t_new <- hmsc_trait_rows(p_new,input)
    stopifnot(!is.null(score_old$traits),max(abs(t_old[,'estimate']-score_old$traits$estimate))<1e-10)
    traits <- data.frame(trait_rows=nrow(t_old),trait_max_abs_change=max(abs(t_new[,'estimate']-t_old[,'estimate'])),
      trait_intervals_overlap=mean(t_new[,'lower']<=t_old[,'upper'] & t_old[,'lower']<=t_new[,'upper']))
  }
  fails <- function(d,block) sum(!d$pass & d$block==block)
  diag <- data.frame(passed_new=isTRUE(refit$diagnostic_pass),n_fail_old=sum(!d_old$pass),n_fail_new=sum(!d_new$pass),
    coef_fail_old=fails(d_old,'coefficient'),coef_fail_new=fails(d_new,'coefficient'),cov_fail_old=fails(d_old,'covariance'),
    cov_fail_new=fails(d_new,'covariance'),prob_fail_old=fails(d_old,'new_site_probability'),prob_fail_new=fails(d_new,'new_site_probability'),
    trait_fail_old=fails(d_old,'trait'),trait_fail_new=fails(d_new,'trait'),
    max_coef_rhat_old=max(d_old$rhat[d_old$block=='coefficient']),max_coef_rhat_new=max(d_new$rhat[d_new$block=='coefficient']),
    min_coef_ess_old=min(pmin(d_old$bulk_ess,d_old$tail_ess)[d_old$block=='coefficient']),
    min_coef_ess_new=min(pmin(d_new$bulk_ess,d_new$tail_ess)[d_new$block=='coefficient']),
    fit_minutes_new=refit$fit_seconds/60,saved_elapsed_min=j$saved_elapsed_min,n_warnings_new=length(refit$warnings))
  row <- cbind(j[c('job','job_number','scenario','n_sites','n_species','response','use_traits','replicate','badly_unmixed')],diag,pred,gridc,coef,traits)
  new_rows <- brow;new_rows$estimate <- pi_new$estimate;new_rows$lower <- pi_new$lower;new_rows$upper <- pi_new$upper
  new_rows$diagnostic_pass <- isTRUE(refit$diagnostic_pass)
  new_pred <- bundle$prediction[bundle$prediction$job==j$job,,drop=FALSE];stopifnot(nrow(new_pred)==1L)
  new_pred$bias <- mean(e_new);new_pred$mae <- mean(abs(e_new));new_pred$rmse <- sqrt(mean(e_new^2))
  new_pred$diagnostic_pass <- isTRUE(refit$diagnostic_pass)
  list(row=row,rows=new_rows,prediction=new_pred)
}

cell_prediction <- function(pred) {
  id <- do.call(paste,c(pred[CELL_KEYS],sep='|'))
  do.call(rbind,lapply(split(pred,id),function(z) cbind(z[1,CELL_KEYS,drop=FALSE],communities=nrow(z),
    flagged_communities=sum(!z$diagnostic_pass),bias=100*mean(z$bias),mae=100*mean(z$mae),rmse=100*sqrt(mean(z$rmse^2)))))
}

compare_all <- function(repo,root,study,out=NULL) {
  code <- load_code(repo);jobs <- read_jobs(repo)
  bundle <- readRDS(file.path(repo,BUNDLE_REL))
  stopifnot(identical(unname(tools::md5sum(file.path(root,'status.csv'))),unname(bundle$provenance$archive_status)))
  grid <- code$fixed_grid();rownames(grid) <- GRID_LABELS
  parts <- lapply(seq_len(nrow(jobs)),function(i) {
    r <- compare_job(jobs[i,],repo,root,study,code,bundle,grid)
    cat(sprintf('%d/%d %s: diagnostics %s; coefficient failures %d -> %d\n',i,nrow(jobs),jobs$job[i],
      if(r$row$passed_new) 'passed' else 'flagged',r$row$coef_fail_old,r$row$coef_fail_new));flush.console()
    r
  })
  by_job <- do.call(rbind,lapply(parts,`[[`,'row'));rownames(by_job) <- NULL
  # Lesson 6's Hmsc cells, saved and with the refits substituted.
  rows_h <- bundle$rows[bundle$rows$package=='Hmsc',,drop=FALSE]
  rows_sub <- rbind(rows_h[!rows_h$job %in% jobs$job,,drop=FALSE],do.call(rbind,lapply(parts,`[[`,'rows')))
  stopifnot(nrow(rows_sub)==nrow(rows_h))
  s_old <- code$calibration_group_summary(rows_h);s_new <- code$calibration_group_summary(rows_sub)
  keys <- c(CELL_KEYS,'target','term','population')
  b_old <- bundle$summary[bundle$summary$package=='Hmsc',,drop=FALSE]
  k_b <- do.call(paste,c(b_old[keys],sep='|'));k_o <- do.call(paste,c(s_old[keys],sep='|'));k_n <- do.call(paste,c(s_new[keys],sep='|'))
  stopifnot(setequal(k_b,k_o),all(abs(b_old$coverage[match(k_o,k_b)]-s_old$coverage)<1e-10,na.rm=TRUE),
    all(abs(b_old$bias[match(k_o,k_b)]-s_old$bias)<1e-10),all(abs(b_old$rmse[match(k_o,k_b)]-s_old$rmse)<1e-10))
  m <- match(k_o,k_n);stopifnot(!anyNA(m))
  prob_cells <- cbind(s_old[keys],communities_old=s_old$communities,communities_new=s_new$communities[m],
    flagged_old=s_old$flagged_communities,flagged_new=s_new$flagged_communities[m],
    coverage_old=100*s_old$coverage,coverage_new=100*s_new$coverage[m],coverage_mcse_old=100*s_old$coverage_mcse,
    bias_old=100*s_old$bias,bias_new=100*s_new$bias[m],rmse_old=100*s_old$rmse,rmse_new=100*s_new$rmse[m],
    width_old=s_old$width,width_new=s_new$width[m])
  prob_cells$coverage_change <- prob_cells$coverage_new-prob_cells$coverage_old
  prob_cells$bias_change <- prob_cells$bias_new-prob_cells$bias_old
  prob_cells$rmse_change <- prob_cells$rmse_new-prob_cells$rmse_old
  prob_cells$touched <- vapply(seq_len(nrow(prob_cells)),function(i) any(do.call(paste,c(jobs[CELL_KEYS],sep='|'))==
    do.call(paste,c(prob_cells[i,CELL_KEYS],sep='|'))),logical(1))
  prob_cells$flag <- with(prob_cells,(is.finite(coverage_change) & abs(coverage_change)>=FLAG_COVERAGE_POINTS) |
    abs(bias_change)>=FLAG_ERROR_POINTS | abs(rmse_change)>=FLAG_ERROR_POINTS)
  rownames(prob_cells) <- NULL
  pred_h <- bundle$prediction[bundle$prediction$package=='Hmsc',,drop=FALSE]
  pred_sub <- rbind(pred_h[!pred_h$job %in% jobs$job,,drop=FALSE],do.call(rbind,lapply(parts,`[[`,'prediction')))
  stopifnot(nrow(pred_sub)==nrow(pred_h))
  c_old <- cell_prediction(pred_h);c_new <- cell_prediction(pred_sub)
  k_o <- do.call(paste,c(c_old[CELL_KEYS],sep='|'));k_n <- do.call(paste,c(c_new[CELL_KEYS],sep='|'));m <- match(k_o,k_n)
  pred_cells <- cbind(c_old[CELL_KEYS],communities=c_old$communities,flagged_old=c_old$flagged_communities,flagged_new=c_new$flagged_communities[m],
    bias_old=c_old$bias,bias_new=c_new$bias[m],mae_old=c_old$mae,mae_new=c_new$mae[m],rmse_old=c_old$rmse,rmse_new=c_new$rmse[m])
  pred_cells$bias_change <- pred_cells$bias_new-pred_cells$bias_old
  pred_cells$mae_change <- pred_cells$mae_new-pred_cells$mae_old
  pred_cells$rmse_change <- pred_cells$rmse_new-pred_cells$rmse_old
  pred_cells$flag <- with(pred_cells,abs(bias_change)>=FLAG_ERROR_POINTS | abs(mae_change)>=FLAG_ERROR_POINTS | abs(rmse_change)>=FLAG_ERROR_POINTS)
  rownames(pred_cells) <- NULL
  summary <- data.frame(jobs=nrow(by_job),now_passing=sum(by_job$passed_new),coefficient_failures_remaining=sum(by_job$coef_fail_new>0),
    badly_unmixed_now_passing=sum(by_job$passed_new & by_job$badly_unmixed),
    max_coef_rhat_old=max(by_job$max_coef_rhat_old),max_coef_rhat_new=max(by_job$max_coef_rhat_new),
    max_abs_change300=max(by_job$max_abs_change300),max_abs_change5=max(by_job$max_abs_change5),
    mean_coef_intervals_overlap=mean(by_job$coef_intervals_overlap),max_coef_mean_change_in_saved_sd=max(by_job$coef_mean_change_in_saved_sd),
    probability_cells=nrow(prob_cells),probability_cells_touched=sum(prob_cells$touched),probability_cells_flagged=sum(prob_cells$flag),
    max_coverage_change=max(abs(prob_cells$coverage_change),na.rm=TRUE),prediction_cells=nrow(pred_cells),prediction_cells_flagged=sum(pred_cells$flag),
    max_prediction_rmse_change=max(abs(pred_cells$rmse_change)),mean_fit_minutes=mean(by_job$fit_minutes_new),
    hmsc_version=HMSC_VERSION,compared=format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),stringsAsFactors=FALSE)
  out <- out %||% file.path(repo,STUDY_REL,'results')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(by_job,file.path(out,'comparison-by-job.csv'),row.names=FALSE)
  utils::write.csv(prob_cells,file.path(out,'lesson6-probability-cells.csv'),row.names=FALSE)
  utils::write.csv(pred_cells,file.path(out,'lesson6-prediction-cells.csv'),row.names=FALSE)
  utils::write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
  changed <- summary$probability_cells_flagged>0 || summary$prediction_cells_flagged>0
  writeLines(c(paste0('lesson6_hmsc_numbers_change: ',if(changed) 'yes' else 'no'),
    paste0('refits_now_passing: ',summary$now_passing,' of ',summary$jobs),
    paste0('flagged_probability_cells: ',summary$probability_cells_flagged,' of ',summary$probability_cells),
    paste0('flagged_prediction_cells: ',summary$prediction_cells_flagged,' of ',summary$prediction_cells)),file.path(out,'decision.txt'))
  cat(sprintf('Refits %d: %d now pass the diagnostics, %d keep a coefficient failure; max coefficient Rhat %.3f -> %.3f; Lesson 6 Hmsc cells flagged: %d of %d probability cells (max coverage change %.2f points), %d of %d prediction cells; results in %s\n',
    summary$jobs,summary$now_passing,summary$coefficient_failures_remaining,summary$max_coef_rhat_old,summary$max_coef_rhat_new,
    summary$probability_cells_flagged,summary$probability_cells,summary$max_coverage_change,summary$prediction_cells_flagged,summary$prediction_cells,out))
  invisible(list(by_job=by_job,probability_cells=prob_cells,prediction_cells=pred_cells,summary=summary))
}

compare_main <- function(args) {
  o <- parse_options(args,known=c('repo','root','study','out'),required=c('repo','root','study'))
  repo <- normalizePath(o$repo,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this compare.R')
  compare_all(repo,normalizePath(o$root,mustWork=TRUE),normalizePath(o$study,mustWork=TRUE),o$out)
  0L
}

if(sys.nframe()==0L) {
  status <- tryCatch(compare_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
