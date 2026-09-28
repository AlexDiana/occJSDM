args<-commandArgs(TRUE)
stopifnot(length(args)>=3L)
repo<-normalizePath(args[1]);study<-normalizePath(args[2]);mode<-args[3]
workers<-if(length(args)>3L)as.integer(args[4]) else 2L
stopifnot(mode %in% c('pilot','all','summarise'),workers %in% 1:2)
folder<-file.path(repo,'dev/simstudy/spatial-amplitude-prior/diagnosis')
out<-file.path(study,'diagnosis-v1');dir.create(out,showWarnings=FALSE)
.libPaths(c(file.path(study,'library'),.libPaths()))
Rcpp::sourceCpp(file.path(folder,'ellipse.cpp'),cacheDir=file.path(out,'cpp-cache'))
source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/robust.R'))
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/score.R'))
validated<-read.csv(file.path(out,'validation-source-md5.csv'))
stopifnot(identical(unname(tools::md5sum(validated$file)),validated$md5))
sources<-c(file.path(folder,c('PLAN.md','ellipse.cpp','run.R','check.R')),
  file.path(repo,'dev/simstudy/spatial-amplitude-prior/robust.R'),
  file.path(repo,'dev/simstudy/spatial-targeted-recheck/score.R'))
hashes<-data.frame(file=sources,md5=unname(tools::md5sum(sources)))
atomic<-function(object,path) {
  temp<-paste0(path,'.tmp');saveRDS(object,temp);stopifnot(file.rename(temp,path))
}
log<-function(...)cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),...,'\n')
prepared<-file.path(out,'prepared');dir.create(prepared,showWarnings=FALSE)
reference<-file.path(dirname(study),'spatial-targeted-20260927')
inputs<-sort(list.files(file.path(reference,'inputs'),pattern='[.]rds$',full.names=TRUE))
stopifnot(length(inputs)==9L)
geometry<-information<-list()
for(path in inputs) {
  input<-readRDS(path);t<-input$truth;n<-nrow(t$field);S<-ncol(t$field)
  community<-sub('[.]rds$','',basename(path));key<-paste0(community,'-binary-k100')
  baseline<-readRDS(file.path(study,'inverse_gamma/baseline',paste0(key,'-initial-result.rds')))
  stopifnot(unname(tools::md5sum(path))==baseline$job$input_md5,
    unname(tools::md5sum(baseline$source_fit))==baseline$source_fit_md5)
  fit<-readRDS(baseline$source_fit)$fit
  stopifnot(max(abs(fit$Xs-t$Xs))<1e-12,fit$infos$ps==100L,
    fit$infos$l_s_grid[t$grid_index]==t$range,identical(input$data$binary$OTU,t$z))
  K<-exp(-as.matrix(dist(t$Xs))^2/(2*t$range^2))
  Q<-K%*%solve(K+diag(1e-5,n),K);Q<-(Q+t(Q))/2
  independent<-independent_bases(fit)[[t$grid_index]]
  native<-suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid,fit$infos$list_Xs))
  H<-occJSDM:::KsBproduct(native$Ks_all[,,t$grid_index],diag(n),fit$infos$list_Xs$Xs_centers)
  checks<-c(native_covariance=max(abs(Q-tcrossprod(H))),
    independent_covariance=max(abs(Q-tcrossprod(independent))),
    generating_covariance=max(abs(Q-K-diag(1e-10,n))),minimum_eigenvalue=min(eigen(Q,symmetric=TRUE,only.values=TRUE)$values))
  stopifnot(checks['native_covariance']<1e-10,checks['independent_covariance']<1e-10,
    checks['minimum_eigenvalue']>0)
  seed<-93000000L+t$grid_index*1000L+input$settings$replicate
  set.seed(seed);extra<-matrix(rbinom(n*S,9L,as.vector(t$psi)),n,S)
  successes<-unname(t$z)+extra
  stopifnot(all(successes>=t$z),all(successes-t$z<=9L),all(successes<=10L))
  object<-list(input=path,input_md5=unname(tools::md5sum(path)),community=community,
    truth=t,replicate=input$settings$replicate,covariance=Q,covariance_checks=checks,
    extra_seed=seed,extra_successes=extra,ten_successes=successes,
    reference_fit=baseline$source_fit,reference_fit_md5=baseline$source_fit_md5)
  destination<-file.path(prepared,paste0(community,'.rds'))
  if(file.exists(destination))stopifnot(identical(readRDS(destination),object)) else atomic(object,destination)
  neighbour<-K;diag(neighbour)<-0
  projection<-qr.fitted(qr(cbind(1,t$X)),t$field)
  geometry[[community]]<-data.frame(community=community,range=t$range,
    mean_neighbours_above_half=mean(rowSums(neighbour>.5)),
    fraction_without_half_neighbour=mean(rowSums(neighbour>.5)==0),
    median_nearest_correlation=median(apply(neighbour,1,max)),
    effective_covariance_rank=sum(diag(K))^2/sum(K^2),
    intercept_environment_energy=mean(colSums(projection^2)/colSums(t$field^2)),
    as.list(checks),input_md5=object$input_md5,extra_seed=seed)
  for(s in seq_len(S)) {
    w<-t$psi[,s]*(1-t$psi[,s]);f<-t$field[,s];fc<-f-mean(f)
    weighted_residual<-sqrt(w)*f-qr.fitted(qr(sqrt(w)*cbind(1,t$X)),sqrt(w)*f)
    information[[paste(community,s)]]<-data.frame(community=community,species=s,
      prevalence=t$target_prevalence[s],occupied=sum(t$z[,s]),
      intercept_information=sum(w),known_pattern_information=sum(w*f^2),
      pattern_information_after_intercept_environment=sum(weighted_residual^2),
      truth_vs_zero_loglik=sum(dbinom(t$z[,s],1,t$psi[,s],log=TRUE)-
        dbinom(t$z[,s],1,plogis(t$B0[s]+t$X*t$B[s]),log=TRUE)),
      zero_raw_rmse=sqrt(mean(f^2)),zero_centred_rmse=sqrt(mean(fc^2)),
      true_intercept=t$B0[s])
  }
  rm(fit,native);gc(verbose=FALSE)
}
write.csv(do.call(rbind,geometry),file.path(out,'geometry.csv'),row.names=FALSE)
write.csv(do.call(rbind,information),file.path(out,'information.csv'),row.names=FALSE)
write.csv(hashes,file.path(out,'source-md5.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))

score<-function(x,t,s,offset,unknown_intercept) {
  n<-nrow(t$field);ni<-dim(x)[2];nc<-dim(x)[3]
  field<-matrix(x[seq_len(n),,,drop=FALSE],n,ni*nc)
  probabilities<-plogis(sweep(field,1,offset,'+')+
    if(unknown_intercept)rep(x[n+1,,],each=n) else 0)
  summaries<-t(apply(field,1,quantile,probs=c(.025,.5,.975),names=FALSE))
  med<-summaries[,2];truth<-t$field[,s];tc<-truth-mean(truth);fc<-med-mean(med)
  diagnostics<-list()
  add<-function(values,label) {
    diagnostics[[length(diagnostics)+1L]]<<-data.frame(quantity=label,
      as.list(robust_trace_diagnostics(matrix(values,ni,nc))))
  }
  for(i in seq_len(n)) {
    add(field[i,],paste0('field_',i));add(probabilities[i,],paste0('probability_',i))
  }
  add(colMeans(field),'field_average');add(sqrt(colMeans(field^2)),'field_rms')
  add(colSums(field*tc)/sum(tc^2),'field_truth_projection')
  if(unknown_intercept)add(x[n+1,,],'intercept')
  diagnostics<-do.call(rbind,diagnostics)
  pmean<-rowMeans(probabilities)
  pq<-t(apply(probabilities,1,quantile,probs=c(.025,.975),names=FALSE))
  metrics<-data.frame(raw_rmse=sqrt(mean((med-truth)^2)),
    centred_rmse=sqrt(mean((fc-tc)^2)),centred_correlation=cor(fc,tc),
    centred_slope=sum(fc*tc)/sum(tc^2),median_field_rms=sqrt(mean(med^2)),
    field_coverage=mean(summaries[,1]<=truth & summaries[,3]>=truth),
    field_interval_width=mean(summaries[,3]-summaries[,1]),
    occupancy_bias=mean(pmean-t$psi[,s]),occupancy_mae=mean(abs(pmean-t$psi[,s])),
    occupancy_coverage=mean(pq[,1]<=t$psi[,s] & pq[,2]>=t$psi[,s]),
    flag_count=sum(robust_flag_rows(diagnostics)),max_rhat=max(diagnostics$rhat),
    min_ess=min(as.matrix(diagnostics[c('ess_bulk','ess_median','ess_q025','ess_q975')])) )
  list(metrics=metrics,diagnostics=diagnostics,field_quantiles=summaries,
    probability_mean=pmean,probability_quantiles=pq)
}
jobs<-expand.grid(community=names(geometry),configuration=LETTERS[1:4],species=1:8,
  stringsAsFactors=FALSE)
work<-function(i) {
  job<-jobs[i,];input_path<-file.path(prepared,paste0(job$community,'.rds'))
  input<-readRDS(input_path);t<-input$truth;s<-job$species;n<-nrow(t$field)
  unknown<-job$configuration=='C';d<-n+unknown
  amplitude<-if(job$configuration=='A')sqrt(1/qgamma(.5,shape=10)) else 1
  prior_covariance<-diag(d);prior_covariance[seq_len(n),seq_len(n)]<-amplitude^2*input$covariance
  lower<-t(chol(prior_covariance))
  offset<-t$X*t$B[s]+if(unknown)0 else t$B0[s]
  trials<-if(job$configuration=='D')10L else 1L
  y<-if(trials==10L)input$ten_successes[,s] else t$z[,s]
  seed<-94000000L+t$grid_index*100000L+input$replicate*10000L+
    match(job$configuration,LETTERS)*1000L+s*10L
  key<-paste(job$community,job$configuration,sprintf('species%02d',s),sep='-')
  perform<-function(phase) {
    dest<-file.path(out,paste0(key,'-',phase,'.rds'))
    nburn<-if(phase=='initial')1000L else 2000L
    niter<-if(phase=='initial')2000L else 8000L
    if(file.exists(dest)) {
      result<-readRDS(dest)
      stopifnot(identical(result$hashes,hashes),identical(result$job,job),
        identical(result$input_md5,unname(tools::md5sum(input_path))),
        result$nburn==nburn,result$niter==niter,result$seed==seed)
      return(result)
    }
    set.seed(seed);initial<-cbind(rep(0,d),lower%*%matrix(rnorm(d*3L),d,3L))
    started<-Sys.time();fit<-ellipse_draws_cpp(lower,y,trials,offset,unknown,initial,nburn,niter)
    scored<-score(fit$draws,t,s,offset,unknown)
    result<-c(list(job=job,input_md5=unname(tools::md5sum(input_path)),hashes=hashes,
      seed=seed,initial=initial,phase=phase,nburn=nburn,niter=niter,nchain=4L,
      amplitude=amplitude,trials=trials,elapsed=as.numeric(difftime(Sys.time(),started,units='secs')),
      draws=fit$draws,mean_proposals=fit$mean_proposals),scored)
    stopifnot(identical(unname(tools::md5sum(hashes$file)),hashes$md5))
    atomic(result,dest)
    log(key,phase,'complete; seconds',round(result$elapsed,1),'flags',result$metrics$flag_count)
    result
  }
  a<-perform('initial');initial_metrics<-a$metrics
  if(a$metrics$flag_count>0L)a<-perform('long')
  data.frame(job,phase=a$phase,amplitude=a$amplitude,trials=a$trials,
    elapsed=a$elapsed,initial_flag_count=initial_metrics$flag_count,a$metrics,
    file=file.path(out,paste0(key,'-',a$phase,'.rds')))
}
if(mode=='pilot') {
  which_jobs<-which(jobs$community=='range6-rep01' & jobs$species==6L & jobs$configuration %in% c('B','D'))
} else which_jobs<-seq_len(nrow(jobs))
if(mode=='summarise') {
  for(i in which_jobs) {
    job<-jobs[i,];key<-paste(job$community,job$configuration,sprintf('species%02d',job$species),sep='-')
    stopifnot(file.exists(file.path(out,paste0(key,'-initial.rds'))))
    a<-readRDS(file.path(out,paste0(key,'-initial.rds')))
    if(a$metrics$flag_count>0L)stopifnot(file.exists(file.path(out,paste0(key,'-long.rds'))))
  }
}
results<-if(workers==1L)lapply(which_jobs,work) else parallel::mclapply(which_jobs,work,
  mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
table<-do.call(rbind,results);table$result_md5<-unname(tools::md5sum(table$file))
write.csv(table,file.path(out,paste0(mode,'-selected.csv')),row.names=FALSE)
if(mode!='pilot')stopifnot(nrow(table)==288L)
log(mode,'finished;',nrow(table),'conditional posteriors;',sum(table$flag_count>0),'flagged')
