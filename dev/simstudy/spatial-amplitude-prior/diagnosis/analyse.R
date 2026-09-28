args<-commandArgs(TRUE);stopifnot(length(args)==3L)
repo<-normalizePath(args[1]);out<-normalizePath(args[2]);workers<-as.integer(args[3])
destination<-file.path(out,'summary');dir.create(destination,showWarnings=FALSE)
success_marker<-file.path(destination,'audit-complete.rds')
if(file.exists(success_marker))unlink(success_marker)
stopifnot(!file.exists(success_marker))
stopifnot(workers %in% 1:2)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/robust.R'))
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
analysis_files<-c(file.path(repo,'dev/simstudy/spatial-amplitude-prior/diagnosis/analyse.R'),
  file.path(repo,'dev/simstudy/spatial-amplitude-prior/robust.R'),
  file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
analysis_hashes<-data.frame(file=analysis_files,md5=unname(tools::md5sum(analysis_files)))
frozen<-read.csv(file.path(out,'source-md5.csv'))
stopifnot(identical(unname(tools::md5sum(frozen$file)),frozen$md5))
selected<-read.csv(file.path(out,'all-selected.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(selected)==288L,!anyDuplicated(selected[c('community','configuration','species')]),
  all(table(selected$community,selected$configuration)==8L),
  setequal(selected$configuration,LETTERS[1:4]),setequal(selected$species,1:8),
  identical(unname(tools::md5sum(selected$file)),selected$result_md5))
audit_one<-function(i) {
  row<-selected[i,];a<-readRDS(row$file)
  prepared_path<-file.path(out,'prepared',paste0(row$community,'.rds'))
  p<-readRDS(prepared_path);truth<-p$truth;s<-row$species;n<-nrow(truth$field)
  original<-readRDS(p$input)
  stopifnot(identical(truth,original$truth),p$replicate==original$settings$replicate,
    sub('[.]rds$','',basename(p$input))==row$community)
  kernel<-exp(-as.matrix(dist(truth$Xs))^2/(2*truth$range^2))
  basis<-t(forwardsolve(t(chol(kernel+diag(1e-5,n))),t(kernel)))
  covariance_difference<-max(abs(tcrossprod(basis)-p$covariance))
  stopifnot(covariance_difference<1e-10)
  stopifnot(identical(a$hashes,frozen),a$input_md5==unname(tools::md5sum(prepared_path)),
    p$input_md5==unname(tools::md5sum(p$input)),
    a$job$community==row$community,a$job$configuration==row$configuration,a$job$species==s,
    a$nchain==4L,a$phase==row$phase,
    a$niter==if(a$phase=='long')8000L else 2000L,
    a$nburn==if(a$phase=='long')2000L else 1000L)
  initial_file<-sub('-(initial|long)[.]rds$','-initial.rds',row$file)
  initial<-readRDS(initial_file)
  initial_flags<-robust_flag_rows(initial$diagnostics)|(initial$diagnostics$mean_required &
    (!is.finite(initial$diagnostics$ess_mean)|initial$diagnostics$ess_mean<100))
  stopifnot(identical(initial$job,a$job),identical(initial$hashes,a$hashes),
    identical(initial$input_md5,a$input_md5),initial$metrics$flag_count==row$initial_flag_count,
    sum(initial_flags)==row$initial_flag_count,initial$nchain==4L,
    initial$niter==2000L,initial$nburn==1000L,
    identical(a$phase,if(initial$metrics$flag_count>0)'long' else 'initial'))
  seed<-93000000L+truth$grid_index*1000L+p$replicate
  set.seed(seed);fresh<-matrix(rbinom(n*ncol(truth$field),9L,as.vector(truth$psi)),n)
  stopifnot(seed==p$extra_seed,identical(fresh,p$extra_successes),
    identical(unname(truth$z)+fresh,p$ten_successes))
  unknown<-row$configuration=='C';d<-n+unknown
  amplitude<-if(row$configuration=='A')sqrt(1/qgamma(.5,shape=10)) else 1
  stopifnot(a$amplitude==amplitude,a$trials==if(row$configuration=='D')10L else 1L,
    identical(dim(a$draws),c(d,a$niter,4L)),all(is.finite(a$draws)))
  fit_seed<-94000000L+truth$grid_index*100000L+p$replicate*10000L+
    match(row$configuration,LETTERS)*1000L+s*10L
  covariance<-diag(d);covariance[seq_len(n),seq_len(n)]<-amplitude^2*p$covariance
  set.seed(fit_seed);starts<-cbind(rep(0,d),t(chol(covariance))%*%matrix(rnorm(d*3L),d,3L))
  stopifnot(fit_seed==a$seed,identical(starts,a$initial))
  fields<-a$draws[seq_len(n),,,drop=FALSE]
  probabilities<-array(NA_real_,dim(fields))
  for(ch in 1:4)for(it in seq_len(a$niter)) {
    intercept<-if(unknown)a$draws[n+1L,it,ch] else truth$B0[s]
    probabilities[,it,ch]<-plogis(fields[,it,ch]+intercept+truth$X*truth$B[s])
  }
  q<-t(vapply(seq_len(n),function(j)quantile(fields[j,,],c(.025,.5,.975),names=FALSE),numeric(3)))
  pq<-t(vapply(seq_len(n),function(j)quantile(probabilities[j,,],c(.025,.975),names=FALSE),numeric(2)))
  pm<-vapply(seq_len(n),function(j)mean(probabilities[j,,]),numeric(1))
  difference<-max(abs(q-a$field_quantiles),abs(pq-a$probability_quantiles),abs(pm-a$probability_mean))
  tc<-truth$field[,s]-mean(truth$field[,s]);fc<-q[,2]-mean(q[,2])
  metrics<-c(raw_rmse=sqrt(mean((q[,2]-truth$field[,s])^2)),
    centred_rmse=sqrt(mean((fc-tc)^2)),centred_correlation=cor(fc,tc),
    centred_slope=sum(fc*tc)/sum(tc^2),median_field_rms=sqrt(mean(q[,2]^2)),
    field_coverage=mean(q[,1]<=truth$field[,s]&q[,3]>=truth$field[,s]),
    field_interval_width=mean(q[,3]-q[,1]),occupancy_bias=mean(pm-truth$psi[,s]),
    occupancy_mae=mean(abs(pm-truth$psi[,s])),
    occupancy_coverage=mean(pq[,1]<=truth$psi[,s]&pq[,2]>=truth$psi[,s]))
  metric_difference<-max(abs(metrics-unlist(a$metrics[names(metrics)])))
  # Check numerical reconstruction independently above; reproduce the recorded
  # arithmetic for rank diagnostics, whose folded ties can change at roundoff.
  field_matrix<-matrix(fields,n,a$niter*4L)
  offset<-truth$X*truth$B[s]+if(unknown)0 else truth$B0[s]
  probability_matrix<-plogis(sweep(field_matrix,1,offset,'+')+
    if(unknown)rep(a$draws[n+1,,],each=n) else 0)
  probability_draw_difference<-max(abs(probability_matrix-as.vector(probabilities)))
  traces<-list()
  for(j in seq_len(n)) {
    traces[[paste0('field_',j)]]<-fields[j,,]
    traces[[paste0('probability_',j)]]<-matrix(probability_matrix[j,],a$niter,4L)
  }
  traces$field_average<-matrix(colMeans(field_matrix),a$niter,4L)
  traces$field_rms<-matrix(sqrt(colMeans(field_matrix^2)),a$niter,4L)
  traces$field_truth_projection<-matrix(colSums(field_matrix*tc)/sum(tc^2),a$niter,4L)
  traces$probability_average<-matrix(colMeans(probability_matrix),a$niter,4L)
  if(unknown)traces$intercept<-a$draws[n+1L,,]
  stopifnot(identical(names(traces),a$diagnostics$quantity))
  rechecked<-do.call(rbind,lapply(names(traces),function(label) {
    z<-traces[[label]];need_mean<-startsWith(label,'probability_')
    data.frame(quantity=label,as.list(robust_trace_diagnostics(z)),mean_required=need_mean,
      ess_mean=if(need_mean)as.numeric(posterior::ess_mean(z)) else NA_real_,
      mcse_mean=if(need_mean)as.numeric(posterior::mcse_mean(z)) else NA_real_)
  }))
  numeric_columns<-names(rechecked)[vapply(rechecked,is.numeric,logical(1))]
  stopifnot(identical(is.na(rechecked[numeric_columns]),is.na(a$diagnostics[numeric_columns])))
  diagnostic_difference<-max(abs(as.matrix(rechecked[numeric_columns])-as.matrix(a$diagnostics[numeric_columns])),na.rm=TRUE)
  flags<-robust_flag_rows(rechecked)|(rechecked$mean_required &
    (!is.finite(rechecked$ess_mean)|rechecked$ess_mean<100))
  stopifnot(sum(flags)==row$flag_count,sum(flags)==a$metrics$flag_count,
    difference<1e-12,metric_difference<1e-12,diagnostic_difference<1e-12,
    probability_draw_difference<1e-12)
  sites<-data.frame(community=row$community,configuration=row$configuration,species=s,
    site=seq_len(n),prevalence=truth$target_prevalence[s],grid_index=truth$grid_index,
    x=truth$Xs[,1],y=truth$Xs[,2],field_truth=truth$field[,s],field_median=q[,2],
    field_lower=q[,1],field_upper=q[,3],probability_truth=truth$psi[,s],
    probability_mean=pm,probability_lower=pq[,1],probability_upper=pq[,2])
  for(ch in 1:4) {
    sites[[paste0('field_chain',ch)]]<-apply(fields[,,ch],1,median)
    sites[[paste0('probability_chain',ch)]]<-rowMeans(probabilities[,,ch])
  }
  audit<-data.frame(community=row$community,configuration=row$configuration,species=s,
    result_md5=row$result_md5,input_md5=a$input_md5,summary_difference=difference,
    metric_difference=metric_difference,diagnostic_difference=diagnostic_difference,
    probability_draw_difference=probability_draw_difference,
    covariance_difference=covariance_difference,
    audited_traces=nrow(rechecked),flag_count=sum(flags))
  list(sites=sites,audit=audit,diagnostics=cbind(row[c('community','configuration','species')],rechecked))
}
results<-parallel::mclapply(seq_len(nrow(selected)),audit_one,mc.cores=workers,
  mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
sites<-do.call(rbind,lapply(results,`[[`,'sites'))
audits<-do.call(rbind,lapply(results,`[[`,'audit'))
diagnostics<-do.call(rbind,lapply(results,`[[`,'diagnostics'))
stopifnot(nrow(sites)==28800L,nrow(audits)==288L)
groups<-chain_groups<-chain_extrema<-list()
for(community in unique(sites$community))for(configuration in LETTERS[1:4]) {
  a<-sites[sites$community==community & sites$configuration==configuration,]
  for(group in c('all',as.character(c(.01,.05,.25,.75)))) {
    b<-if(group=='all')a else a[a$prevalence==as.numeric(group),]
    fc<-b$field_median-ave(b$field_median,b$species,FUN=mean)
    tc<-b$field_truth-ave(b$field_truth,b$species,FUN=mean)
    groups[[length(groups)+1L]]<-data.frame(community=community,configuration=configuration,
      group=group,grid_index=b$grid_index[1],raw_rmse=sqrt(mean((b$field_median-b$field_truth)^2)),
      centred_rmse=sqrt(mean((fc-tc)^2)),centred_correlation=cor(fc,tc),
      centred_slope=sum(fc*tc)/sum(tc^2),zero_centred_rmse=sqrt(mean(tc^2)),
      field_coverage=mean(b$field_lower<=b$field_truth & b$field_upper>=b$field_truth),
      field_interval_width=mean(b$field_upper-b$field_lower),
      occupancy_bias=mean(b$probability_mean-b$probability_truth),
      occupancy_mae=mean(abs(b$probability_mean-b$probability_truth)),
      occupancy_coverage=mean(b$probability_lower<=b$probability_truth & b$probability_upper>=b$probability_truth))
    for(ch in 1:4) {
      med<-b[[paste0('field_chain',ch)]];centred<-med-ave(med,b$species,FUN=mean)
      chain_groups[[length(chain_groups)+1L]]<-data.frame(community=community,configuration=configuration,
        group=group,grid_index=b$grid_index[1],chain=ch,
        raw_rmse=sqrt(mean((med-b$field_truth)^2)),centred_rmse=sqrt(mean((centred-tc)^2)),
        occupancy_mae=mean(abs(b[[paste0('probability_chain',ch)]]-b$probability_truth)))
    }
    # Species were sampled independently: chain labels need not align across
    # species. Add species-wise extrema before taking the community RMSE.
    extremes<-vapply(split(seq_len(nrow(b)),b$species),function(ids) {
      raw_error<-centred_error<-absolute_error<-numeric(4)
      for(ch in 1:4) {
        med<-b[[paste0('field_chain',ch)]][ids]
        raw_error[ch]<-sum((med-b$field_truth[ids])^2)
        centred_error[ch]<-sum((med-mean(med)-tc[ids])^2)
        absolute_error[ch]<-sum(abs(b[[paste0('probability_chain',ch)]][ids]-b$probability_truth[ids]))
      }
      c(raw_min=min(raw_error),raw_max=max(raw_error),centred_min=min(centred_error),
        centred_max=max(centred_error),mae_min=min(absolute_error),mae_max=max(absolute_error))
    },numeric(6))
    bounds<-rowSums(extremes)/nrow(b)
    chain_extrema[[length(chain_extrema)+1L]]<-data.frame(community=community,configuration=configuration,
      group=group,grid_index=b$grid_index[1],raw_rmse_min=sqrt(bounds['raw_min']),
      raw_rmse_max=sqrt(bounds['raw_max']),centred_rmse_min=sqrt(bounds['centred_min']),
      centred_rmse_max=sqrt(bounds['centred_max']),occupancy_mae_min=bounds['mae_min'],
      occupancy_mae_max=bounds['mae_max'])
  }
}
groups<-do.call(rbind,groups)
chain_groups<-do.call(rbind,chain_groups)
chain_extrema<-do.call(rbind,chain_extrema)
metrics<-setdiff(names(groups),c('community','configuration','group','grid_index'))
aggregate<-do.call(rbind,lapply(split(groups,interaction(groups$configuration,groups$group,drop=TRUE)),
  function(a)do.call(rbind,lapply(metrics,function(m)data.frame(configuration=a$configuration[1],
    group=a$group[1],metric=m,as.list(stratified_mean(a[[m]],a$grid_index)))))))
paired<-do.call(rbind,lapply(c('A','C','D'),function(other) {
  a<-groups[groups$configuration=='B',];b<-groups[groups$configuration==other,]
  z<-merge(a,b,by=c('community','group','grid_index'),suffixes=c('_B','_other'))
  do.call(rbind,lapply(split(z,z$group),function(v)do.call(rbind,lapply(metrics,function(m) {
    delta<-v[[paste0(m,'_other')]]-v[[paste0(m,'_B')]]
    data.frame(comparison=paste0(other,' minus B'),group=v$group[1],metric=m,
      decreased=sum(delta<0),increased=sum(delta>0),as.list(stratified_mean(delta,v$grid_index)))
  }))))
}))
# Descriptive envelope over all observed independent-species chain choices,
# not an MC confidence bound. Pooled-chain medians need not lie inside it.
chain_sensitivity<-list()
for(other in c('A','C','D'))for(group in unique(chain_groups$group))
  for(metric in c('raw_rmse','centred_rmse','occupancy_mae')) {
    envelopes<-lapply(unique(chain_groups$community),function(community) {
      z<-chain_extrema[chain_extrema$community==community & chain_extrema$group==group,]
      base<-z[z$configuration=='B',];alternative<-z[z$configuration==other,]
      stopifnot(nrow(base)==1L,nrow(alternative)==1L)
      c(lower=alternative[[paste0(metric,'_min')]]-base[[paste0(metric,'_max')]],
        upper=alternative[[paste0(metric,'_max')]]-base[[paste0(metric,'_min')]])
    })
    bounds<-do.call(rbind,envelopes)
    chain_sensitivity[[length(chain_sensitivity)+1L]]<-data.frame(comparison=paste0(other,' minus B'),
      group=group,metric=metric,observed_chain_min_change=mean(bounds[,'lower']),
      observed_chain_max_change=mean(bounds[,'upper']))
  }
chain_sensitivity<-do.call(rbind,chain_sensitivity)
stopifnot(identical(unname(tools::md5sum(analysis_hashes$file)),analysis_hashes$md5),
  identical(unname(tools::md5sum(selected$file)),selected$result_md5))
write.csv(audits,file.path(destination,'audit.csv'),row.names=FALSE)
write.csv(diagnostics,file.path(destination,'diagnostics.csv'),row.names=FALSE)
write.csv(sites,file.path(destination,'sites.csv'),row.names=FALSE)
write.csv(groups,file.path(destination,'community-groups.csv'),row.names=FALSE)
write.csv(aggregate,file.path(destination,'aggregate.csv'),row.names=FALSE)
write.csv(paired,file.path(destination,'paired.csv'),row.names=FALSE)
write.csv(chain_groups,file.path(destination,'chain-community-groups.csv'),row.names=FALSE)
write.csv(chain_extrema,file.path(destination,'chain-community-extrema.csv'),row.names=FALSE)
write.csv(chain_sensitivity,file.path(destination,'chain-sensitivity.csv'),row.names=FALSE)
write.csv(analysis_hashes,file.path(destination,'analysis-source-md5.csv'),row.names=FALSE)
artifacts<-c(list.files(destination,pattern='[.]csv$',full.names=TRUE),
  file.path(out,c('all-selected.csv','geometry.csv','information.csv','source-md5.csv','quadrature-checks.csv','validation-source-md5.csv')))
artifact_hashes<-data.frame(file=artifacts,md5=unname(tools::md5sum(artifacts)))
saveRDS(list(selected=selected,analysis_hashes=analysis_hashes,artifact_hashes=artifact_hashes,
  finished=Sys.time()),success_marker)
cat('Audited all 288 conditional posteriors; produced pooled community and prevalence results.\n')
