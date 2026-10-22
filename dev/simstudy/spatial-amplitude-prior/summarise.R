#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
option <- function(n,default=NULL) {
  x<-args[startsWith(args,paste0('--',n,'='))];stopifnot(length(x)<=1L)
  if(!length(x)){if(is.null(default))stop('Missing --',n);return(default)}
  substring(x,nchar(n)+4L)
}
repo<-normalizePath(option('repo'));study<-normalizePath(option('study'))
phase<-option('phase','select');arms<-strsplit(option('arms','binary'),',',TRUE)[[1]]
stopifnot(phase %in% c('select','initial','final'),all(arms %in% c('binary','low','high')))
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/metrics.R'))
source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/analysis.R'))
new_fit_hashes<-readRDS(file.path(study,'half_cauchy/initial/settings.rds'))$fit_hashes
stopifnot(identical(unname(tools::md5sum(names(new_fit_hashes))),unname(new_fit_hashes)))
keys<-unlist(lapply(arms,function(arm)unlist(lapply(c(4,6,8),function(g)
  sprintf('range%d-rep%02d-%s-k100',g,1:3,arm)))))
rows<-list();selected<-list();selection<-list();groups<-list();diagnostics<-list()
for(key in keys) {
  igfile<-file.path(study,'inverse_gamma/baseline',paste0(key,'-initial-result.rds'))
  hcfile<-file.path(study,'half_cauchy/initial',paste0(key,'-result.rds'))
  ig<-read_checked_result(igfile);hc<-read_checked_result(hcfile)
  initial_ig<-ig;initial_hc<-hc
  expected_community<-sub('-(binary|low|high)-k100$','',key)
  expected_arm<-sub('.*-(binary|low|high)-k100$','\\1',key)
  stopifnot(ig$job$key==key,hc$job$key==key,ig$job$community==expected_community,
    hc$job$community==expected_community,ig$job$arm==expected_arm,hc$job$arm==expected_arm,
    identical(ig$job[c('grid_index','replicate')],hc$job[c('grid_index','replicate')]))
  validate_choice(ig,initial_ig,'inverse_gamma','initial',new_fit_hashes)
  validate_choice(hc,initial_hc,'half_cauchy','initial',new_fit_hashes)
  stopifnot(identical(ig$job$input_md5,hc$job$input_md5),identical(ig$mcmc,hc$mcmc),
    identical(hc$priors,list(sigma_bs_prior='half_cauchy',sigma_bs_scale=1)),
    is.null(hc$starts))
  need<-pair_needs_long(ig,hc);chosen<-'initial'
  longig<-file.path(study,'inverse_gamma/baseline',paste0(key,'-long-result.rds'))
  if(!file.exists(longig))longig<-file.path(study,'inverse_gamma/long',paste0(key,'-result.rds'))
  longhc<-file.path(study,'half_cauchy/long',paste0(key,'-result.rds'))
  selection[[key]]<-data.frame(key=key,needs_long=need,
    initial_ig_reasons=paste(ig$reasons,collapse='; '),
    initial_hc_reasons=paste(hc$reasons,collapse='; '),
    run_ig_long=need && !file.exists(longig),run_hc_long=need && !file.exists(longhc))
  if(phase=='select')next
  if(phase=='final' && need) {
    ig<-read_checked_result(longig);hc<-read_checked_result(longhc)
    igfile<-longig;hcfile<-longhc;chosen<-'long'
  }
  budget<-if(chosen=='long')list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L) else
    list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L)
  stopifnot(identical(ig$mcmc,budget),identical(hc$mcmc,budget),
    identical(ig$job$input_md5,hc$job$input_md5),is.null(hc$starts))
  validate_choice(ig,initial_ig,'inverse_gamma',chosen,new_fit_hashes)
  validate_choice(hc,initial_hc,'half_cauchy',chosen,new_fit_hashes)
  for(prior in c('inverse_gamma','half_cauchy')) {
    a<-if(prior=='inverse_gamma')ig else hc;file<-if(prior=='inverse_gamma')igfile else hcfile
    id<-paste(key,prior,sep=':');selected[[id]]<-a
    rows[[id]]<-cbind(result_row(a,prior,chosen),result_file=normalizePath(file),
      result_md5=unname(tools::md5sum(file)),fit_file=a$source_fit,fit_md5=a$source_fit_md5)
    prefix<-rows[[id]][c('key','community','grid_index','replicate','arm','prior','phase')]
    groups[[id]]<-cbind(prefix[rep(1,nrow(a$groups)),],a$groups)
    diagnostics[[id]]<-cbind(prefix[rep(1,nrow(a$spatial$diagnostics)),],a$spatial$diagnostics)
  }
}
out<-file.path(study,paste0('summary-',paste(arms,collapse='-'),'-',phase));dir.create(out,showWarnings=FALSE)
provenance<-list(args=args,created=Sys.time(),session=sessionInfo(),
  script_hashes=tools::md5sum(c(file.path(repo,'dev/simstudy/spatial-amplitude-prior',
    c('analysis.R','summarise.R','metrics.R','PLAN.md')),
    file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))))
saveRDS(provenance,file.path(out,'provenance.rds'))
selection<-do.call(rbind,selection);write.csv(selection,file.path(out,'selection.csv'),row.names=FALSE)
cat('Paired longer schedules:',sum(selection$needs_long),'; new IG:',sum(selection$run_ig_long),
  '; new HC:',sum(selection$run_hc_long),'\n')
if(phase=='select')quit(status=0)
rows<-do.call(rbind,rows);pairs<-pair_rows(rows)
write.csv(rows,file.path(out,'fits.csv'),row.names=FALSE)
write.csv(pairs,file.path(out,'paired-communities.csv'),row.names=FALSE)
write.csv(paired_summary(pairs),file.path(out,'paired-summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,groups),file.path(out,'groups.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'spatial-diagnostics.csv'),row.names=FALSE)
saveRDS(selected,file.path(out,'selected-results.rds'))
if(phase=='final' && 'binary' %in% arms) {
  probe_key<-'range6-rep01-binary-k100'
  probefile<-file.path(study,'half_cauchy/starts-initial',paste0(probe_key,'-result.rds'))
  probe<-read_checked_result(probefile)
  probe_initial<-read_checked_result(file.path(study,'half_cauchy/initial',paste0(probe_key,'-result.rds')))
  validate_choice(probe,probe_initial,'half_cauchy','starts-initial',new_fit_hashes)
  if(length(spatial_flags(probe$spatial$diagnostics))) {
    probefile<-file.path(study,'half_cauchy/starts-long',paste0(probe_key,'-result.rds'))
    probe<-read_checked_result(probefile)
    validate_choice(probe,probe_initial,'half_cauchy','starts-long',new_fit_hashes)
  }
  stopifnot(identical(probe$starts,c(.1,.3,1,3)))
  gate<-extension_gate(pairs[pairs$arm=='binary',],
    Filter(function(a)a$job$arm=='binary',selected),probe)
  write.csv(gate,file.path(out,'extension-gate.csv'),row.names=FALSE)
  write.csv(data.frame(file=normalizePath(probefile),md5=unname(tools::md5sum(probefile)),
    fit=probe$source_fit,fit_md5=probe$source_fit_md5,
    spatial_flags=paste(spatial_flags(probe$spatial$diagnostics),collapse='; ')),
    file.path(out,'initialization-selection.csv'),row.names=FALSE)
  write.csv(data.frame(key=probe$job$key,nchain=probe$mcmc$nchain,nburn=probe$mcmc$nburn,
    niter=probe$mcmc$niter,starts=paste(probe$starts,collapse=', '),
    spatial_flags=paste(spatial_flags(probe$spatial$diagnostics),collapse='; ')),
    file.path(out,'initialization-summary.csv'),row.names=FALSE)
  cat('Two-stage extension gate:',if(all(gate$pass))'PASS' else 'FAIL','\n')
}
