#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)==3L)
repo<-normalizePath(args[1]);study<-normalizePath(args[2]);phase<-args[3]
stopifnot(phase %in% c('select','initial','final'))
scripts<-file.path(repo,'dev/simstudy/spatial-amplitude-prior')
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
source(file.path(scripts,'metrics.R'));source(file.path(scripts,'analysis.R'));source(file.path(scripts,'robust.R'))
root<-file.path(study,'robust-v1')
read_robust<-function(path) {
  if(!file.exists(path))stop('Required robust result missing: ',path)
  a<-readRDS(path);validate_robust_spatial(a)
  stopifnot(identical(a$estimand_version,'robust-v1'),
    identical(unname(tools::md5sum(a$source_fit)),a$source_fit_md5),
    identical(unname(tools::md5sum(a$legacy_result)),a$legacy_result_md5),
    identical(unname(tools::md5sum(names(a$scoring_hashes))),unname(a$scoring_hashes)),
    identical(a$reasons,robust_reasons(a)))
  a
}
new_fit_hashes<-readRDS(file.path(study,'half_cauchy/initial/settings.rds'))$fit_hashes
stopifnot(identical(unname(tools::md5sum(names(new_fit_hashes))),unname(new_fit_hashes)))
keys<-unlist(lapply(c(4,6,8),function(g)sprintf('range%d-rep%02d-binary-k100',g,1:3)))
selection<-rows<-selected<-groups<-diagnostics<-field_diagnostics<-list()
for(key in keys) {
  igfile<-file.path(root,'inverse_gamma/baseline',paste0(key,'-initial-result.rds'))
  hcfile<-file.path(root,'half_cauchy/initial',paste0(key,'-result.rds'))
  ig<-read_robust(igfile);hc<-read_robust(hcfile);initial_ig<-ig;initial_hc<-hc
  stopifnot(ig$job$key==key,hc$job$key==key,
    identical(ig$job[c('community','grid_index','replicate','arm','knots','input_md5')],
      hc$job[c('community','grid_index','replicate','arm','knots','input_md5')]))
  validate_choice(ig,initial_ig,'inverse_gamma','initial',new_fit_hashes)
  validate_choice(hc,initial_hc,'half_cauchy','initial',new_fit_hashes)
  legacy_need<-length(ig$legacy_reasons)>0L || length(hc$legacy_reasons)>0L
  robust_need<-length(ig$reasons)>0L || length(hc$reasons)>0L
  need<-legacy_need || robust_need;chosen<-'initial'
  longig<-file.path(root,'inverse_gamma/baseline',paste0(key,'-long-result.rds'))
  if(!file.exists(longig))longig<-file.path(root,'inverse_gamma/long',paste0(key,'-result.rds'))
  longhc<-file.path(root,'half_cauchy/long',paste0(key,'-result.rds'))
  selection[[key]]<-data.frame(key=key,needs_long=need,legacy_needs_long=legacy_need,
    robust_needs_long=robust_need,initial_ig_reasons=paste(ig$reasons,collapse='; '),
    initial_hc_reasons=paste(hc$reasons,collapse='; '),
    run_ig_long=need && !file.exists(longig),run_hc_long=need && !file.exists(longhc))
  if(phase=='select')next
  if(phase=='final' && need){ig<-read_robust(longig);hc<-read_robust(longhc);igfile<-longig;hcfile<-longhc;chosen<-'long'}
  validate_choice(ig,initial_ig,'inverse_gamma',chosen,new_fit_hashes)
  validate_choice(hc,initial_hc,'half_cauchy',chosen,new_fit_hashes)
  for(prior in c('inverse_gamma','half_cauchy')) {
    a<-if(prior=='inverse_gamma')ig else hc;file<-if(prior=='inverse_gamma')igfile else hcfile
    id<-paste(key,prior,sep=':');selected[[id]]<-a
    rows[[id]]<-cbind(robust_result_row(a,prior,chosen),result_file=normalizePath(file),
      result_md5=unname(tools::md5sum(file)),fit_file=a$source_fit,fit_md5=a$source_fit_md5,
      legacy_result_file=a$legacy_result,legacy_result_md5=a$legacy_result_md5)
    prefix<-rows[[id]][c('key','community','grid_index','replicate','arm','prior','phase')]
    for(table in c('groups','diagnostics','field_diagnostics')) {
      value<-if(table=='groups')a$groups else a$spatial[[table]]
      result<-cbind(prefix[rep(1,nrow(value)),],value)
      if(table=='groups')groups[[id]]<-result else if(table=='diagnostics')diagnostics[[id]]<-result else
        field_diagnostics[[id]]<-result
    }
  }
}
out<-file.path(root,paste0('summary-binary-',phase));dir.create(out,showWarnings=FALSE)
provenance<-list(args=args,created=Sys.time(),session=sessionInfo(),
  script_hashes=tools::md5sum(file.path(scripts,c('robust.R','robust-summarise.R','rescore.R','ESTIMAND-AMENDMENT.md'))))
saveRDS(provenance,file.path(out,'provenance.rds'))
selection<-do.call(rbind,selection);write.csv(selection,file.path(out,'selection.csv'),row.names=FALSE)
cat('Amended paired longer schedules:',sum(selection$needs_long),'; pending IG:',sum(selection$run_ig_long),
  '; pending HC:',sum(selection$run_hc_long),'\n')
if(phase=='select')quit(status=0)
rows<-do.call(rbind,rows);pairs<-pair_rows(rows)
write.csv(rows,file.path(out,'fits.csv'),row.names=FALSE)
write.csv(pairs,file.path(out,'paired-communities.csv'),row.names=FALSE)
write.csv(robust_paired_summary(pairs),file.path(out,'paired-summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,groups),file.path(out,'groups.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'spatial-diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,field_diagnostics),file.path(out,'field-diagnostics.csv'),row.names=FALSE)
saveRDS(selected,file.path(out,'selected-results.rds'))
if(phase=='final') {
  key<-'range6-rep01-binary-k100'
  initial<-read_robust(file.path(root,'half_cauchy/initial',paste0(key,'-result.rds')))
  probefile<-file.path(root,'half_cauchy/starts-initial',paste0(key,'-result.rds'))
  probe<-read_robust(probefile)
  validate_choice(probe,initial,'half_cauchy','starts-initial',new_fit_hashes)
  if(length(probe$legacy_reasons)>0L || any(robust_flag_rows(probe$spatial$diagnostics)) ||
     any(robust_flag_rows(probe$spatial$field_diagnostics))) {
    probefile<-file.path(root,'half_cauchy/starts-long',paste0(key,'-result.rds'))
    probe<-read_robust(probefile);validate_choice(probe,initial,'half_cauchy','starts-long',new_fit_hashes)
  }
  gate<-robust_extension_gate(pairs,selected,probe)
  write.csv(gate,file.path(out,'extension-gate.csv'),row.names=FALSE)
  flags<-paste(robust_reasons(probe),collapse='; ')
  write.csv(data.frame(file=normalizePath(probefile),md5=unname(tools::md5sum(probefile)),
    fit=probe$source_fit,fit_md5=probe$source_fit_md5,spatial_flags=flags),
    file.path(out,'initialization-selection.csv'),row.names=FALSE)
  write.csv(data.frame(key=probe$job$key,nchain=probe$mcmc$nchain,nburn=probe$mcmc$nburn,
    niter=probe$mcmc$niter,starts=paste(probe$starts,collapse=', '),spatial_flags=flags),
    file.path(out,'initialization-summary.csv'),row.names=FALSE)
  write.csv(probe$spatial$diagnostics,file.path(out,'initialization-diagnostics.csv'),row.names=FALSE)
  write.csv(probe$spatial$field_diagnostics,file.path(out,'initialization-field-diagnostics.csv'),row.names=FALSE)
  cat('Amended gate before numerical audit:',if(all(gate$pass))'PASS' else 'FAIL','\n')
}
