#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
option <- function(name,default=NULL) {
  hit <- args[startsWith(args,paste0('--',name,'='))]
  if(length(hit)>1L) stop('Repeated option: ',name)
  if(!length(hit)) {if(is.null(default))stop('Missing --',name);return(default)}
  substring(hit,nchar(name)+4L)
}
repo <- normalizePath(option('repo','.'));study <- normalizePath(option('study'))
mode <- option('mode','select');stopifnot(mode %in% c('select','initial','final'))
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
selection <- read_study_selection(study,final=mode=='final')
if(mode=='select') {
  path <- file.path(study,'long-selection.csv')
  write.csv(selection$manifest,path,row.names=FALSE)
  keys <- selection$manifest$key[selection$manifest$needs_long]
  writeLines(keys,file.path(study,'long-keys.txt'))
  cat(length(keys),'of 81 fits require the prespecified longer checks.\n')
  quit(status=0)
}
out <- option('out',file.path(study,paste0('summary-',mode)))
if(dir.exists(out) && length(list.files(out))) stop('Use an empty output directory: ',out)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tables <- lapply(c('groups','elements','species','field','range'),function(t) collect_study_rows(selection$selected,t))
names(tables) <- c('groups','elements','species','field','range')
for(t in names(tables)) write.csv(tables[[t]],file.path(out,paste0(t,'.csv')),row.names=FALSE)
overall <- aggregate_scores(tables$groups)
paired <- paired_effects(tables$groups)
write.csv(overall,file.path(out,'aggregate.csv'),row.names=FALSE)
write.csv(paired,file.path(out,'paired.csv'),row.names=FALSE)
write.csv(selection$manifest,file.path(out,'selected-fits.csv'),row.names=FALSE)
range_groups <- split(tables$groups,interaction(tables$groups$arm,tables$groups$knots,
  tables$groups$metric,tables$groups$group,tables$groups$grid_index,drop=TRUE))
by_range <- do.call(rbind,lapply(range_groups,function(a) do.call(rbind,
  lapply(c('bias','mae','rmse','coverage'),function(q) data.frame(
    a[1,c('arm','knots','metric','group','grid_index')],quantity=q,n=nrow(a),
    mean=mean(a[[q]]),se=sd(a[[q]])/sqrt(nrow(a)),minimum=min(a[[q]]),maximum=max(a[[q]]))))))
write.csv(by_range,file.path(out,'by-range.csv'),row.names=FALSE)
if(mode=='final') {
  initial_groups <- collect_study_rows(selection$initial,'groups')
  community_changes <- community_sensitivity(initial_groups,tables$groups)
  write.csv(community_changes,file.path(out,'long-run-community-sensitivity.csv'),row.names=FALSE)
  original <- aggregate_scores(initial_groups)
  changes <- merge(overall,original,by=c('arm','knots','metric','group','quantity'),suffixes=c('_selected','_initial'))
  changes$change <- changes$mean_selected-changes$mean_initial
  write.csv(changes,file.path(out,'long-run-sensitivity.csv'),row.names=FALSE)
  initial_pairs <- paired_effects(initial_groups)
  pair_changes <- merge(paired,initial_pairs,by=c('comparison','group','quantity'),suffixes=c('_selected','_initial'))
  pair_changes$change <- pair_changes$mean_selected-pair_changes$mean_initial
  write.csv(pair_changes,file.path(out,'long-run-paired-sensitivity.csv'),row.names=FALSE)
}
saveRDS(selection,file.path(out,'selection.rds'))
scripts <- file.path(repo,'dev/simstudy/spatial-targeted-recheck',c('analysis.R','summarise.R'))
saveRDS(list(scripts=tools::md5sum(scripts),results=selection$manifest,
  session=sessionInfo(),time=Sys.time()),file.path(out,'provenance.rds'))
cat('Wrote',mode,'summary to',out,'\n')
