#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args))args[1] else 'results'
out <- if(length(args)>1L)args[2] else file.path(input,'summary')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files <- unlist(lapply(strsplit(input, ',', fixed=TRUE)[[1]], list.files, pattern='-result.rds$', full.names=TRUE))
stopifnot(length(files)>0)
records <- lapply(files,readRDS)
tag <- function(a,x) { if(is.null(x))return(NULL); data.frame(scenario=readRDS(a$job$input_file)$scenario$label,replicate=readRDS(a$job$input_file)$replicate,prior=a$job$prior_name,x) }
groups <- do.call(rbind,lapply(records,function(a)tag(a,a$groups)))
elements <- do.call(rbind,lapply(records,function(a)tag(a,a$elements)))
species <- do.call(rbind,lapply(records,function(a)tag(a,a$species)))
oracles <- do.call(rbind,lapply(records,function(a)tag(a,a$known_w)))
if(nrow(oracles))oracles$bias <- oracles$estimate-oracles$truth
se <- function(x)if(length(x)>1)sd(x)/sqrt(length(x)) else NA_real_
keys <- c('scenario','prior','metric','group')
chunks <- split(groups,interaction(groups[keys],drop=TRUE))
summary <- do.call(rbind,lapply(chunks,function(x)data.frame(x[1,keys],datasets=nrow(x),truth=mean(x$truth),estimate=mean(x$estimate),bias=mean(x$bias),between_dataset_se=se(x$bias),mae=mean(x$mae),rmse=sqrt(mean(x$rmse^2)),combined_mcse=sqrt(sum(x$mcse^2))/nrow(x),max_rhat=max(x$rhat),min_ess_mean=min(x$ess_mean),max_mcse=max(x$mcse),max_chain_gap=max(x$chain_gap))))
rownames(summary)<-NULL
comparison <- list()
for(sc in unique(groups$scenario)) {
 d <- groups[groups$scenario==sc,]
 ref <- if(sc %in% c('traits_isolated','occupancy'))'var2' else 'q20'
 for(prior in setdiff(unique(d$prior),ref)) {
  paired <- merge(d[d$prior==ref,],d[d$prior==prior,],by=c('scenario','replicate','metric','group'),suffixes=c('_reference','_alternative'))
  paired$change_bias <- paired$bias_alternative-paired$bias_reference
  paired$change_mae <- paired$mae_alternative-paired$mae_reference
  paired$change_mse <- paired$rmse_alternative^2-paired$rmse_reference^2
  for(x in split(paired,interaction(paired[c('metric','group')],drop=TRUE))) comparison[[length(comparison)+1L]] <- data.frame(scenario=sc,reference=ref,alternative=prior,metric=x$metric[1],group=x$group[1],datasets=nrow(x),change_bias=mean(x$change_bias),paired_se=se(x$change_bias),change_mae=mean(x$change_mae),change_mse=mean(x$change_mse))
 }
}
comparisons <- do.call(rbind,comparison)
warning_rows <- do.call(rbind,lapply(records,function(a)data.frame(key=a$job$key,warning=if(length(a$warnings))a$warnings else NA_character_)))
for(pair in list(c('groups','dataset-groups.csv'),c('elements','parameter-elements.csv'),c('species','species-occupancy.csv'),c('oracles','known-w.csv'),c('summary','summary.csv'),c('comparisons','paired-prior-comparisons.csv'),c('warning_rows','warnings.csv'))) write.csv(get(pair[1]),file.path(out,pair[2]),row.names=FALSE)
cat(length(records),'fits;',sum(!is.na(warning_rows$warning)),'fitting warnings. Summary in',out,'\n')
options(width=200)
show <- summary[(summary$metric %in% c('B0','p','q','q_nominal') & summary$group=='all') | (summary$metric=='collection_slope' & summary$group %in% c('negative','positive')) | (summary$metric=='occupancy' & summary$group %in% c('all','low','medium','high')),]
print(show[,c('scenario','prior','metric','group','datasets','bias','between_dataset_se','combined_mcse','max_rhat')],row.names=FALSE,digits=3)
