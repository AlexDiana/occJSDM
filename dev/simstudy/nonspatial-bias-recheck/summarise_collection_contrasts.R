#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else '.')
x <- read.csv(file.path(root,'q-audit/completed-fits/collection_groups.csv'))
out <- file.path(root,'main-results/summary')
se <- function(x) if(length(x)>1)sd(x)/sqrt(length(x)) else NA_real_
keys <- c('scenario','prior','metric','slope_sign')
a <- do.call(rbind,lapply(split(x,interaction(x[keys],drop=TRUE)),function(z) data.frame(z[1,keys],datasets=nrow(z),species=sum(z$n_species),truth=mean(z$truth),estimate=mean(z$estimate),bias=mean(z$bias),between_dataset_se=se(z$bias),mae=mean(z$mae),rmse=sqrt(mean(z$rmse^2)),combined_mcse=sqrt(sum(z$mcse^2))/nrow(z),max_rhat=max(z$rhat,na.rm=TRUE))))
write.csv(a,file.path(out,'collection-contrast-summary.csv'),row.names=FALSE)
paired <- list()
for(sc in unique(x$scenario)) {
 z<-x[x$scenario==sc,];ref<-if(sc %in% c('occupancy','traits_isolated'))'var2' else 'q20'
 for(pr in setdiff(unique(z$prior),ref)) {
  p<-merge(z[z$prior==ref,],z[z$prior==pr,],by=c('scenario','replicate','metric','slope_sign'),suffixes=c('_ref','_alt'))
  for(g in split(p,interaction(p[c('metric','slope_sign')],drop=TRUE))) paired[[length(paired)+1L]]<-data.frame(scenario=sc,reference=ref,alternative=pr,metric=g$metric[1],slope_sign=g$slope_sign[1],datasets=nrow(g),change_bias=mean(g$bias_alt-g$bias_ref),paired_se=se(g$bias_alt-g$bias_ref),change_mae=mean(g$mae_alt-g$mae_ref))
 }
}
write.csv(do.call(rbind,paired),file.path(out,'collection-contrast-paired-comparison.csv'),row.names=FALSE)
cat('Summarised collection contrasts for',length(unique(x$key)),'fits.\n')
