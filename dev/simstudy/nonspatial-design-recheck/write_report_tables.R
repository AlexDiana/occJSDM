#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/nonspatial-design-recheck-20260919')
s<-read.csv(file.path(here,'summary-long/summary.csv'));p<-read.csv(file.path(here,'summary-long/paired-comparisons.csv'));m<-read.csv(file.path(here,'summary-long/fit-manifest.csv'))
stopifnot(nrow(m)==60L,sum(m$longer_fit)==8L)
labels<-c(baseline='100 sites, 2 field samples',field4='100 sites, 4 field samples',sites300='300 sites, 2 field samples',knownU='Hidden site conditions supplied')
get<-function(arm,sc,group='all',metric='occupancy',scope=if(metric=='occupancy')'original_sites' else if(metric=='collection_probability')'original_samples' else 'all_fitted_sites'){
 ix<-s$arm==arm & s$scenario==sc & s$metric==metric & s$group==group & s$scope==scope;z<-s[ix,];stopifnot(nrow(z)==1L);z
}
lines<-c('## How much occupancy error remains?','',
'All numbers in this section are measured simulation results after substituting the eight longer fits. Errors are in **percentage points**. Each community receives equal weight.','',
'| Design | Average absolute error: low contamination | Average absolute error: high contamination |',
'|---|---:|---:|')
for(arm in names(labels))lines<-c(lines,sprintf('| %s | %.1f | %.1f |',labels[arm],100*get(arm,'qnear_K6')$mae,100*get(arm,'qfar_K6')$mae))
lines<-c(lines,'',
'![Average absolute occupancy errors at the same original 100 sites](nonspatial-design-recheck/results/occupancy-absolute-error-by-design.png)','',
'## How consistent are the improvements?','',
'For each community, subtract the new absolute error from its baseline absolute error. A positive reduction means the new design did better. The figure summarizes those ten paired changes. Its bars describe uncertainty in the average improvement across communities; they do not describe uncertainty around a particular species-at-site probability.','',
'![Paired reductions in average absolute occupancy errors](nonspatial-design-recheck/results/occupancy-error-reduction-by-design.png)','',
'| Change from baseline | Reduction under low contamination (95% interval) | Reduction under high contamination (95% interval) |',
'|---|---:|---:|')
for(arm in names(labels)[-1]){
 vals<-vapply(c('qnear_K6','qfar_K6'),function(sc){z<-p[p$arm==arm & p$scenario==sc & p$metric=='occupancy' & p$group=='all',];stopifnot(nrow(z)==1L);sprintf('%.2f (%.2f to %.2f)',-100*z$change_mae,100*(-z$change_mae-qt(.975,9)*z$change_mae_se),100*(-z$change_mae+qt(.975,9)*z$change_mae_se))},character(1))
 lines<-c(lines,sprintf('| %s | %s | %s |',labels[arm],vals[1],vals[2]))
}
lines<-c(lines,'','These are paired t intervals based on ten communities. They do not include every source of ecological or numerical uncertainty. The full [paired comparisons](nonspatial-design-recheck/results/paired-comparisons.csv) retain standard errors and the number of communities that improved.','',
'## Are low probabilities still too high, and high probabilities too low?','',
'The table below gives **average signed error**. Values closer to zero are better. Positive values mean the model estimates occupancy too high; negative values mean it estimates too low.','',
'| Contamination | Design | True probability below 20% | True probability 20% to 80% | True probability above 80% |',
'|---|---|---:|---:|---:|')
for(sc in c('qnear_K6','qfar_K6'))for(arm in names(labels)){
 vals<-vapply(c('low','medium','high'),function(g)sprintf('%+.1f',100*get(arm,sc,g)$bias),character(1))
 lines<-c(lines,sprintf('| %s | %s | %s | %s | %s |',if(sc=='qnear_K6')'Low' else 'High',labels[arm],vals[1],vals[2],vals[3]))
}
lines<-c(lines,'','![Signed occupancy errors in the three true-probability groups](nonspatial-design-recheck/results/occupancy-bias-by-design.png)','',
'The middle group still illustrates cancellation: small signed averages coexist with much larger individual errors. Its mean absolute errors are:','',
'| Design | Middle-group absolute error: low contamination | Middle-group absolute error: high contamination |',
'|---|---:|---:|')
for(arm in names(labels))lines<-c(lines,sprintf('| %s | %.1f | %.1f |',labels[arm],100*get(arm,'qnear_K6','medium')$mae,100*get(arm,'qfar_K6','medium')$mae))
writeLines(lines,file.path(here,'report-tables.md'))
# Compact secondary results for the report writer; the complete CSVs are retained.
sec<-subset(s,metric %in% c('p','q','collection_probability'));sec$mae_pp<-100*sec$mae;sec$bias_pp<-100*sec$bias
write.csv(sec,file.path(here,'summary-long/secondary-results.csv'),row.names=FALSE)
cat('Wrote report tables from the eight-longer-fit comparison.\n')
d<-read.csv(file.path(here,'summary-long/dataset-groups.csv'));pairs<-list()
for(sc in c('qnear_K6','qfar_K6')){
 a<-subset(d,scenario==sc & arm=='field4' & scope=='original_sites' & metric=='occupancy' & group=='all');b<-subset(d,scenario==sc & arm=='sites300' & scope=='original_sites' & metric=='occupancy' & group=='all')
 a<-a[order(a$replicate),];b<-b[order(b$replicate),];stopifnot(identical(a$replicate,b$replicate));diff<-a$mae-b$mae
 pairs[[sc]]<-data.frame(scenario=sc,datasets=length(diff),extra_mae_reduction_sites300_over_field4=mean(diff),se=sd(diff)/sqrt(length(diff)),lower=mean(diff)-qt(.975,9)*sd(diff)/sqrt(10),upper=mean(diff)+qt(.975,9)*sd(diff)/sqrt(10),sites300_better_datasets=sum(diff>0))
}
write.csv(do.call(rbind,pairs),file.path(here,'summary-long/design-alternative-comparison.csv'),row.names=FALSE)
