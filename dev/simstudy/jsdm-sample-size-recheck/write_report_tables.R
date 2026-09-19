args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/jsdm-sample-size-20260919');out<-file.path(here,'summary')
s<-read.csv(file.path(out,'summary.csv'));p<-read.csv(file.path(out,'paired-summary.csv'));g<-read.csv(file.path(out,'dataset-groups.csv'));f<-read.csv(file.path(out,'fit-manifest.csv'))
lines<-c('| Sites used for fitting | Average absolute error | Average signed error | All fitted sites: absolute error |','|---:|---:|---:|---:|')
for(n in c(100,300,1000)) {
 a<-s[s$n==n&s$scope=='original100'&s$group=='all',];b<-s[s$n==n&s$scope=='allsites'&s$group=='all',]
 lines<-c(lines,sprintf('| %s | %.1f points | %+.1f points | %.1f points |',format(n,big.mark=',',trim=TRUE),100*a$mae,100*a$bias,100*b$mae))
}
lines<-c(lines,'','| Added sites | Reduction in absolute error | 95% interval across ten communities | Communities that improved |','|---|---:|---:|---:|')
for(k in list(c(100,300),c(100,1000),c(300,1000))) {
 a<-p[p$from==k[1]&p$to==k[2]&p$group=='all',]
 lines<-c(lines,sprintf('| %s to %s | %.2f points | %.2f to %.2f points | %d of 10 |',format(k[1],big.mark=',',trim=TRUE),format(k[2],big.mark=',',trim=TRUE),100*a$reduction_mae,100*a$lower,100*a$upper,a$improved))
}
lines<-c(lines,'','| True probability | Sites used for fitting | Average true probability | Average estimate | Signed error | Absolute error |','|---|---:|---:|---:|---:|---:|')
for(b in c('low','middle','high'))for(n in c(100,300,1000)) {
 a<-s[s$n==n&s$scope=='original100'&s$group==b,]
 lines<-c(lines,sprintf('| %s | %s | %.1f%% | %.1f%% | %+.1f points | %.1f points |',switch(b,low='Below 20%',middle='20-80%',high='Above 80%'),format(n,big.mark=',',trim=TRUE),100*a$truth,100*a$estimate,100*a$bias,100*a$mae))
}
lines<-c(lines,'',sprintf('Diagnostics: maximum probability-group Rhat %.6f; maximum element Rhat %.6f; minimum element ESS %.1f; maximum combined mean-probability MCSE %.5f percentage points; selected-fit warnings %d.',max(f$max_group_rhat),max(f$max_element_rhat),min(f$min_element_ess),100*max(s$combined_mean_mcse),sum(f$warnings)))
writeLines(lines,file.path(here,'report-tables.md'));cat(paste(lines,collapse='\n'),'\n')
