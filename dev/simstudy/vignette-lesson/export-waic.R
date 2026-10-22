# Rscript dev/simstudy/vignette-lesson/export-waic.R ARCHIVE
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
a <- normalizePath(args[1], mustWork = TRUE)
small <- readRDS(file.path(a,'waic-1000/validation.rds'))
large <- readRDS(file.path(a,'waic-4000/validation.rds'))
summaries <- do.call(rbind,lapply(c('1000','4000'),function(n) read.csv(file.path(a,paste0('waic-',n),'summary.csv'))))
comparisons <- do.call(rbind,lapply(list(small,large),function(x) {
 data.frame(draws=nrow(x$results$one_factor$score$log_lik), difference_one_minus_two=x$comparison$difference, paired_site_SE=x$comparison$SE,
 approximate_difference_MCSE=sqrt(sum(vapply(x$results,function(z) z$waic_mcse^2,numeric(1)))))
}))
stopifnot(nrow(summaries)==4,
 identical(small$results$one_factor$score$observations,large$results$one_factor$score$observations))
source_paths<-c('R/site-waic.R','src/site_waic.cpp','dev/simstudy/vignette-lesson/validate-waic.R')
stopifnot(identical(large$source_md5,tools::md5sum(source_paths)),identical(small$source_md5,large$source_md5))
bundle <- list(date='2026-10-07',source_commit='a56f5548e68e09d244513d09fcc04ad52f0caa19',
 summaries=summaries,comparisons=comparisons,
 pointwise=lapply(large$results,function(x)x$score$pointwise),
 source_md5=large$source_md5,
 exporter_md5=tools::md5sum("dev/simstudy/vignette-lesson/export-waic.R"),
 validation_md5=tools::md5sum(file.path(a,c('waic-1000/validation.rds','waic-4000/validation.rds'))),
 survey=list(sites=100L,species=10L,samples_per_site=3L,primers=2L,pcrs_per_primer=6L,threshold=1),
 retained_draws_per_fit=24000L,
 conclusion='Numerical checks pass; WAIC reliability remains unvalidated. Do not select a factor count from this example.')
saveRDS(bundle,'vignettes/teaching-data/site-waic-lesson.rds',compress='xz')
p <- 'dev/simstudy/vignette-lesson/post-pr14-validation'
write.csv(summaries,file.path(p,'waic-summary.csv'),row.names=FALSE)
write.csv(comparisons,file.path(p,'waic-comparisons.csv'),row.names=FALSE)
write.csv(data.frame(path=names(large$source_md5),md5=unname(large$source_md5)),file.path(p,'source-hashes.csv'),row.names=FALSE)
for(n in c('1000','4000'))file.copy(file.path(a,paste0('waic-',n),'summary.csv'),file.path(p,paste0('waic-',n,'.csv')),overwrite=TRUE)
print(comparisons)
cat('Strict total-score changes:',max(abs(summaries$WAIC-summaries$stricter_WAIC)),'\n')
