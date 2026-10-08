#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
repo<-normalizePath(if(length(args))args[1] else '.')
out<-file.path(repo,'dev/simstudy/interval-calibration/results')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(repo,'dev/simstudy/interval-calibration/helpers.R'))
source(file.path(repo,'tests/testthat/helper-simstudy.R'))
history<-file.path(repo,'dev/simstudy/results/simstudy-20260821-135619.rds')
x<-readRDS(history)
old<-subset(x$rows,block=='q');z<-old
z$nominal_truth<-z$truth
z$nominal_covered<-z$covered
retention<-pnorm(log(1.5),mean=1.5,sd=1,lower.tail=FALSE)
z$truth<-z$truth*retention
corrected<-write_interval_tables(z,'historical-q',out)
stopifnot(nrow(z)==2400L,all(table(z$scenario)==600L),x$R==30L)
write.csv(aggregate(cbind(nominal_covered,covered)~scenario,corrected$rows,mean),file.path(out,'historical-q-correction.csv'),row.names=FALSE)

# Selected fits include all chains of the prespecified longer runs.
archive<-file.path(repo,'dev/simstudy/results/pr11-current-20260927')
selection<-readRDS(file.path(archive,'selection.rds'))
selection<-subset(selection,family=='design' & arm!='knownU')
stopifnot(nrow(selection)==60L)
rows<-list();hashes<-list()
for(i in seq_len(nrow(selection))) {
  job<-selection[i,];input<-readRDS(job$input_file)
  fitfile<-sub('-result.rds$','-fit.rds',job$result_file)
  saved<-readRDS(fitfile);fit<-saved$fit
  stopifnot(unname(tools::md5sum(job$input_file))==job$input_md5,fit$infos$ps==0L)
  blocks<-simstudy_param_blocks(fit,input$sim,input$truth)
  scenario<-paste(job$scenario,job$arm,sep=':')
  bt<-blocks$beta_theta
  slope<-interval_rows(bt$post[2,,,,drop=FALSE],bt$truth[2,],'collection_slope',scenario,job$replicate)
  slope$group<-'all'
  signs<-slope;signs$group<-c('negative','zero','positive')[sign(signs$truth)+2]
  intercept<-interval_rows(bt$post[1,,,,drop=FALSE],bt$truth[1,],'collection_intercept',scenario,job$replicate)
  q<-interval_rows(blocks$q$post,blocks$q$truth,'q',scenario,job$replicate)
  intercept$group<-q$group<-'all'
  rows[[i]]<-rbind(slope,signs,intercept,q)
  hashes[[i]]<-data.frame(key=job$key,schedule=job$schedule,input=job$input_file,input_md5=job$input_md5,
    fit=fitfile,fit_md5=unname(tools::md5sum(fitfile)),warnings=length(saved$warnings),
    original_max_rhat=job$max_element_rhat)
  cat(i,'/60',job$key,'scored\n');flush.console()
  rm(saved,fit);gc(FALSE)
}
write_interval_tables(do.call(rbind,rows),'selected-nonspatial',out)
write.csv(do.call(rbind,hashes),file.path(out,'selected-fit-provenance.csv'),row.names=FALSE)
source_paths<-c(history,file.path(archive,'selection.rds'),file.path(repo,'tests/testthat/helper-simstudy.R'),
  file.path(repo,'dev/simstudy/interval-calibration',c('helpers.R','rescore.R')))
write.csv(data.frame(file=source_paths,md5=unname(tools::md5sum(source_paths))),file.path(out,'rescore-source-hashes.csv'),row.names=FALSE)
cat('Historical q and 60 selected nonspatial fits rescored.\n')
