#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
rows<-list()
for(K in c(3L,6L))for(prior in c('q20','q9','p32')) {
 folder<-if(K==3L)'main-results' else 'operational-k6';key<-sprintf('qfar_K%d-%s-02',K,prior)
 a<-readRDS(file.path(root,folder,paste0(key,'-fit.rds')));r<-readRDS(file.path(root,folder,paste0(key,'-result.rds')))
 f<-a$fit;sp<-match('OTU_3',f$infos$speciesNames);primer<-match('2',as.character(f$infos$primerNames));stopifnot(!is.na(sp),!is.na(primer))
 p<-f$results_output$p_output[primer,sp,,];q<-f$results_output$q_output[primer,sp,,];gap<-p-q;index<-primer+(sp-1L)*length(f$infos$primerNames)
 tp<-r$elements$truth[r$elements$metric=='p' & r$elements$element==index];tq<-r$elements$truth[r$elements$metric=='q' & r$elements$element==index]
 rows[[length(rows)+1L]]<-data.frame(key=key,replicates_per_primer=K,prior=prior,species='OTU_3',primer='2',true_p=tp,true_q=tq,estimated_p=mean(p),estimated_q=mean(q),probability_p_le_q=mean(gap<=0),chain1_probability=mean(gap[,1]<=0),chain2_probability=mean(gap[,2]<=0),chain1_order_transitions=sum(diff(gap[,1]<=0)!=0),chain2_order_transitions=sum(diff(gap[,2]<=0)!=0),gap_rhat=posterior::rhat(gap),gap_ess_mean=posterior::ess_mean(gap),gap_mcse=posterior::mcse_mean(gap))
}
x<-do.call(rbind,rows);write.csv(x,file.path(root,'main-results/summary/detection-ordering-case.csv'),row.names=FALSE)
options(width=170);print(x,row.names=FALSE,digits=5)
