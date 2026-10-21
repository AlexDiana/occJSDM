#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
repo<-normalizePath(if(length(args))args[1]else'.')
out<-file.path(repo,'dev/simstudy/interval-calibration/results')
suppressPackageStartupMessages(library(ggplot2))
historical<-read.csv(file.path(out,'historical-q-summary.csv'))
historical$contamination<-ifelse(grepl('qnear',historical$scenario),'Low contamination','High contamination')
historical$PCRs<-ifelse(grepl('K30',historical$scenario),30,3)
p1<-ggplot(historical,aes(PCRs,coverage,color=contamination))+geom_hline(yintercept=.95,linetype=2,color='grey40')+
  geom_errorbar(aes(ymin=coverage_mc_lower,ymax=coverage_mc_upper),width=.8)+geom_line()+geom_point(size=2.3)+
  scale_x_continuous(breaks=c(3,30))+coord_cartesian(ylim=c(.88,1))+
  labs(title='Historical q after correcting truth',subtitle='30 communities per cell; K30 is diagnostic only',
    x='PCR replicates per primer',y='Coverage of 95% intervals',color=NULL)+theme_bw(base_size=11)+theme(legend.position='bottom')
current<-read.csv(file.path(out,'selected-nonspatial-summary.csv'))
current<-subset(current,block=='collection_slope'&group=='all')
current$contamination<-ifelse(grepl('qnear',current$scenario),'Low contamination','High contamination')
current$design<-factor(sub('^.*:','',current$scenario),c('baseline','field4','sites300'),
  c('100 sites, 2 samples','100 sites, 4 samples','300 sites, 2 samples'))
p2<-ggplot(current,aes(design,coverage,color=contamination,group=contamination))+
  geom_hline(yintercept=.95,linetype=2,color='grey40')+
  geom_errorbar(aes(ymin=coverage_mc_lower,ymax=coverage_mc_upper),position=position_dodge(.25),width=.12)+
  geom_point(position=position_dodge(.25),size=2.3)+coord_cartesian(ylim=c(.88,1))+
  labs(title='Collection slopes in selected nonspatial fits',subtitle='10 communities per design; convergence flags retained',
    x=NULL,y='Coverage of 95% intervals',color=NULL)+theme_bw(base_size=11)+theme(legend.position='bottom',axis.text.x=element_text(angle=15,hjust=1))
continuous<-subset(read.csv(file.path(out,'continuous-selected-summary.csv')),block=='B0')
continuous$design<-factor(continuous$scenario,c('d0','d2'),c('No hidden factors','Two hidden factors'))
p3<-ggplot(continuous,aes(design,coverage))+geom_hline(yintercept=.95,linetype=2,color='grey40')+
  geom_errorbar(aes(ymin=coverage_mc_lower,ymax=coverage_mc_upper),width=.12,color='#245881')+
  geom_point(size=2.7,color='#245881')+
  coord_cartesian(ylim=c(max(0,min(c(.88,continuous$coverage,continuous$coverage_mc_lower),na.rm=TRUE)-.01),1))+
  labs(title='Continuous B0 without spatial effects',subtitle='100 paired communities; diagnostic-selected longer chains',x=NULL,y='Coverage of 95% intervals')+theme_bw(base_size=11)
p4<-ggplot(current,aes(design,width,color=contamination,group=contamination))+geom_line()+geom_point(size=2.3)+
  labs(title='Collection slope intervals become narrower',subtitle='The means above assess coverage separately',x=NULL,y='Mean interval width',color=NULL)+
  theme_bw(base_size=11)+theme(legend.position='bottom',axis.text.x=element_text(angle=15,hjust=1))
fig<-gridExtra::arrangeGrob(p1,p2,p3,p4,ncol=2)
ggsave(file.path(out,'interval-calibration.png'),fig,width=13,height=8.5,dpi=180,bg='white')
