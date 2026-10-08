#!/usr/bin/env Rscript
args<-commandArgs(TRUE);repo<-normalizePath(if(length(args))args[1]else '.')
out<-file.path(repo,'dev/simstudy/theta0-prior-mean/results')
suppressPackageStartupMessages(library(ggplot2))
s<-read.csv(file.path(out,'importance-summary.csv'));s<-subset(s,contamination=='low' & theta_b==20)
p<-read.csv(file.path(out,'paired-summary.csv'));p<-subset(p,contamination=='low' & scope=='all' & metric=='width_percent' & grepl('theta b=20',contrast,fixed=TRUE))
labels<-c(baseline='100 sites, 2 field samples',field4='100 sites, 4 field samples',sites300='300 sites, 2 field samples')
s$design<-factor(s$design,levels=names(labels),labels=labels);p$design<-factor(p$design,levels=names(labels),labels=labels)
z<-rbind(data.frame(design=s$design,arm=paste0('Collection slope mean ',s$collection_mean),metric='Mean theta0 interval width',value=s$width),
  data.frame(design=s$design,arm=paste0('Collection slope mean ',s$collection_mean),metric='Mean signed theta0 error',value=s$bias))
figure<-ggplot(z,aes(value,design,color=arm))+geom_vline(xintercept=0,color='grey65',linewidth=.4)+
  geom_point(data=subset(z,arm=='Collection slope mean 0'),position=position_nudge(y=-.08),size=3)+
  geom_point(data=subset(z,arm=='Collection slope mean 1'),position=position_nudge(y=.08),size=3)+
  facet_wrap(~metric,ncol=1,scales='free_x')+
  scale_color_manual(values=c('#24699b','#bb5a20'))+theme_minimal(base_size=12)+
  theme(legend.position='bottom',panel.grid.minor=element_blank())+labs(y=NULL,x='Probability units',color=NULL,
    title='Collection-prior mean has little effect on theta0 in these fits',
    subtitle='10 low-contamination communities per design; theta0 prior fixed at Beta(1,20)',
    caption='Full-joint posterior reweighting. All fits retained; 17/30 have source diagnostic flags.\nThese targeted comparisons do not establish universal interval calibration.')
ggsave(file.path(out,'theta0-mean-sensitivity.png'),figure,width=9,height=6.5,dpi=180,bg='white')
change<-ggplot(p,aes(design,change))+geom_hline(yintercept=0,color='grey65')+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,color='#24699b')+geom_point(size=3,color='#24699b')+
  coord_flip()+theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank())+
  labs(x=NULL,y='Paired change in interval width (%)',title='Restoring mean zero does not reproduce the historical 25% widening',
    subtitle='Mean 1 to mean 0; theta0 prior fixed at Beta(1,20)',
    caption='95% t intervals across communities; source sampling uncertainty is additional.\nThe historical experiments also differed in other settings and lack frozen loaded-source provenance.')
ggsave(file.path(out,'theta0-paired-width-change.png'),change,width=9,height=4.7,dpi=180,bg='white')
cat('Two figures saved.\n')


# Historical theta-prior contrast at the current collection mean; keep old figures.
all_pairs<-read.csv(file.path(out,'paired-summary.csv'))
h<-subset(all_pairs,contamination=='low' & scope=='all' & contrast=='theta b 30 to 20, collection mean=0' &
  metric %in% c('bias','width_percent','covered'))
if(nrow(h)) {
  h$design<-factor(h$design,levels=names(labels),labels=labels)
  probability<-h$metric %in% c('bias','covered')
  for(v in c('change','lower','upper'))h[[v]][probability]<-100*h[[v]][probability]
  h$metric<-factor(h$metric,levels=c('bias','width_percent','covered'),
    labels=c('Change in mean error (percentage points)','Change in interval width (%)','Change in coverage (percentage points)'))
  history_plot<-ggplot(h,aes(change,design))+geom_vline(xintercept=0,color='grey65')+
    geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.15,color='#24699b')+
    geom_point(size=3,color='#24699b')+facet_wrap(~metric,ncol=1,scales='free_x')+
    theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank())+
    labs(x=NULL,y=NULL,title='Effect of the historical theta0-prior change',
      subtitle='Beta(1,30) to Beta(1,20); collection slope mean fixed at zero',
      caption='10 low-contamination communities per design; 95% t intervals across communities.\nSource posterior sampling uncertainty is additional. This does not reconstruct the old experiments.')
  ggsave(file.path(out,'theta0-historical-prior-change.png'),history_plot,width=9,height=7.5,dpi=180,bg='white')
  cat('Historical-prior figure saved.\n')
}
