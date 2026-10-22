#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/jsdm-sample-size-20260919');out<-file.path(here,'summary')
s<-read.csv(file.path(out,'summary.csv'));d<-read.csv(file.path(out,'dataset-groups.csv'));cal<-read.csv(file.path(out,'calibration.csv'))
suppressPackageStartupMessages(library(ggplot2))
colours<-c('100'='#a75126','300'='#137f91','1000'='#495db3')
base_theme<-theme_minimal(base_size=12)+theme(plot.title=element_text(size=19,face='bold',colour='#182e3b'),plot.subtitle=element_text(size=11,margin=margin(b=12)),strip.text=element_text(face='bold',size=12),panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),plot.caption=element_text(hjust=0,size=10,lineheight=1.15,margin=margin(t=12)),plot.margin=margin(15,20,15,15))
z<-s[s$scope=='original100' & s$group=='all',];raw<-d[d$scope=='original100' & d$group=='all',]
z$size<-factor(z$n,levels=c(100,300,1000));raw$size<-factor(raw$n,levels=c(100,300,1000))
g<-ggplot()+geom_line(data=raw,aes(n,100*mae,group=replicate),colour='#c3cbd0',linewidth=.7)+geom_point(data=raw,aes(n,100*mae),colour='#a3adb5',size=1.8)+
 geom_errorbar(data=z,aes(n,ymin=100*(mae-qt(.975,9)*mae_se),ymax=100*(mae+qt(.975,9)*mae_se),colour=size),width=20,linewidth=.9)+geom_point(data=z,aes(n,100*mae,colour=size),size=4)+
 geom_text(data=z,aes(n,100*mae,label=sprintf('%.1f points',100*mae)),nudge_x=35,hjust=0,size=4)+scale_colour_manual(values=colours)+scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.08)))+scale_x_continuous(breaks=c(100,300,1000),labels=c('100','300','1,000'),limits=c(0,1200),expand=expansion(mult=0))+
 labs(title='Does adding sites improve the JSDM alone?',subtitle='True presence/absence supplied. Errors assessed at the same original 100 sites in every fit.',x='Number of sites used to fit the model (10 species throughout)',y='Average absolute probability error\n(percentage points; smaller is better)',caption=paste('Each grey line follows one simulated community. Coloured points average ten communities.',
 'Bars are 95% t intervals for the mean across communities, not uncertainty for an individual occupancy estimate.',
 'Errors compare estimated probabilities with simulated probabilities. Knowing a presence or absence does not reveal its probability.',sep='\n'))+base_theme+theme(legend.position='none')
ggsave(file.path(out,'absolute-error-by-sites.png'),g,width=11,height=7,dpi=170,bg='white')
cal$size<-factor(cal$n,levels=c(100,300,1000))
g<-ggplot(cal,aes(100*truth,100*estimate,colour=size))+geom_abline(slope=1,intercept=0,linetype=2,colour='#707e88')+geom_line(linewidth=.85)+geom_point(size=2.5)+
 scale_colour_manual(values=colours,labels=c('100 sites','300 sites','1,000 sites'))+scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20),expand=expansion(mult=0))+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),expand=expansion(mult=0))+coord_equal()+
 labs(title='Are low and high probabilities pulled towards the middle?',subtitle='The same original 100 sites, grouped into ten bands of true probability.',x='True occupancy probability (%)',y='Average estimated occupancy probability (%)',caption=paste('The dashed line represents perfect average recovery. Above it: overestimation. Below it: underestimation.',
 'Each point first averages within a community and probability band, then across the ten communities.',
 'This plot shows the direction of error. Individual errors can be larger because opposite errors cancel in an average.',sep='\n'))+base_theme
 ggsave(file.path(out,'true-versus-estimated.png'),g,width=10,height=9,dpi=170,bg='white')
z<-s[s$scope=='original100' & s$group!='all',];z$band<-factor(z$group,levels=c('low','middle','high'),labels=c('Below 20%','20-80%','Above 80%'));z$size<-factor(z$n,levels=c(100,300,1000))
both<-rbind(transform(z,measure='Direction of error',value=bias,err=bias_se),transform(z,measure='Size of error',value=mae,err=mae_se))
both$measure<-factor(both$measure,levels=c('Direction of error','Size of error'))
g<-ggplot(both,aes(band,100*value,colour=size))+geom_hline(yintercept=0,colour='#7b8b94',linetype=2)+geom_errorbar(aes(ymin=100*(value-qt(.975,9)*err),ymax=100*(value+qt(.975,9)*err)),position=position_dodge(width=.65),width=.12)+geom_point(position=position_dodge(width=.65),size=3)+facet_wrap(~measure,scales='free_y')+scale_colour_manual(values=colours,labels=c('100 sites','300 sites','1,000 sites'))+
 labs(title='Separate the direction of error from its size',subtitle='Ten paired communities; identical original sites used for every comparison.',x='True occupancy probability',y='Probability error (percentage points)',caption=paste('Left: signed error. Positive means overestimated; negative means underestimated. Opposite errors can cancel.',
 'Right: absolute error. Count the size of each error before averaging, so opposite errors cannot cancel.',
 'Points average ten communities; bars show 95% t intervals for those means.',sep='\n'))+base_theme
 ggsave(file.path(out,'errors-by-probability-band.png'),g,width=12,height=7,dpi=170,bg='white')
cat('Rendered three figures from the completed summaries.\n')
