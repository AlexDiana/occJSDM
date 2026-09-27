#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/nonspatial-design-recheck-20260919')
out<-file.path(here,if(length(args)>1L && args[2]=='long')'summary-long' else 'summary');s<-read.csv(file.path(out,'summary.csv'));p<-read.csv(file.path(out,'paired-comparisons.csv'))
suppressPackageStartupMessages(library(ggplot2))
colours<-c(baseline='#3f4d58',field4='#c67313',sites300='#2377af',knownU='#8b509f')
labels<-c(baseline='100 sites, 2 field samples',field4='100 sites, 4 field samples',sites300='300 sites, 2 field samples',knownU='100 sites, 2 samples; hidden site conditions known')
prepare<-function(x){x$arm<-factor(x$arm,levels=names(labels));x$contamination<-factor(x$scenario,levels=c('qnear_K6','qfar_K6'),labels=c('Low contamination','High contamination'));x}
s<-prepare(s);p<-prepare(p)
base_theme<-theme_minimal(base_size=12)+theme(plot.title=element_text(size=18,face='bold',colour='#142839'),plot.subtitle=element_text(size=11,margin=margin(b=14)),strip.text=element_text(face='bold',size=12),panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),plot.caption=element_text(hjust=0,size=9,lineheight=1.2,margin=margin(t=12)),plot.margin=margin(15,20,15,15))
z<-subset(s,scope=='original_sites' & metric=='occupancy' & group %in% c('low','medium','high'))
z$group<-factor(z$group,levels=c('low','medium','high'),labels=c('Below 20%','20-80%','Above 80%'))
z$error<-100*z$bias;z$lower<-100*(z$bias-qt(.975,9)*z$bias_se);z$upper<-100*(z$bias+qt(.975,9)*z$bias_se)
g<-ggplot(z,aes(group,error,colour=arm))+geom_hline(yintercept=0,colour='#778899',linetype=2)+
 geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=.65),width=.12,linewidth=.5)+geom_point(position=position_dodge(width=.65),size=2.8)+
 facet_wrap(~contamination,nrow=1)+scale_colour_manual(values=colours,labels=labels)+
 labs(title='Do the designs reduce systematic occupancy errors?',subtitle='All designs assessed at the same original 100 sites. Six PCR replicates per primer; default priors.',x='True occupancy probability',y='Average signed occupancy error\n(percentage points)',caption=paste('Above zero: occupancy overestimated. Below zero: underestimated. Points average ten paired simulated communities.',
 'Bars are 95% t intervals for the between-dataset mean; they are not uncertainty intervals for individual occupancy estimates.',
 'The control is given the true hidden site conditions (U), but still learns species responses, occupancy and detection. This is an idealized diagnostic.',sep='\n'))+
 base_theme+guides(colour=guide_legend(nrow=2,byrow=TRUE))
ggsave(file.path(out,'occupancy-bias-by-design.png'),g,width=12,height=7,dpi=160,bg='white')
z<-subset(s,scope=='original_sites' & metric=='occupancy' & group=='all')
z$error<-100*z$mae;z$lower<-pmax(0,100*(z$mae-qt(.975,9)*z$mae_se));z$upper<-100*(z$mae+qt(.975,9)*z$mae_se)
z$design<-factor(as.character(z$arm),levels=rev(names(labels)),labels=rev(c('100 sites, 2 samples\nBaseline','100 sites, 4 samples\nMore field samples','300 sites, 2 samples\nMore sites','100 sites, 2 samples\nHidden site conditions known')))
g<-ggplot(z,aes(error,design,colour=arm))+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.15,linewidth=.6)+geom_point(size=3)+
 geom_text(aes(label=sprintf('%.1f points',error)),nudge_y=.2,size=3.4,show.legend=FALSE)+facet_wrap(~contamination,nrow=1)+
 scale_colour_manual(values=colours,labels=labels)+scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.12)))+
 labs(title='How far do individual occupancy estimates miss?',subtitle='Mean absolute error at the same original 100 sites. Smaller values are better.',x='Average absolute occupancy error (percentage points)',y=NULL,caption=paste('The size of each error is counted before averaging, so overestimates and underestimates cannot cancel.',
 'Each point averages ten paired communities. Bars show 95% t intervals for those means.',
 'Extra sampling has a real cost: 200, 400 and 600 field samples, respectively. The control receives hidden site conditions that a normal survey would not know exactly.',sep='\n'))+
 base_theme+theme(legend.position='none')
ggsave(file.path(out,'occupancy-absolute-error-by-design.png'),g,width=12,height=6.8,dpi=160,bg='white')
cat('Rendered two figures from the completed design summaries.\n')
z<-subset(p,metric=='occupancy' & group=='all')
z$reduction<- -100*z$change_mae;z$lower<-100*(-z$change_mae-qt(.975,9)*z$change_mae_se);z$upper<-100*(-z$change_mae+qt(.975,9)*z$change_mae_se)
z$design<-factor(as.character(z$arm),levels=c('knownU','sites300','field4'),labels=c('Hidden site conditions known','300 sites, 2 samples','100 sites, 4 samples'))
g<-ggplot(z,aes(reduction,design,colour=arm))+geom_vline(xintercept=0,linetype=2,colour='#778899')+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.16)+geom_point(size=3)+facet_wrap(~contamination,nrow=1)+
 scale_colour_manual(values=colours)+labs(title='How much does each change reduce the error?',subtitle='Each comparison uses the same community and the same original 100 sites.',x='Reduction in average absolute occupancy error\n(percentage points; positive means improvement)',y=NULL,caption=paste('Illustration: a reduction of 3 points would turn an 18-point error into a 15-point error.',
 'Points and 95% t intervals summarize the ten within-community changes, preserving the pairing.',
 'Four field samples double the baseline sampling/PCR effort; 300 sites triple it. The control adds idealized information.',sep='\n'))+
 base_theme+theme(legend.position='none')
ggsave(file.path(out,'occupancy-error-reduction-by-design.png'),g,width=12,height=6,dpi=160,bg='white')
cat('Rendered paired improvement figure.\n')
