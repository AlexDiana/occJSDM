#!/usr/bin/env Rscript
a<-commandArgs(TRUE);stopifnot(length(a)==2L)
repo<-normalizePath(a[1]);summary<-normalizePath(a[2])
library(ggplot2)
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
fits<-read.csv(file.path(summary,'fits.csv'),stringsAsFactors=FALSE)
groups<-read.csv(file.path(summary,'groups.csv'),stringsAsFactors=FALSE)
labels<-c(inverse_gamma='Current inverse-gamma',half_cauchy='Half-Cauchy(1)')
colours<-c('Current inverse-gamma'='#3874A5','Half-Cauchy(1)'='#D76B26')
theme_set(theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
metrics<-c(raw_rmse='Spatial field: total error',centred_rmse='Spatial pattern: centred error',
  centred_correlation='Spatial pattern: correlation',amplitude_mean='Spatial amplitude estimate')
long<-do.call(rbind,lapply(names(metrics),function(m)data.frame(fits[c('community','grid_index','arm','prior')],
  metric=metrics[[m]],value=fits[[m]])))
long$prior<-factor(labels[long$prior],levels=labels)
p<-ggplot(long[long$arm=='binary',],aes(prior,value,group=community,colour=factor(grid_index)))+
  geom_line(alpha=.55)+geom_point(size=2.4)+facet_wrap(~metric,scales='free_y',ncol=2)+
  labs(x=NULL,y=NULL,colour='True range grid position',
    title='Does a broader spatial prior recover the field better?',
    subtitle='Each line follows the same simulated community under both priors.',
    caption='Lower error and higher correlation are better. Amplitude alone is not a recovery criterion.')+
  scale_colour_manual(values=c('4'='#386CB0','6'='#7F5D9C','8'='#308778'))
ggsave(file.path(summary,'paired-field-recovery.png'),p,width=10,height=7,dpi=180)
# The prior is on the coefficient SD; a realized field has its own spatial variation.
x<-seq(.001,3.5,length.out=1500)
prior<-rbind(data.frame(x=x,density=2/x^3*dgamma(1/x^2,shape=10,rate=1),prior=labels[1]),
  data.frame(x=x,density=2*dcauchy(x,scale=1),prior=labels[2]))
p<-ggplot(prior,aes(x,density,colour=prior))+geom_line(linewidth=.9)+
  scale_colour_manual(values=colours)+labs(x='Spatial coefficient standard deviation',y='Prior density',colour=NULL,
    title='The experiment weakens a concentrated spatial-amplitude prior',
    subtitle='Current prior: variance ~ inverse-gamma(10, 1). Experimental prior: SD ~ half-Cauchy(1).',
    caption='Display ends at 3.5; the half-Cauchy retains a long right tail beyond the plot.')
ggsave(file.path(summary,'spatial-amplitude-priors.png'),p,width=10,height=4.5,dpi=180)
# Truth bands remain separate from whole-species prevalence groups.
d<-groups[groups$arm=='binary' & groups$metric=='occupancy' & groups$group %in% c('low','medium','high'),]
stats<-do.call(rbind,lapply(split(d,interaction(d$prior,d$group,drop=TRUE)),function(z)
  data.frame(prior=z$prior[1],group=z$group[1],as.list(stratified_mean(z$bias,z$grid_index)))))
d$group<-factor(d$group,c('low','medium','high'),c('Below 20%','20% to 80%','Above 80%'))
stats$group<-factor(stats$group,c('low','medium','high'),c('Below 20%','20% to 80%','Above 80%'))
d$prior<-factor(labels[d$prior],levels=labels);stats$prior<-factor(labels[stats$prior],levels=labels)
p<-ggplot(d,aes(prior,100*bias,colour=prior))+geom_hline(yintercept=0,colour='grey55',linetype=2)+
  geom_point(alpha=.5,position=position_jitter(width=.06,height=0,seed=11))+
  geom_errorbar(data=stats,aes(y=100*mean,ymin=100*lower,ymax=100*upper),width=.08,linewidth=.8)+
  geom_point(data=stats,aes(y=100*mean),shape=18,size=4)+facet_wrap(~group)+
  scale_colour_manual(values=colours)+labs(x=NULL,y='Occupancy bias (percentage points)',colour=NULL,
    title='Occupancy bias at the fitted sites',
    subtitle='Dots are communities; diamonds and bars show the stratified mean and descriptive 95% interval.')+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave(file.path(summary,'occupancy-bias.png'),p,width=10,height=4.5,dpi=180)
# One prespecified middle-range community, one species per prevalence group.
selected<-readRDS(file.path(summary,'selected-results.rds'))
key<-'range6-rep01-binary-k100';ig<-selected[[paste0(key,':inverse_gamma')]];hc<-selected[[paste0(key,':half_cauchy')]]
if(!is.null(ig) && !is.null(hc)) {
  input<-readRDS(ig$job$input_file);field<-list(Truth=input$truth$field,
    'Current inverse-gamma'=ig$spatial$field_mean,'Half-Cauchy(1)'=hc$spatial$field_mean)
  xy<-input$truth$Xs;species<-c(2,4,6,8)
  maps<-do.call(rbind,lapply(names(field),function(prior)do.call(rbind,lapply(species,function(s){
    z<-field[[prior]][,s];data.frame(x=xy[,1],y=xy[,2],value=z-mean(z),prior=prior,
      species=paste0('Species ',s,' (',input$truth$target_prevalence[s]*100,'% prevalence)'))
  }))))
  maps$prior<-factor(maps$prior,names(field))
  p<-ggplot(maps,aes(x,y,colour=value))+geom_point(size=2.5)+coord_equal()+
    facet_grid(species~prior)+scale_colour_gradient2(low='#2B6CA3',mid='white',high='#C34732',midpoint=0)+
    labs(x='Standardized longitude',y='Standardized latitude',colour='Centred field',
      title='Spatial patterns for one prespecified community',
      subtitle='Range grid position 6, replicate 1; each species is centred over sites.')
  ggsave(file.path(summary,'illustrative-spatial-patterns.png'),p,width=9,height=10,dpi=180)
}
