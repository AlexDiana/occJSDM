#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
option <- function(name,default=NULL) {
  hit <- grep(paste0('^--',name,'='),args,value=TRUE)
  if(length(hit)>1L)stop('Repeated option')
  if(!length(hit)){if(is.null(default))stop('Missing --',name);return(default)}
  substring(hit,nchar(name)+4L)
}
input <- normalizePath(option('summary'));out <- option('out',input)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(ggplot2))
a <- read.csv(file.path(input,'aggregate.csv'),stringsAsFactors=FALSE)
groups <- read.csv(file.path(input,'groups.csv'),stringsAsFactors=FALSE)
labels <- c(binary='True occupancy states supplied',low='Two-stage: low contamination',high='Two-stage: high contamination')
colors <- c('True occupancy states supplied'='#343d4d','Two-stage: low contamination'='#147d92','Two-stage: high contamination'='#c56328')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
  legend.position='bottom',plot.title.position='plot',strip.text=element_text(face='bold'),
  plot.caption=element_text(hjust=0,size=10)))
decorate <- function(d) {
  d$arm <- factor(labels[d$arm],levels=unname(labels))
  d$support <- factor(d$knots,levels=c(20,50,100),labels=c('20 (default)','50','100 (all sites)'))
  d
}
b <- decorate(a[a$metric=='occupancy' & a$quantity=='bias' & a$group %in% c('low','medium','high'),])
b$band <- factor(b$group,levels=c('low','medium','high'),
  labels=c('Truth below 20%','Truth 20% to 80%','Truth above 80%'))
g <- decorate(groups[groups$metric=='occupancy' & groups$group %in% c('low','medium','high'),])
g$band <- factor(g$group,levels=c('low','medium','high'),labels=levels(b$band))
pd <- position_dodge(width=.38)
p <- ggplot(b,aes(support,100*mean,color=arm,group=arm))+
  annotate('rect',xmin=-Inf,xmax=Inf,ymin=-5,ymax=5,fill='#e1eddf',alpha=.6)+
  geom_hline(yintercept=0,color='grey45',linewidth=.4)+
  geom_point(data=g,aes(support,100*bias,color=arm,group=arm),position=pd,alpha=.25,size=1.3)+
  geom_line(position=pd,linewidth=.6)+
  geom_errorbar(aes(ymin=100*lower,ymax=100*upper),position=pd,width=.12)+
  geom_point(position=pd,size=2.5)+facet_wrap(~band,nrow=1)+scale_color_manual(values=colors)+
  labs(title='Spatial occupancy bias under feasible sampling',
    subtitle='Nine independent communities; 100 sites, two field samples, two primers, six PCR replicates per primer',
    x='Spatial support points',y='Occupancy bias (percentage points)',color=NULL,
    caption='Faint points: individual communities. Bars: approximate 95% intervals using variation within each spatial range.\nGreen band: provisional five-point target for average bias; it is not an individual-prediction or rare-species guarantee.')
ggsave(file.path(out,'occupancy-bias.png'),p,width=12,height=5.3,dpi=180,bg='white')
r <- decorate(a[a$metric=='occupancy' & a$quantity=='bias' & a$group %in% c('prevalence_1pct','prevalence_5pct'),])
r$target <- ifelse(r$group=='prevalence_1pct',.01,.05)
r$prevalence <- factor(r$target,levels=c(.01,.05),labels=c('True mean occupancy: 1%','True mean occupancy: 5%'))
lines <- unique(r[c('prevalence','target')])
p <- ggplot(r,aes(support,100*(mean+target),color=arm,group=arm))+
  geom_hline(data=lines,aes(yintercept=100*target),inherit.aes=FALSE,linetype=2,color='grey35')+
  geom_line(position=pd)+geom_errorbar(aes(ymin=100*(lower+target),ymax=100*(upper+target)),position=pd,width=.12)+
  geom_point(position=pd,size=2.5)+facet_wrap(~prevalence)+scale_color_manual(values=colors)+
  labs(title='Are deliberately rare species estimated as rare?',x='Spatial support points',
    y='Estimated mean occupancy (%)',color=NULL,
    caption='Each community contains two species at each rare prevalence. Dashed lines show truth.\nZero-occupancy and zero-detection species are retained. Bars describe simulation uncertainty, not posterior intervals.')
ggsave(file.path(out,'rare-species.png'),p,width=10.5,height=5.3,dpi=180,bg='white')
c <- decorate(a[a$metric=='occupancy' & a$quantity=='coverage' & a$group %in%
  c('low','medium','high','prevalence_1pct','prevalence_5pct'),])
c$band <- factor(c$group,levels=c('low','medium','high','prevalence_1pct','prevalence_5pct'),
  labels=c('Truth below 20%','Truth 20% to 80%','Truth above 80%','Species prevalence 1%','Species prevalence 5%'))
p <- ggplot(c,aes(support,100*mean,color=arm,group=arm))+
  geom_hline(yintercept=95,linetype=2,color='grey40')+
  geom_line(position=pd)+geom_point(position=pd,size=2.3)+
  facet_wrap(~band,nrow=2)+scale_color_manual(values=colors)+coord_cartesian(ylim=c(0,100))+
  labs(title='How often did the saved 95% intervals contain the true probability?',
    subtitle='Descriptive containment across nine independent communities',x='Spatial support points',
    y='Intervals containing truth (%)',color=NULL,
    caption='Dashed line: nominal 95%. Site/species intervals within a community are dependent.\nThis small study does not establish precise repeated-sampling coverage; flagged fits require caution.')
ggsave(file.path(out,'interval-containment.png'),p,width=12,height=7,dpi=180,bg='white')

range_points <- decorate(groups[groups$metric=='range',])
range_cells <- split(range_points,interaction(range_points$arm,range_points$support,range_points$grid_index,drop=TRUE))
range_means <- do.call(rbind,lapply(range_cells,function(d) {
  stopifnot(nrow(d)==3L)
  data.frame(arm=d$arm[1],support=d$support[1],truth=d$truth[1],
    estimate=mean(d$estimate),se=sd(d$estimate)/sqrt(nrow(d)))
}))
p <- ggplot(range_means,aes(truth,estimate,color=support,group=support))+
  geom_abline(slope=1,intercept=0,linetype=2,color='grey40')+
  geom_point(data=range_points,alpha=.35,size=1.5)+
  geom_line(linewidth=.6)+
  geom_errorbar(aes(ymin=estimate-se,ymax=estimate+se),width=.005)+
  geom_point(size=2.5)+facet_wrap(~arm,nrow=1)+
  scale_color_manual(values=c('20 (default)'='#777777','50'='#147d92','100 (all sites)'='#783f89'))+
  scale_x_continuous(breaks=sort(unique(range_points$truth)),labels=function(x)sprintf('%.3f',x))+
  coord_cartesian(ylim=range(c(range_points$truth,range_points$estimate)))+
  labs(title='Did estimated spatial ranges track the three different truths?',
    x='True spatial range (standardized coordinates)',y='Estimated spatial range',color='Support points',
    caption='Dashed line: exact recovery. Faint points: independent communities; larger points: three-community means.\nBars show one between-community standard error at each fixed range, not posterior credible intervals.')
ggsave(file.path(out,'spatial-range.png'),p,width=12,height=5.3,dpi=180,bg='white')
