#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
library(grid)
# Defaults read the committed summary and write figures beside it. Optional:
# Rscript plot_occupancy_error_explainer.R SUMMARY_CSV OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
script_dir <- dirname(normalizePath(script))
input <- if(length(args)>=1L) args[1] else file.path(script_dir,'results','summary.csv')
root <- if(length(args)>=2L) args[2] else file.path(script_dir,'results')
dir.create(root,recursive=TRUE,showWarnings=FALSE)
s <- read.csv(input)
d <- s[s$scenario %in% c('qnear_K6','qfar_K6') & s$prior=='q20' & s$metric=='occupancy' & s$group %in% c('low','medium','high'),]
stopifnot(nrow(d)==6L, all(d$datasets==10L), max(abs(d$estimate-d$truth-d$bias))<1e-12)
d$case <- factor(d$scenario, levels=c('qnear_K6','qfar_K6'), labels=c('Low contamination','High contamination'))
d$y <- match(d$group,c('high','medium','low'))
d$truth_percent <- 100*d$truth;d$estimate_percent <- 100*d$estimate
blue <- '#146b99'; black <- '#252b31'; teal <- '#21806e'; orange <- '#bd6429'
theme_clear <- theme_minimal(base_size=14) + theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
  plot.title=element_text(face='bold',size=20,margin=margin(b=9)), plot.subtitle=element_text(size=13,margin=margin(b=20)),
  strip.text=element_text(face='bold',size=16,margin=margin(b=16)), axis.title=element_text(size=14), axis.text=element_text(size=13),
  plot.caption=element_text(hjust=0,size=11,lineheight=1.2,margin=margin(t=16)),plot.margin=margin(20,24,15,15))
fig1 <- ggplot(d,aes(y=y)) +
  geom_segment(aes(x=truth_percent,xend=estimate_percent,yend=y),linewidth=1.5,colour=blue,arrow=arrow(length=unit(.1,'inches'),type='closed')) +
  geom_point(aes(x=truth_percent),shape=21,size=4,stroke=1.3,fill='white',colour=black)+
  geom_point(aes(x=estimate_percent),shape=16,size=3.5,colour=blue)+
  geom_text(aes(x=truth_percent,y=y+.20,label=sprintf('True %.1f%%',truth_percent)),size=4.2,colour=black)+
  geom_text(aes(x=estimate_percent,y=y-.20,label=sprintf('Estimated %.1f%%',estimate_percent)),size=4.2,colour=blue)+
  scale_x_continuous(limits=c(-10,110),breaks=c(0,20,40,60,80,100),labels=function(x)paste0(x,'%'),expand=c(0,0))+
  scale_y_continuous(limits=c(.55,3.45),breaks=1:3,labels=c('High\n(true >80%)','Middle\n(true 20-80%)','Low\n(true <20%)'))+
  facet_wrap(~case,nrow=1)+
  labs(title='Low probabilities move up; high probabilities move down',
       subtitle='Real simulation results: six PCR replicates per primer, default priors, no spatial effects',
       x='Occupancy probability (this axis shows probabilities, not errors)',y=NULL,
       caption='Each arrow goes from the average true probability to the average estimated probability for that group.\nAverages across 10 simulated datasets; these are group means, not individual sites or uncertainty intervals.')+theme_clear

ggsave(file.path(root,'true-and-estimated-probabilities.png'),fig1,width=12,height=6.0,dpi=170,bg='white')

# A transparent teaching example, not observations selected from the study.
toy <- data.frame(site=c('Site A','Site B'),y=c(2,1),truth=c(30,70),estimate=c(40,60),error=c(10,-10))
fig2 <- ggplot(toy,aes(y=y))+
  geom_segment(aes(x=truth,xend=estimate,yend=y,colour=site),linewidth=1.8,arrow=arrow(length=unit(.12,'inches'),type='closed'))+
  geom_point(aes(x=truth),shape=21,size=4,stroke=1.2,fill='white',colour=black)+
  geom_point(aes(x=estimate,colour=site),size=4)+
  geom_text(aes(x=truth,y=y+.20,label=paste0('True ',truth,'%')),colour=black,size=5)+
  geom_text(aes(x=estimate,y=y-.21,label=paste0('Estimated ',estimate,'%'),colour=site),size=5)+
  geom_text(aes(x=85,y=y,label=ifelse(error>0,'10 points\ntoo high','10 points\ntoo low'),colour=site),size=5,fontface='bold')+
  annotate('label',x=27,y=.06,label='Average direction of error\n(+10 - 10) / 2 = 0 points',size=4.6,colour=black,fill='#f3f5f6',linewidth=0,lineheight=1.2)+
  annotate('label',x=70,y=.06,label='Average size of error\n(10 + 10) / 2 = 10 points',size=4.6,colour=black,fill='#f3f5f6',linewidth=0,lineheight=1.2)+
  scale_colour_manual(values=c('Site A'=orange,'Site B'=teal),guide='none')+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20),labels=function(x)paste0(x,'%'),expand=c(0,0))+
  scale_y_continuous(limits=c(-.3,2.5),breaks=c(1,2),labels=c('Site B','Site A'))+
  labs(title='Two wrong estimates can have an average error of zero',
       subtitle='Made-up example for one species at two sites; both true probabilities are between 20% and 80%',
       x='Occupancy probability',y=NULL,
       caption='The opposite errors cancel when we keep their signs. They cannot cancel when we count how far each estimate missed.\nThis example illustrates the calculation; these two sites are not results from the simulation study.')+theme_clear+
  theme(panel.grid.major.x=element_blank())
ggsave(file.path(root,'how-errors-cancel.png'),fig2,width=11,height=5.7,dpi=170,bg='white')

# The actual middle-group numbers, in two side-by-side comparisons.
m <- d[d$group=='medium',]
bars <- rbind(data.frame(case=m$case,value=100*m$bias,measure='Average direction of error\n(positive and negative errors can cancel)'),
              data.frame(case=m$case,value=100*m$mae,measure='Average size of error\n(each miss counts as positive)'))
bars$measure <- factor(bars$measure,levels=unique(bars$measure))
bars$case <- factor(bars$case,levels=c('High contamination','Low contamination'))
bars$label <- sprintf('%.1f points',bars$value)
fig3 <- ggplot(bars,aes(x=value,y=case,fill=measure))+
  geom_vline(xintercept=0,colour=black,linewidth=.5)+
  geom_col(width=.4)+
  geom_text(aes(x=ifelse(value<0,value-.8,value+.7),label=label,hjust=ifelse(value<0,1,0)),size=4.2,lineheight=1.1)+
  facet_wrap(~measure,nrow=1)+
  scale_fill_manual(values=c(blue,orange),guide='none')+
  scale_x_continuous(limits=c(-12,25),breaks=c(-10,0,10,20))+
  labs(title='The actual middle group: small average bias, much larger misses',
       subtitle='True occupancy probabilities 20-80%; six PCR replicates per primer; default priors',
       x='Error in percentage points (both panels use the same scale)',y=NULL,
       caption='Average size of error is the mean absolute error (MAE). It does not mean every estimate misses by that amount.\nCalculated from the same saved fits as the first figure, averaging over 10 simulated datasets per scenario.')+theme_clear+
  theme(strip.text=element_text(face='bold',size=13),panel.spacing=unit(1.2,'lines'))
ggsave(file.path(root,'middle-group-bias-and-error-size.png'),fig3,width=12,height=4.8,dpi=170,bg='white')
write.csv(d[,c('scenario','prior','group','datasets','truth','estimate','bias','between_dataset_se','mae','rmse')],file.path(root,'occupancy-error-explainer-values.csv'),row.names=FALSE)
