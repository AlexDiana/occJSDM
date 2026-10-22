args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L,startsWith(args,'--results='))
out <- normalizePath(sub('^--results=','',args))
suppressPackageStartupMessages(library(ggplot2))
d <- read.csv(file.path(out,'oracle-aggregate.csv'))
x <- rbind(d[d$target==.01&d$arm=='binary',],d[d$target==.01&d$arm=='field'&d$nuisance=='known',],
 d[d$target==.01&d$nuisance=='collection_intercept_unknown'&d$arm=='field',],
 d[d$target==.01&d$nuisance=='collection_intercept_unknown'&d$arm=='low',],
 d[d$target==.01&d$nuisance=='collection_intercept_unknown'&d$arm=='high',])
labels <- c('True site states','True field-sample states\nCollection probability known',
 'True field-sample states\nCollection intercept estimated','Low-contamination PCR\nCollection intercept estimated',
 'High-contamination PCR\nCollection intercept estimated')
x$control <- factor(rep(labels,each=2),levels=rev(labels))
x$prior <- factor(x$prior_sd,levels=c(1,2.5),labels=c('Default: SD 1','Sensitivity probe: SD 2.5'))
g <- ggplot(x,aes(x=estimate*100,y=control,fill=prior))+
 geom_col(position=position_dodge(width=.75),width=.66)+
 geom_text(aes(label=sprintf('%.1f%%',100*estimate)),position=position_dodge(width=.75),hjust=-.12,size=3.6)+
 geom_vline(xintercept=1,linetype=2,color='#222222',linewidth=.5)+
 scale_fill_manual(values=c('#B56340','#3C718A'))+scale_x_continuous(limits=c(0,36),breaks=seq(0,35,5),expand=expansion(mult=c(0,.01)))+
 labs(title='Why a 1%-occupancy species can look common',
 subtitle='Controlled calculations on the nine saved communities; dashed line = 1% truth',
 x='Mean estimated occupancy (%)',y=NULL,fill='Occupancy-intercept prior',
 caption='Spatial effects and other nuisance parameters are supplied at truth.\nThese are conditional diagnostics, not full-model refits. The wider prior does not eliminate bias.')+
 theme_minimal(base_size=12)+theme(legend.position='bottom',plot.title=element_text(face='bold',size=16),
 panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0),
 plot.margin=margin(12,18,12,12))
ggsave(file.path(out,'prior-information-diagnosis.png'),g,width=10.5,height=6.4,dpi=180,bg='white')
