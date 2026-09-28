args<-commandArgs(TRUE);stopifnot(length(args)==2L)
repo<-normalizePath(args[1]);raw<-normalizePath(args[2]);folder<-file.path(raw,'summary')
marker<-readRDS(file.path(folder,'audit-complete.rds'))
stopifnot(identical(unname(tools::md5sum(marker$artifact_hashes$file)),marker$artifact_hashes$md5),
  identical(unname(tools::md5sum(marker$analysis_hashes$file)),marker$analysis_hashes$md5),
  identical(unname(tools::md5sum(marker$selected$file)),marker$selected$result_md5))
library(ggplot2)
aggregate<-read.csv(file.path(folder,'aggregate.csv'));sites<-read.csv(file.path(folder,'sites.csv'))
selected<-read.csv(file.path(raw,'all-selected.csv'))
flag_caption<-paste('Unresolved convergence flags:',paste(vapply(LETTERS[1:4],function(cfg) {
  rows<-selected[selected$configuration==cfg,]
  sprintf('%s %d/%d',cfg,sum(rows$flag_count>0),nrow(rows))
},character(1)),collapse='; '),'. Flagged results are provisional.')
labels<-c(A='A: amplitude 0.32',B='B: true amplitude 1',C='C: estimate intercept',D='D: ten ecological states')
group_labels<-c(all='All species','0.01'='1% occupancy','0.05'='5% occupancy',
  '0.25'='25% occupancy','0.75'='75% occupancy')
colours<-c(A='#666666',B='#225EA8',C='#D95F0E',D='#238B45')
base_theme<-theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=10),
  legend.position='none',strip.text=element_text(face='bold'))
draw_comparison<-function(metric,title,axis,filename,scale=1,zero=FALSE) {
  a<-aggregate[aggregate$metric==metric,]
  a$configuration<-factor(a$configuration,levels=rev(LETTERS[1:4]))
  a$group<-factor(a$group,levels=names(group_labels),labels=group_labels)
  a[c('mean','lower','upper')]<-a[c('mean','lower','upper')]*scale
  g<-ggplot(a,aes(y=configuration,x=mean,colour=configuration))+
    geom_segment(aes(x=lower,xend=upper,yend=configuration),linewidth=.75)+
    geom_point(size=2.8)+facet_grid(~group)+scale_y_discrete(labels=labels)+
    scale_colour_manual(values=colours)+labs(title=title,x=axis,y=NULL,
      subtitle='Nine saved communities; range and environmental coefficient supplied at truth',
      caption=paste0('B is the reference: A changes amplitude; C estimates the intercept under N(0,1); D adds nine independent occupancy states per site.\n',
        'Intervals describe variation across communities, not posterior uncertainty. ',
        if(zero)'Dashed lines show error from estimating a zero spatial field.' else 'Probability errors are in percentage points.',
        '\n',flag_caption))+base_theme
  if(zero) {
    z<-aggregate[aggregate$metric=='zero_centred_rmse' & aggregate$configuration=='B',]
    z$group<-factor(z$group,levels=names(group_labels),labels=group_labels)
    g<-g+geom_vline(data=z,aes(xintercept=mean),linetype=2,colour='#444444')
  }
  ggsave(file.path(folder,filename),g,width=15,height=5.1,dpi=160,bg='white')
}
draw_comparison('centred_rmse','What limits recovery of the spatial pattern?',
  'Centred field RMSE (lower is better)','conditional-field-error.png',zero=TRUE)
draw_comparison('occupancy_mae','How do the conditional changes affect occupancy estimates?',
  'Occupancy mean absolute error (percentage points)','conditional-occupancy-error.png',scale=100)

# Same prespecified illustrative community/species as the main prior report.
m<-sites[sites$community=='range6-rep01' & sites$species %in% c(2,4,6,8),]
m$value<-m$field_median-ave(m$field_median,interaction(m$species,m$configuration),FUN=mean)
m$method<-m$configuration
truth<-m[m$configuration=='B',];truth$value<-truth$field_truth-ave(truth$field_truth,truth$species,FUN=mean)
truth$method<-'Truth';m<-rbind(truth,m)
m$method<-factor(m$method,levels=c('Truth',LETTERS[1:4]),
  labels=c('Truth','A: amplitude 0.32','B: true amplitude 1','C: estimate intercept','D: ten ecological states'))
m$species_label<-factor(m$species,levels=c(2,4,6,8),labels=c('1% occupancy','5% occupancy','25% occupancy','75% occupancy'))
limit<-max(abs(m$value))
g<-ggplot(m,aes(x,y,fill=value))+geom_point(shape=21,size=2.8,colour='transparent',stroke=0)+
  facet_grid(species_label~method)+coord_equal()+
  scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-limit,limit),name='Centred field')+
  labs(title='Known parameters and extra ecological observations: an illustrative community',
    subtitle='Middle spatial range, community 1; one prespecified species from each prevalence group',
    x='Standardized x coordinate',y='Standardized y coordinate',
    caption=paste('Each map is centred across sites. The shared colour scale shows both pattern and magnitude. D adds independent ecological states, not PCR replicates.',
      flag_caption,sep='\n'))+
  theme_minimal(base_size=11)+theme(panel.grid=element_blank(),plot.title=element_text(face='bold'),
    strip.text=element_text(face='bold'),plot.caption=element_text(hjust=0),legend.position='right')
ggsave(file.path(folder,'conditional-field-maps.png'),g,width=15,height=11,dpi=160,bg='white')

# Export compact audited evidence; full draws and diagnostics remain in raw.
output<-file.path(repo,'dev/simstudy/spatial-amplitude-prior/diagnosis/results')
dir.create(output,showWarnings=FALSE)
copy_csv<-function(path) {
  x<-read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  for(nm in names(x))if(is.character(x[[nm]])) {
    x[[nm]]<-ifelse(startsWith(x[[nm]],paste0(raw,'/')),
      substring(x[[nm]],nchar(raw)+2L),x[[nm]])
    x[[nm]]<-ifelse(startsWith(x[[nm]],paste0(repo,'/')),
      substring(x[[nm]],nchar(repo)+2L),x[[nm]])
  }
  write.csv(x,file.path(output,basename(path)),row.names=FALSE)
}
for(name in c('aggregate.csv','paired.csv','community-groups.csv','chain-community-groups.csv',
  'chain-community-extrema.csv','chain-sensitivity.csv','audit.csv','analysis-source-md5.csv'))
  copy_csv(file.path(folder,name))
for(name in c('all-selected.csv','geometry.csv','information.csv','source-md5.csv','quadrature-checks.csv','validation-source-md5.csv'))
  copy_csv(file.path(raw,name))
for(name in c('conditional-field-error.png','conditional-occupancy-error.png','conditional-field-maps.png'))
  stopifnot(file.copy(file.path(folder,name),file.path(output,name),overwrite=TRUE))
write.csv(m,file.path(output,'illustrative-sites.csv'),row.names=FALSE)
scripts<-list.files(dirname(output),pattern='[.](R|cpp|md)$',full.names=TRUE)
write.csv(data.frame(file=substring(scripts,nchar(repo)+2L),md5=unname(tools::md5sum(scripts))),
  file.path(output,'research-source-md5.csv'),row.names=FALSE)
writeLines(c('# Compact conditional spatial diagnosis','',
  paste('The full draw archive is',raw),
  'The A/B/C/D targets are specified in ../PLAN.md. Read ../REPORT.md for conclusions and limitations.',
  'Every row of all-selected.csv identifies one of 288 independently sampled conditional posteriors; audit.csv verifies all of them.',
  'Result paths are relative to the raw archive. Source paths are repository-relative. Inputs retain their original nine-community identities.',
  'The paired table reports A, C or D minus B within each community, with the three generating ranges as fixed strata.',
  'Probability errors are proportions in the tables; the figure multiplies them by 100. Containment from nine communities does not establish calibration.',
  'Full pointwise diagnostics, all site summaries, initial/longer draws, settings and seeds are preserved in the raw archive.'),
  file.path(output,'README.md'))
files<-list.files(output,full.names=TRUE);files<-files[basename(files)!='artifact-md5.csv']
write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),
  file.path(output,'artifact-md5.csv'),row.names=FALSE)
cat('Plotted and exported compact conditional diagnosis.\n')
