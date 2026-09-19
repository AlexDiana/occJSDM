here<-'work/jsdm-sample-size-20260919';source(file.path(here,'helpers.R'));frozen<-load_frozen('work/nonspatial-bias-recheck-20260914')
m<-error_metrics(c(.1,.9),c(.3,.7));stopifnot(abs(m['bias'])<1e-12,abs(m['mae']-.2)<1e-12,abs(m['rmse']-.2)<1e-12)
for(n in c(100L,300L,1000L)) {
 key<-sprintf('n%04d-01',n);fit<-readRDS(file.path(here,'pilot',paste0(key,'-fit.rds')))$fit
 input<-readRDS(file.path(here,'inputs',paste0(key,'.rds')));r<-score_fit(fit,input)
 # Binary fits do not expose a public psi mean. Missing is explicit, not -Inf.
 stopifnot(is.null(fit$results_output$psi_output),is.na(r$public_mean_difference))
 # Independently reconstruct every pilot draw, then check the original-site errors.
 j<-fit$results_output$jsdm_output;pd<-array(NA_real_,c(n,10L,80L,2L))
 for(ch in 1:2)for(it in 1:80)for(s in 1:10)pd[,s,it,ch]<-plogis(j$B0_output[s,it,ch]+drop(fit$X_psi%*%j$B_output[,s,it,ch])+drop(j$U_output[,,it,ch]%*%j$L_output[,s,it,ch]))
 mu<-apply(pd,c(1,2),mean);truth<-plogis(input$truth$eta)
 stopifnot(max(abs(mu-r$estimate))<1e-12)
 expected<-mean(abs(mu[1:100,]-truth[1:100,]));actual<-r$groups$mae[r$groups$scope=='original100' & r$groups$group=='all']
 stopifnot(abs(expected-actual)<1e-12)
 for(k in seq_len(nrow(r$groups))) {
  g<-r$groups[k,];ix<-matrix(FALSE,n,10);ix[if(g$scope=='original100')1:100 else 1:n,]<-TRUE
  if(g$group=='low')ix<-ix & truth<.2
  if(g$group=='middle')ix<-ix & truth>=.2 & truth<=.8
  if(g$group=='high')ix<-ix & truth>.8
  tr<-apply(pd,c(3,4),function(p)mean(p[ix]))
  stopifnot(max(abs(tr-r$traces[paste(g$scope,g$group,sep=':'),,]))<1e-12)
 }
}
cat('PASS: signed versus absolute errors, all pilot posterior means and every group trace independently reconstructed.\n')
