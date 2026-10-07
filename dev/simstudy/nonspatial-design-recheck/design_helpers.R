load_frozen <- function(root) {
  .libPaths(c(file.path(root,'library'),.libPaths()))
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  suppressPackageStartupMessages(library(occJSDM))
  RcppParallel::setThreadOptions(numThreads=1)
  source(file.path(root,'source-main/tests/testthat/helper-simstudy.R'),local=.GlobalEnv)
  original<-readRDS(file.path(root,'operational-k6/settings.rds'))
  stopifnot(identical(unname(tools::md5sum(names(original$hashes))),unname(original$hashes)))
  expr<-as.list(parse(file.path(root,'run_nonspatial_recheck_balanced.R')))
  for(nm in c('score_fit','trace_summary')) {
    hit<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name(nm)),expr)
    stopifnot(length(hit)==1);eval(hit[[1]],envir=.GlobalEnv)
  }
  invisible(original)
}

make_design_input <- function(a,arm) {
  stopifnot(arm %in% c('field4','sites300','knownU'),a$scenario$n==100,a$scenario$M==2,a$scenario$K==6,!a$scenario$useSpatField)
  b<-a;b$design<-list(arm=arm,original_sites=100L,base_seed=a$seed)
  if(arm=='knownU') {
    b$design$original_rows<-seq_len(nrow(a$sim$data_list$info));b$design$original_samples<-seq_len(200L)
    return(b)
  }
  seed_offset<-if(arm=='field4')1000000L else 2000000L
  b$design$extension_seeds<-a$seed+seed_offset+c(sites=10000L,field=20000L,reads=30000L)
  n<-if(arm=='sites300')300L else 100L;M<-if(arm=='field4')4L else 2L;S<-a$scenario$S;P<-2L;K<-6L
  oldinfo<-a$sim$data_list$info;siteinfo<-oldinfo[!duplicated(oldinfo$Site),];siteinfo<-siteinfo[order(siteinfo$Site),]
  cols<-paste0('X_psi.EnvCov.',1:2);xold<-as.matrix(siteinfo[,cols]);x<-xold
  xs<-as.matrix(siteinfo[,c('Xs.1','Xs.2')]);tp<-a$sim$true_params;jp<-tp$jsdmParams_true
  if(n>100L) {
    set.seed(b$design$extension_seeds['sites'])
    x<-rbind(xold,matrix(rnorm((n-100L)*2L,sd=10),n-100L,2L))
    xs<-rbind(xs,matrix(runif((n-100L)*2L),n-100L,2L))
    jp$U<-rbind(jp$U,matrix(rnorm((n-100L)*2L,sd=a$truth$jsdmParams$sigma_h),n-100L,2L))
    oldsd<-apply(xold,2,sd);newsd<-apply(x,2,sd);ratio<-newsd/oldsd
    jp$B0<-jp$B0+as.vector(((colMeans(x)-colMeans(xold))/oldsd)%*%jp$B)
    jp$B<-sweep(jp$B,1,ratio,'*');jp$G<-sweep(jp$G,2,ratio,'*');jp$C<-sweep(jp$C,2,ratio,'*');jp$Bt<-sweep(jp$Bt,2,ratio,'*')
    jp$eta<-sweep(scale(x)%*%jp$B+jp$U%*%jp$L,2,jp$B0,'+')
    tp$z_true<-rbind(tp$z_true,matrix(rbinom((n-100L)*S,1,plogis(jp$eta[101:n,])),n-100L,S))
    jp$spatField<-matrix(0,n,S);jp$SE<-matrix(0,n,S);jp$varPart<-NULL
  }
  b$scenario$n<-n;b$scenario$M<-M
  b$truth$datasettings$n<-n;b$truth$datasettings$M<-rep(M,n);b$truth$datasettings$K<-rep(K,n*M*P)
  sample_site<-rep(seq_len(n),each=M);slot<-rep(seq_len(M),times=n)
  oldsample<-ifelse(sample_site<=100 & slot<=2,(sample_site-1)*2+slot,0L)
  original_samples<-match(seq_len(200L),oldsample);new_samples<-which(oldsample==0L)
  sampleinfo<-oldinfo[!duplicated(oldinfo$Sample),];sampleinfo<-sampleinfo[order(sampleinfo$Sample),]
  xt<-numeric(n*M);xt[original_samples]<-sampleinfo$X_theta
  w<-matrix(0,n*M,S);w[original_samples,]<-tp$w_true
  set.seed(b$design$extension_seeds['field'])
  xt[new_samples]<-rnorm(length(new_samples))
  theta<-plogis(cbind(1,xt)%*%tp$beta_theta_true)
  znew<-tp$z_true[sample_site[new_samples],,drop=FALSE]
  pr<-znew*theta[new_samples,,drop=FALSE]+(1-znew)*matrix(a$truth$params$theta0,length(new_samples),S,byrow=TRUE)
  w[new_samples,]<-matrix(rbinom(length(pr),1,pr),length(new_samples),S)
  sample_row<-rep(seq_len(n*M),each=P*K);primer<-rep(rep(seq_len(P),each=K),times=n*M)
  site_row<-sample_site[sample_row];pcr<-rep(seq_len(K),times=n*M*P)
  info<-data.frame(Site=site_row,Sample=sample_row,Primer=primer)
  info[,cols]<-x[site_row,,drop=FALSE];info[,c('Xs.1','Xs.2')]<-xs[site_row,,drop=FALSE];info$X_theta<-xt[sample_row]
  oldpcr<-ave(seq_len(nrow(oldinfo)),interaction(oldinfo$Sample,oldinfo$Primer),FUN=seq_along)
  map_sample<-(oldinfo$Site-1)*M+((oldinfo$Sample-1)%%2)+1
  original_rows<-match(paste(map_sample,oldinfo$Primer,oldpcr),paste(sample_row,primer,pcr))
  stopifnot(!anyNA(original_rows),!anyDuplicated(original_rows))
  y<-matrix(0,nrow(info),S,dimnames=list(NULL,colnames(a$sim$data_list$OTU)))
  y[original_rows,]<-a$sim$data_list$OTU
  new_rows<-which(oldsample[sample_row]==0L);nr<-length(new_rows)
  set.seed(b$design$extension_seeds['reads'])
  wr<-w[sample_row[new_rows],,drop=FALSE]
  pr<-wr*tp$p_true[primer[new_rows],,drop=FALSE]+(1-wr)*tp$q_true[primer[new_rows],,drop=FALSE]
  event<-matrix(runif(nr*S),nr,S)<pr
  logpos<-matrix(rnorm(nr*S,5,1),nr,S);logfalse<-matrix(rnorm(nr*S,1.5,1),nr,S)
  logreads<-event*(wr*logpos+(1-wr)*logfalse)
  y[new_rows,]<-pmax(0,round(exp(logreads)-1))
  tp$jsdmParams_true<-jp;tp$w_true<-w;b$sim$true_params<-tp;b$sim$data_list$info<-info;b$sim$data_list$OTU<-y
  b$design$original_rows<-original_rows;b$design$original_samples<-original_samples
  b
}

check_nested_input <- function(a,b,arm) {
  i<-b$sim$data_list$info;old<-a$sim$data_list$info;rows<-b$design$original_rows;samples<-b$design$original_samples
  n<-if(arm=='sites300')300L else 100L;M<-if(arm=='field4')4L else 2L
  stopifnot(nrow(i)==n*M*2L*6L,nrow(b$sim$data_list$OTU)==nrow(i),all(table(i$Sample,i$Primer)==6),
            all(table(unique(i[c('Site','Sample')])$Site)==M),
            identical(unname(b$sim$data_list$OTU[rows,,drop=FALSE]),unname(a$sim$data_list$OTU)),
            max(abs(b$sim$true_params$z_true[1:100,]-a$sim$true_params$z_true))==0,
            max(abs(b$sim$true_params$w_true[samples,]-a$sim$true_params$w_true))==0,
            max(abs(b$sim$true_params$jsdmParams_true$eta[1:100,]-a$sim$true_params$jsdmParams_true$eta))<1e-12,
            identical(b$sim$data_list$traits,a$sim$data_list$traits),identical(b$truth$params,a$truth$params),identical(b$fit_rng,a$fit_rng))
  for(col in setdiff(names(old),'Sample'))stopifnot(all(i[rows,col]==old[,col]))
  for(col in c('p_true','q_true','beta_theta_true'))stopifnot(identical(b$sim$true_params[[col]],a$sim$true_params[[col]]))
  invisible(TRUE)
}

known_u_fitter <- function(U_true) {
  ns<-asNamespace('occJSDM');env<-new.env(parent=ns);env$known_U<-U_true
  rewrite<-function(expr,mode) {
    count<-0L
    visit<-function(e) {
      if(!is.call(e))return(e)
      if(identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('U'))) {
        rhs<-e[[3]]
        yes<-if(mode=='run') identical(rhs,quote(matrix(0,n,d))) else
          identical(rhs,quote(list_params$U)) || (is.call(rhs)&&identical(rhs[[1]],as.name('sample_U_cpp')))
        if(yes) {count<<-count+1L;e[[3]]<-as.name('known_U');return(e)}
      }
      for(j in seq_along(e))if(is.call(e[[j]]))e[[j]]<-visit(e[[j]])
      e
    }
    result<-visit(expr);stopifnot(count==if(mode=='run')1L else 2L);result
  }
  f<-get('runOccJSDM',ns);body(f)<-rewrite(body(f),'run');environment(f)<-env
  u<-get('update_jSDMcoef',ns);body(u)<-rewrite(body(u),'update');environment(u)<-env
  env$update_jSDMcoef<-u
  f
}

make_design_scorer <- function() {
  f<-score_fit;count<-0L
  extra<-quote({
    n0<-input$design$original_sites
    base_idx<-unlist(lapply(seq_len(S),function(j)((j-1L)*n)+seq_len(n0)))
    base_truth<-as.vector(true_psi)[base_idx]
    add_group('occupancy_original_sites','all',base_truth,psi_draws[base_idx,,,drop=FALSE])
    for(band in c('low','medium','high')) {
      ix0<-which(switch(band,low=base_truth<.2,medium=base_truth>=.2 & base_truth<=.8,high=base_truth>.8))
      if(length(ix0))add_group('occupancy_original_sites',band,base_truth[ix0],psi_draws[base_idx[ix0],,,drop=FALSE])
    }
    prevalence<-colMeans(true_psi)
  })
  visit<-function(e) {
    if(!is.call(e))return(e)
    if(identical(e,quote(prevalence <- colMeans(true_psi)))) {count<<-count+1L;return(extra)}
    for(j in seq_along(e))if(is.call(e[[j]]))e[[j]]<-visit(e[[j]])
    e
  }
  body(f)<-visit(body(f));stopifnot(count==1L);f
}

check_known_u_output <- function(fit,U_true) {
  U<-fit$results_output$jsdm_output$U_output
  flat<-matrix(U,nrow=nrow(U_true));q<-qr(U_true)
  residual<-max(abs(flat-qr.fitted(q,flat)))
  orthogonal<-0
  for(ch in seq_len(dim(U)[4]))for(it in unique(c(1L,ceiling(dim(U)[3]/2),dim(U)[3]))) {
    Q<-qr.coef(q,U[,,it,ch]);orthogonal<-max(orthogonal,max(abs(crossprod(Q)-diag(ncol(U_true)))))
  }
  stopifnot(residual<1e-10,orthogonal<1e-10)
  c(projection=residual,orthogonal=orthogonal)
}
