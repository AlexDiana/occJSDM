# Independent checks; run from the repository root. Does not source producers or rerun MCMC.
# Rscript prediction-verify.R ORIGINAL_ARCHIVE NEW_ARCHIVE [--score-mcse]
# The optional pass checks every posterior mean/interval and probability diagnostic,
# and estimates score Monte Carlo error using all 24000 saved draws per fit.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 2:3)
full <- length(args) == 3L && identical(args[3L], "--score-mcse")
if (length(args) == 3L) stopifnot(full)
archive <- normalizePath(args[1L], mustWork=TRUE)
new_archive <- normalizePath(args[2L], mustWork=TRUE)
.libPaths(c(file.path(archive,"library"), .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(identical(normalizePath(find.package("occJSDM")),
                    normalizePath(file.path(archive,"library/occJSDM"))))
md5 <- function(path) unname(tools::md5sum(path))
near <- function(a,b,tol=1e-11,label="numeric check") {
  if(length(a)!=length(b) || anyNA(a) || anyNA(b) ||
     any(!is.finite(a)) || any(!is.finite(b)) ||
     max(c(0,abs(as.numeric(a)-as.numeric(b)))) > tol) stop(label," failed")
}
key <- function(x) paste(x$Site,x$species,sep="|")
op <- file.path(archive,"input.rds")
ip <- file.path(new_archive,"prediction-input.rds")
o <- readRDS(op); input <- readRDS(ip)
b <- readRDS("vignettes/teaching-data/prediction-lesson.rds")
old_lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
paths <- c(two_factors=file.path(archive,"default-fit.rds"),
           one_factor=file.path(new_archive,"one-factor-fit.rds"))
saved <- lapply(paths,readRDS); fits <- lapply(saved,`[[`,"fit")
stopifnot(identical(input,b$input),identical(o,old_lesson$input),
  identical(b$prediction_input_md5,md5(ip)),identical(input$original_input_md5,md5(op)),
  identical(input$original_lesson_md5,md5("vignettes/teaching-data/nonspatial-lesson.rds")),
  identical(input$baseline_fit_md5,md5(paths[["two_factors"]])),
  identical(input$generator_md5,md5("dev/simstudy/vignette-lesson/prediction-build.R")),
  identical(b$source_md5,tools::md5sum(names(b$source_md5))),
  identical(o$source_hashes,b$source_hashes),identical(input$source_hashes,b$source_hashes),
  identical(b$source_hashes,tools::md5sum(names(b$source_hashes))),
  identical(saved$two_factors$mcmc,saved$one_factor$mcmc),
  identical(saved$two_factors$priors,saved$one_factor$priors),
  identical(saved$one_factor$prediction_input_md5,md5(ip)))
for(arm in names(saved)) {
  r <- saved[[arm]]; f <- r$fit; m <- b$manifests[[arm]]
  stopifnot(identical(r$source_hashes,b$source_hashes),identical(r$input_md5,md5(op)),
    identical(m$md5,md5(paths[[arm]])),identical(m$seed,r$seed),identical(m$mcmc,r$mcmc),
    identical(m$priors,r$priors),identical(m$warnings,r$warnings),identical(m$session,r$session),
    identical(m$factors,f$infos$n_factors),identical(m$waic,f$results_output$WAIC),
    identical(input$scaling,f$infos$list_X_psi_mat),f$infos$model=="two_stage",
    f$infos$jsdmModel=="binary",f$infos$ps==0,
    identical(dim(f$results_output$jsdm_output$B0_output),c(10L,6000L,4L)))
  dat <- o$sim$data_list
  rr <- order(dat$info$Site,dat$info$Sample,dat$info$Primer)
  ss <- match(f$infos$speciesNames,colnames(dat$OTU))
  stopifnot(!anyNA(ss),isTRUE(all.equal(f$infos$data_info[,names(dat$info),drop=FALSE],
                                     dat$info[rr,,drop=FALSE],check.attributes=FALSE)))
  near(f$infos$OTU,dat$OTU[rr,ss,drop=FALSE],0,paste(arm,"original PCR reads"))
}
stopifnot(fits$one_factor$infos$n_factors==1L,fits$two_factors$infos$n_factors==2L)
for(nm in c("OTU","data_info","list_idx","siteNames","speciesNames","primerNames",
             "list_X_psi_mat","list_X_theta_mat","list_Tr_mat"))
  stopifnot(identical(fits$one_factor$infos[[nm]],fits$two_factors$infos[[nm]]))
for(nm in c("X_psi","X_theta","Xs","Tr"))
  stopifnot(identical(fits$one_factor[[nm]],fits$two_factors[[nm]]))
cat("PASS: hashes, manifests, complete posterior dimensions, shared training arrays\n")

# Regenerate genuinely new sites from fixed community parameters and frozen training scale.
stopifnot(input$new_site_seed==20260925L,input$fitting_seed==20260926L,
  input$public_prediction_seed==20260927L,identical(input$selected_sites,as.character(101:110)))
truth <- o$sim$true_params$jsdmParams_true
training_info <- o$sim$data_list$info
train <- training_info[!duplicated(training_info$Site),input$scaling$names_df,drop=FALSE]
centers <- vapply(train,mean,numeric(1)); scales <- vapply(train,sd,numeric(1))
near(centers,input$scaling$mean_df,0,"training centers")
near(scales,input$scaling$sd_df,0,"training scales")
near(sweep(fits$two_factors$X_psi%*%truth$B+truth$U%*%truth$L,2,truth$B0,"+"),
     truth$eta,1e-12,"original truth reconstruction")
set.seed(input$new_site_seed)
raw <- matrix(rnorm(600L,mean=0,sd=10),300L,2L,
              dimnames=list(as.character(101:400),input$scaling$names_df))
hidden <- matrix(rnorm(600L,sd=o$jsdm$sigma_h),300L,2L)
x <- sweep(sweep(raw,2,centers,"-"),2,scales,"/")
eta <- sweep(x%*%truth$B,2,truth$B0,"+")
conditional <- plogis(eta+hidden%*%truth$L)
z <- matrix(rbinom(length(conditional),1,as.vector(conditional)),300L,10L)
dimnames(conditional) <- dimnames(z) <- list(rownames(raw),fits$two_factors$infos$speciesNames)
stopifnot(identical(raw,input$raw_covariates),identical(hidden,input$hidden),
  identical(x,input$standardized),identical(eta,input$environmental_eta),
  identical(conditional,input$conditional_probability),identical(z,input$actual_presence),
  !any(rownames(raw)%in%as.character(fits$two_factors$infos$siteNames)))
outside <- apply(sapply(seq_len(ncol(raw)),function(k)
  raw[,k]<min(train[[k]]) | raw[,k]>max(train[[k]])),1,any)
stopifnot(identical(unname(outside),unname(input$outside_training_range)))

# Independent adaptive integration, explicitly including the standard-normal density.
integrated_mean <- function(a,s) integrate(function(v)
  vapply(v,function(u)mean(plogis(a+s*u)),numeric(1))*dnorm(v),-Inf,Inf,
  rel.tol=1e-11,abs.tol=1e-11,subdivisions=1000L,stop.on.error=TRUE)$value
true_sd <- o$jsdm$sigma_h*sqrt(colSums(truth$L^2))
true_mean <- matrix(NA_real_,300L,10L,dimnames=dimnames(z))
for(s in seq_len(10L)) for(i in seq_len(300L))
  true_mean[i,s] <- integrated_mean(eta[i,s],true_sd[s])
tkeys <- expand.grid(Site=rownames(z),species=colnames(z),stringsAsFactors=FALSE)
ti <- match(key(tkeys),key(b$truth))
stopifnot(!anyNA(ti),!anyDuplicated(key(b$truth)),nrow(b$truth)==3000L)
near(b$truth$truth[ti],true_mean,2e-10,"all marginal truth integrals")
near(b$truth$conditional_truth[ti],conditional,0,"conditional truth")
near(b$truth$z[ti],z,0,"new occupancy realizations")
cat("PASS: independent-site seed, frozen scaling, all 3000 marginal truth integrals\n")

# Separate quadrature implementation for posterior quantiles; adaptive means above
# and below provide numerical references that do not use Gaussian quadrature.
reference_rule <- function(n) {
  j <- matrix(0,n,n); j[row(j)==col(j)+1L] <- sqrt(seq_len(n-1L))
  e <- eigen(j+t(j),symmetric=TRUE)
  list(x=e$values,w=e$vectors[1L,]^2)
}
integrate_draws <- function(a,s,rule) {
  answer <- numeric(length(a))
  for(k in seq_along(rule$x)) answer <- answer+rule$w[k]*plogis(a+s*rule$x[k])
  answer
}
ref <- reference_rule(121L)
near(sum(ref$w),1,1e-12,"quadrature normalization")
near(sum(ref$w*ref$x^2),1,1e-12,"quadrature variance")
extract <- function(f,species) {
  j <- f$results_output$jsdm_output; ni <- dim(j$B0_output)[2]; nc <- dim(j$B0_output)[3]
  # Assemble by chain, not with producer's array flattening; retain every iteration.
  a <- unlist(lapply(seq_len(nc),function(ch)j$B0_output[species,,ch]),use.names=FALSE)
  slopes <- do.call(cbind,lapply(seq_len(nc),function(ch)
    matrix(j$B_output[,species,,ch],nrow=ncol(f$X_psi))))
  s <- unlist(lapply(seq_len(nc),function(ch) {
    loads <- matrix(j$L_output[,species,,ch],nrow=f$infos$n_factors)
    j$sigmah_output[,ch]*sqrt(colSums(loads^2))
  }),use.names=FALSE)
  list(a=a,slopes=slopes,s=s,ni=ni,nc=nc)
}
max_mean_diff <- max_interval_diff <- max_draw_diff <- 0
score_mcse <- list()
for(arm in names(fits)) {
  f <- fits[[arm]]; params <- lapply(seq_len(10L),function(s)extract(f,s))
  stopifnot(all(vapply(params,function(p)length(p$a)==24000L,logical(1))))
  extreme_species <- which.max(vapply(params,function(p)max(p$s),numeric(1)))
  # Cases chosen by indices/residual SD, never by observed outcomes or scores.
  for(s in unique(c(1L,5L,10L,extreme_species))) {
    p <- params[[s]]
    for(i in c(1L,150L,300L)) {
      a <- as.vector(x[i,,drop=FALSE]%*%p$slopes)+p$a
      m <- integrated_mean(a,p$s); draws <- integrate_draws(a,p$s,ref)
      near(mean(draws),m,5e-9,"high-order versus adaptive all-draw mean")
      row <- b$cells[b$cells$arm==arm & b$cells$Site==rownames(x)[i] &
                       b$cells$species==f$infos$speciesNames[s],]
      stopifnot(nrow(row)==1L); q <- quantile(draws,c(.025,.975),names=FALSE)
      near(row$estimate,m,1e-7,"selected-cell exported mean")
      near(c(row$lower,row$upper),q,1e-6,"selected-cell marginal intervals")
      max_mean_diff <- max(max_mean_diff,abs(row$estimate-m))
      max_interval_diff <- max(max_interval_diff,abs(c(row$lower,row$upper)-q))
      for(tt in unique(c(1L,6000L,6001L,12000L,18001L,24000L,which.max(p$s)))) {
        adaptive <- integrated_mean(a[tt],p$s[tt])
        near(adaptive,draws[tt],1e-8,"individual high-SD draw")
        max_draw_diff <- max(max_draw_diff,abs(adaptive-draws[tt]))
      }
    }
  }
  cat("PASS:",arm,"selected-cell all-draw means/quantiles and adaptive high-SD checks\n")
  set.seed(input$public_prediction_seed)
  native <- predictNewSites(f,X_psi=as.data.frame(raw[input$selected_sites,,drop=FALSE]),
                           useSpatial=FALSE,confidence=.95,verbose=FALSE)
  stopifnot(identical(dim(native),c(3L,10L,10L)))
  nk <- expand.grid(Site=input$selected_sites,species=f$infos$speciesNames,stringsAsFactors=FALSE)
  public <- b$public[b$public$arm==arm,]; ii <- match(key(nk),key(public))
  stopifnot(!anyNA(ii),!anyDuplicated(key(public)),nrow(public)==100L)
  near(public$lower[ii],native[1,,],1e-12,"public lower")
  near(public$median[ii],native[2,,],1e-12,"public median")
  near(public$upper[ii],native[3,,],1e-12,"public upper")
  near(public$conditional_truth[ii],conditional[input$selected_sites,],0,"conditional interval truth")
  stopifnot(all(public$lower<=public$median),all(public$median<=public$upper))
  reported <- b$diagnostics[b$diagnostics$arm==arm,]; reported$arm <- NULL
  current <- as.data.frame(returnConvergenceDiagnostics(f))
  rownames(reported) <- rownames(current) <- NULL
  stopifnot(isTRUE(all.equal(reported,current,tolerance=1e-12)))
  cat("PASS:",arm,"seeded public prediction and parameter diagnostics\n")
  if(full) {
    rule <- reference_rule(b$quadrature[[arm]]$nodes)
    ib <- il <- numeric(24000L)
    for(s in seq_len(10L)) {
      p <- params[[s]]; avg <- numeric(24000L)
      rows <- b$cells[b$cells$arm==arm & b$cells$species==f$infos$speciesNames[s],]
      rows <- rows[match(rownames(x),rows$Site),]
      stopifnot(nrow(rows)==300L,!anyNA(rows$Site))
      # Score derivatives at posterior means; preserve joint draw dependence.
      gb <- 2*(rows$estimate-rows$z)/3000
      gl <- (rows$estimate-rows$z)/(rows$estimate*(1-rows$estimate))/3000
      for(i in seq_len(300L)) {
        a <- as.vector(x[i,,drop=FALSE]%*%p$slopes)+p$a
        draws <- integrate_draws(a,p$s,rule)
        near(mean(draws),rows$estimate[i],1e-11,"full all-draw mean")
        near(quantile(draws,c(.025,.975),names=FALSE),c(rows$lower[i],rows$upper[i]),
             1e-11,"full all-draw quantiles")
        avg <- avg+draws/300
        ib <- ib+gb[i]*(draws-rows$estimate[i]); il <- il+gl[i]*(draws-rows$estimate[i])
      }
      dg <- b$probability_diagnostics[b$probability_diagnostics$arm==arm &
                                      b$probability_diagnostics$species==f$infos$speciesNames[s],]
      series <- matrix(avg,p$ni,p$nc)
      near(c(dg$rhat,dg$ess,dg$mcse),
           c(posterior::rhat(series),posterior::ess_mean(series),posterior::mcse_mean(series)),
           1e-6,"probability Monte Carlo diagnostics")
      cat("PASS:",arm,f$infos$speciesNames[s],"all 300 sites, 24000 draws, diagnostics\n")
    }
    score_mcse[[arm]] <- data.frame(arm=arm,metric=c("brier","negative_log_score"),
      mcse=c(posterior::mcse_mean(matrix(ib,6000L,4L)),posterior::mcse_mean(matrix(il,6000L,4L))))
  }
}
# Independently recompute score columns, summaries and uncertainty at the site level.
stopifnot(nrow(b$cells)==6000L,!anyDuplicated(paste(b$cells$arm,key(b$cells))),
          all(b$cells$estimate>0 & b$cells$estimate<1))
for(arm in names(fits)) {
  cells <- b$cells[b$cells$arm==arm,]; ii <- match(key(cells),key(b$truth))
  near(cells$truth,b$truth$truth[ii],0,"cell truth alignment")
  near(cells$z,b$truth$z[ii],0,"cell state alignment")
  error <- cells$estimate-cells$truth; brier <- (cells$estimate-cells$z)^2
  ls <- -dbinom(cells$z,1L,cells$estimate,log=TRUE)
  near(cells$error,error,0,"cell error"); near(cells$absolute_error,abs(error),0,"absolute error")
  near(cells$brier,brier,0,"Brier score"); near(cells$negative_log_score,ls,1e-14,"log score")
  row <- b$scores[b$scores$arm==arm,]
  near(unlist(row[c("cells","signed_error_pp","mean_absolute_error_pp","rmse_pp",
                    "brier","negative_log_score")]),
       c(3000L,mean(error)*100,mean(abs(error))*100,sqrt(mean(error^2))*100,
         mean(brier),mean(ls)),1e-12,"overall score summaries")
}
for(metric in c("brier","negative_log_score")) {
  groups <- lapply(names(fits),function(arm) {
    cells <- b$cells[b$cells$arm==arm,]; tapply(cells[[metric]],cells$Site,mean)
  }); names(groups) <- names(fits)
  stopifnot(identical(names(groups$one_factor),names(groups$two_factors)))
  delta <- groups$one_factor-groups$two_factors; se <- sd(delta)/sqrt(length(delta))
  row <- b$paired_site_differences[b$paired_site_differences$metric==metric,]
  near(unlist(row[c("sites","difference","site_SE","lower","upper")]),
       c(300L,mean(delta),se,mean(delta)+c(-1,1)*1.96*se),1e-12,"paired site uncertainty")
}
cat("PASS: independently recomputed scores and paired site SEs\n")
cat("Maximum selected-cell mean difference:",max_mean_diff,"\n")
cat("Maximum selected-cell interval difference:",max_interval_diff,"\n")
cat("Maximum individual-draw quadrature/adaptive difference:",max_draw_diff,"\n")
cat("New sites outside at least one training covariate range:",sum(outside),"of 300\n")
print(b$scores); print(b$paired_site_differences)
if(full) {
  mcse <- do.call(rbind,score_mcse); print(mcse)
  diffs <- b$paired_site_differences
  diffs$MCSE <- vapply(diffs$metric,function(m)sqrt(sum(mcse$mcse[mcse$metric==m]^2)),numeric(1))
  diffs$difference_over_MCSE <- diffs$difference/diffs$MCSE
  cat("First-order Monte Carlo uncertainty of score differences for independent fits:\n")
  print(diffs)
}
cat("Interpretation: site SE conditions on the saved probability estimates. It excludes\n",
    "MCMC error, uncertainty across training datasets, and uncertainty across communities.\n",
    "Optional score MCSE is a first-order delta approximation, not a refit or formal\n",
    "model-selection decision. Marginal scores do not evaluate joint dependence.\n",
    "Current augmented WAIC is reproduced only as extraction, never ranked.\n",sep="")

# Inspect displayed calls without executing either fit or its prediction again.
# Replace only named API calls in parsed chunk expressions with recording stubs.
read_chunk <- function(label) {
  lines <- readLines("vignettes/occJSDM-lesson-3.Rmd",warn=FALSE)
  first <- grep(paste0("^```\\{r ",label,"[,}]"),lines)
  stopifnot(length(first)==1L)
  last <- first+which(lines[(first+1L):length(lines)]=="```")[1L]
  parse(text=lines[seq.int(first+1L,last-1L)])
}
replace_api <- function(expr, api, replacement) {
  # These inspected chunks contain literal namespaced API calls, not generated code.
  # Reparse only this call-name substitution; NULL and missing subscript arguments
  # must remain intact when capturing the displayed calls.
  parse(text=gsub(paste0("occJSDM::",api),replacement,
                  paste(deparse(expr),collapse="\n"),fixed=TRUE))[[1L]]
}

call_env <- new.env(parent=globalenv())
call_env$lesson <- old_lesson; call_env$prediction_examples <- b
call_env$new_habitat <- as.data.frame(raw); call_env$fitmodel <- fits$two_factors
call_env$tibble <- function(...) data.frame(...,check.names=FALSE)
call_env$capture_fit <- function(...) {
  call_env$fit_call <- list(...); call_env$fit_rng <- .Random.seed; NULL
}
call_env$capture_prediction <- function(...) {
  call_env$prediction_call <- list(...); call_env$prediction_rng <- .Random.seed
  array(0,c(3L,length(input$selected_sites),10L))
}
for(expr in read_chunk("prediction-fit-one-factor"))
  eval(replace_api(expr,"runOccJSDM","capture_fit"),envir=call_env)
for(expr in read_chunk("prediction-native-call"))
  eval(replace_api(expr,"predictNewSites","capture_prediction"),envir=call_env)
expected_fit <- list(data=o$sim$data_list,occCovariates=colnames(raw),collCovariates="X_theta",
  spatCovariates=NULL,threshold=1,listParams=list(n_factors=1,n_lattrait=1),
  MCMCparams=saved$two_factors$mcmc,listPriors=saved$two_factors$priors,
  summarisedLatentPresences=TRUE)
stopifnot(identical(call_env$fit_call,expected_fit))
set.seed(input$fitting_seed); stopifnot(identical(call_env$fit_rng,.Random.seed))
expected_prediction <- list(fits$two_factors,X_psi=as.data.frame(raw[input$selected_sites,,drop=FALSE]),
                            useSpatial=FALSE,confidence=.95,verbose=FALSE)
stopifnot(identical(call_env$prediction_call,expected_prediction))
set.seed(input$public_prediction_seed); stopifnot(identical(call_env$prediction_rng,.Random.seed))
cat("PASS: displayed fit/prediction chunk arguments and seeds captured without fitting\n")
cat("All prediction verification checks passed.\n")
