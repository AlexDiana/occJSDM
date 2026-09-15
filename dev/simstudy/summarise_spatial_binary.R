#!/usr/bin/env Rscript
# Read-only analysis of saved full package fits. Never launches a fit.
.sourced_file <- if(sys.nframe()>0L)sys.frame(1)$ofile else NULL
.analysis_file <- if(!is.null(.sourced_file)) normalizePath(.sourced_file,mustWork=TRUE) else {
  f <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")])
  if(length(f)!=1L)stop("Cannot identify this analysis script")
  normalizePath(f,mustWork=TRUE)
}

trace_summary <- function(x) {
  x <- as.matrix(x)
  stopifnot(all(is.finite(x)))
  constant <- apply(x,2,function(v) length(unique(v))==1L)
  finite_or_na <- function(v) if(length(v)==1L&&is.finite(v)) unname(v) else NA_real_
  notes <- character()
  diagnostic <- function(f) if(any(constant)) NA_real_ else withCallingHandlers(finite_or_na(f(x)),
    warning=function(w){notes<<-c(notes,conditionMessage(w));invokeRestart("muffleWarning")})
  rhat <- if(ncol(x)<2L) NA_real_ else diagnostic(posterior::rhat)
  ess_mean <- diagnostic(posterior::ess_mean)
  ess_bulk <- diagnostic(posterior::ess_bulk)
  mcse <- diagnostic(posterior::mcse_mean)
  notes <- c(if(any(constant)) "Constant chain: Rhat, ESS and MCSE unavailable" else
    if(ncol(x)<2L) "Single-chain timing pilot: Rhat unavailable",unique(notes))
  data.frame(mean=mean(x),sd=sd(as.vector(x)),
    lower=unname(quantile(x,.025)),median=median(x),upper=unname(quantile(x,.975)),
    rhat=rhat,ess_mean=ess_mean,ess_bulk=ess_bulk,mcse=mcse,
    chain_gap=max(colMeans(x))-min(colMeans(x)),
    chain2_minus_chain1=if(ncol(x)==2L) mean(x[,2])-mean(x[,1]) else NA_real_,
    constant_chains=sum(constant),chains=ncol(x),retained=nrow(x),
    diagnostic_note=paste(notes,collapse="; "))
}

probability_groups <- function(p) {
  ans <- list(low=p<.2,medium=p>=.2&p<=.8,high=p>.8,
    bin_0_0.2=p<=.2,bin_0.2_0.4=p>.2&p<=.4,bin_0.4_0.6=p>.4&p<=.6,
    bin_0.6_0.8=p>.6&p<=.8,bin_0.8_1=p>.8,
    rare_below_0.05=p<.05,rare_below_0.1=p<.1,
    rare_above_0.9=p>.9,rare_above_0.95=p>.95)
  for(s in seq_len(ncol(p))) {
    mask <- matrix(FALSE,nrow(p),ncol(p));mask[,s] <- TRUE
    ans[[paste0("species",s)]] <- mask
  }
  ans
}

point_score <- function(x,y) c(bias=mean(x-y),rmse=sqrt(mean((x-y)^2)),
  correlation=cor(as.vector(x),as.vector(y)),
  slope=sum((x-mean(x))*(y-mean(y)))/sum((y-mean(y))^2),
  estimate_sd=sd(as.vector(x)),truth_sd=sd(as.vector(y)),
  sd_ratio=sd(as.vector(x))/sd(as.vector(y)),
  rms_ratio=sqrt(mean(x^2)/mean(y^2)))

analyze_case <- function(root,label,pilot=FALSE,analysis_out) {
  src <- file.path(root,label)
  required <- file.path(src,c("fit.rds","validation.rds","data-truth.rds","provenance.rds"))
  if(!all(file.exists(required))) return(NULL)
  fit <- tryCatch(readRDS(required[1]),error=function(e)NULL)
  saved <- tryCatch(readRDS(required[2]),error=function(e)NULL)
  if(is.null(fit)||is.null(saved)) return(NULL)
  data <- readRDS(required[3]);truth <- data$truth;provenance <- readRDS(required[4])
  js <- fit$results_output$jsdm_output
  n <- nrow(truth$psi);S <- ncol(truth$psi)
  iters <- nrow(js$idx_ls_output);chains <- ncol(js$idx_ls_output)
  if(!pilot) stopifnot(iters==600L,chains==2L,saved$nburn==300L,
                      fit$infos$ps==80L,truth$repeats==60L)
  stopifnot(fit$infos$model=="binary",max(abs(fit$Xs-truth$standardized_coordinates))<1e-12,
            max(abs(fit$X_psi-matrix(truth$standardized_environment,n,1)))<1e-12,
            identical(unname(fit$infos$OTU),unname(data$data$OTU)))
  out <- file.path(analysis_out,label)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  source_file <- .analysis_file
  hashes <- tools::md5sum(c(required,source_file))
  cache <- file.path(out,"analysis.rds")
  if(file.exists(cache)) {
    old <- readRDS(cache)
    if(identical(old$hashes,hashes)) return(old)
  }
  xs <- fit$infos$list_Xs
  index <- xs$Xs_index;first <- match(seq_len(nrow(xs$X_s)),index)
  stopifnot(!anyNA(first),all(table(index)==60L),
            max(abs(fit$Xs-xs$X_s[index,,drop=FALSE]))<1e-12)
  kernel <- function(a,b,l) exp(-(outer(a[,1],b[,1],"-")^2+outer(a[,2],b[,2],"-")^2)/(2*l^2))
  H <- lapply(fit$infos$l_s_grid,function(l) {
    LL <- t(chol(kernel(xs$X_tilde,xs$X_tilde,l)+diag(1e-5,fit$infos$ps)))
    t(solve(LL,t(kernel(xs$X_s,xs$X_tilde,l))))
  })
  # Full-support inputs have the canonical complete coefficient layout.
  stopifnot(identical(unname(xs$Xs_centers),matrix(rep(seq_len(fit$infos$ps),each=n),n,fit$infos$ps)))
  groups <- probability_groups(truth$psi)
  counts <- vapply(groups,sum,numeric(1))
  group_index <- lapply(groups,which)
  mean_truth <- vapply(seq_along(groups),function(g)
    if(counts[g]) mean(truth$psi[group_index[[g]]]) else NA_real_,numeric(1))
  errors <- array(NA_real_,c(iters,chains,length(groups)),dimnames=list(NULL,NULL,names(groups)))
  field_traces <- array(NA_real_,c(iters,chains,4),dimnames=list(NULL,NULL,
    c("field_slope","field_sd_ratio","field_rms_ratio","field_mean_error")))
  probability_mean <- field_mean <- array(0,c(n,S,chains))
  field_truth <- truth$field[first,,drop=FALSE]
  truth_centered <- field_truth-mean(field_truth)
  truth_ss <- sum(truth_centered^2)
  for(ch in seq_len(chains)) for(it in seq_len(iters)) {
    f <- H[[js$idx_ls_output[it,ch]]] %*% matrix(js$Bs_output[,,it,ch],nrow=fit$infos$ps)
    field <- f[index,,drop=FALSE]
    eta <- sweep(field+fit$X_psi %*% matrix(js$B_output[,,it,ch],nrow=ncol(fit$X_psi)),
                 2,js$B0_output[,it,ch],"+")
    if(fit$infos$n_factors>0L)
      eta <- eta+matrix(js$U_output[,,it,ch],nrow=n)%*%matrix(js$L_output[,,it,ch],ncol=S)
    probability <- plogis(eta)
    for(g in seq_along(groups)) if(counts[g])
      errors[it,ch,g] <- mean(probability[group_index[[g]]])-mean_truth[g]
    field_traces[it,ch,] <- c(sum((f-mean(f))*truth_centered)/truth_ss,
      sd(as.vector(f))/sd(as.vector(field_truth)),sqrt(mean(f^2)/mean(field_truth^2)),
      mean(f-field_truth))
    probability_mean[,,ch] <- probability_mean[,,ch]+probability/iters
    field_mean[,,ch] <- field_mean[,,ch]+field/iters
  }
  q <- apply(probability_mean,c(1,2),mean)
  f <- apply(field_mean,c(1,2),mean)
  verification <- c(probability_mean=max(abs(q-saved$psi_mean)),
    probability_chains=max(abs(probability_mean-saved$psi_chain)),
    field_mean=max(abs(f-saved$field_mean)),field_chains=max(abs(field_mean-saved$field_chain)))
  if(!is.null(fit$results_output$psi_output)) {
    stored <- fit$results_output$psi_output
    if(length(dim(stored))==4L) stored <- apply(stored,c(1,2),mean)
    verification <- c(verification,fit_probability=max(abs(q-stored)))
  }
  stopifnot(max(verification)<1e-11)
  strata <- do.call(rbind,lapply(seq_along(groups),function(g) {
    summary <- if(counts[g]) trace_summary(matrix(errors[,,g],iters,chains)) else {
      z <- trace_summary(matrix(0,iters,chains));z[] <- NA;z$diagnostic_note <- "Subset absent";z
    }
    data.frame(case=label,true_range=truth$range,stratum=names(groups)[g],n=counts[g],
               true_probability=mean_truth[g],summary,check.names=FALSE)
  }))
  stopifnot(max(abs(strata$mean[counts>0]-vapply(group_index[counts>0],function(i)mean(q[i]-truth$psi[i]),numeric(1))))<1e-12)
  chain_rows <- do.call(rbind,lapply(seq_along(groups),function(g) do.call(rbind,lapply(seq_len(chains),function(ch) {
    x <- errors[,ch,g]
    data.frame(case=label,stratum=names(groups)[g],n=counts[g],chain=ch,
      mean_error=mean(x),mcse=if(counts[g]) trace_summary(matrix(x,ncol=1))$mcse else NA_real_)
  }))))
  parameter_traces <- list(range=matrix(fit$infos$l_s_grid[js$idx_ls_output],iters,chains),
    spatial_sd=js$sigmabs_output)
  for(term in dimnames(field_traces)[[3]]) parameter_traces[[term]] <- matrix(field_traces[,,term],iters,chains)
  parameters <- do.call(rbind,lapply(names(parameter_traces),function(term)
    data.frame(case=label,parameter=term,truth=if(term=="range")truth$range else
      if(term=="field_mean_error")0 else 1,trace_summary(parameter_traces[[term]]))))
  chain_parameters <- do.call(rbind,lapply(names(parameter_traces),function(term)
    do.call(rbind,lapply(seq_len(chains),function(ch) {
      x <- parameter_traces[[term]][,ch]
      data.frame(case=label,parameter=term,chain=ch,mean=mean(x),median=median(x),
        lower=unname(quantile(x,.025)),upper=unname(quantile(x,.975)),
        transitions=if(term=="range")sum(diff(x)!=0) else NA_integer_)
    }))))
  points <- rbind(data.frame(case=label,quantity="field",chain="combined",t(point_score(f,truth$field))),
                  data.frame(case=label,quantity="probability",chain="combined",t(point_score(q,truth$psi))))
  for(ch in seq_len(chains)) points <- rbind(points,
    data.frame(case=label,quantity="field",chain=as.character(ch),t(point_score(field_mean[,,ch],truth$field))),
    data.frame(case=label,quantity="probability",chain=as.character(ch),t(point_score(probability_mean[,,ch],truth$psi))))
  species <- data.frame(case=label,species=fit$infos$speciesNames,true_prevalence=colMeans(truth$psi),
    observed_prevalence=colMeans(data$data$OTU),fitted_prevalence=colMeans(q))
  result <- list(label=label,pilot=pilot,hashes=hashes,settings=provenance$settings,
    seeds=truth[c("data_seed","fit_seed","environment_seed","outcome_seed")],
    true_range=truth$range,strata=strata,chain_errors=chain_rows,parameters=parameters,
    chain_parameters=chain_parameters,point_metrics=points,species=species,
    verification=verification,errors=errors,parameter_traces=parameter_traces,
    posterior_probability_mean=q,posterior_field_mean=f,
    fit_psi_output_present=!is.null(fit$results_output$psi_output),
    warnings=readLines(file.path(src,"fit-warnings.txt"),warn=FALSE),session=sessionInfo())
  saveRDS(result,cache)
  for(term in c("strata","chain_errors","parameters","chain_parameters","point_metrics","species"))
    write.csv(result[[term]],file.path(out,paste0(term,".csv")),row.names=FALSE)
  write.csv(data.frame(check=names(verification),maximum_difference=verification),file.path(out,"reconstruction.csv"),row.names=FALSE)
  result
}

aggregate_three_cases <- function(results) {
  stopifnot(length(results)==3L,!any(vapply(results,`[[`,logical(1),"pilot")))
  seed <- vapply(results,function(z)z$seeds$data_seed,numeric(1))
  stopifnot(length(unique(seed))==3L,length(unique(vapply(results,`[[`,numeric(1),"true_range")))==3L)
  terms <- results[[1]]$strata$stratum
  do.call(rbind,lapply(terms,function(term) {
    rows <- lapply(results,function(z)z$strata[z$strata$stratum==term,])
    bias <- vapply(rows,function(z)z$mean,numeric(1))
    mcse <- vapply(rows,function(z)z$mcse,numeric(1))
    n_valid <- sum(is.finite(bias))
    data.frame(stratum=term,datasets=n_valid,
      mean_bias=if(n_valid==3L)mean(bias) else NA_real_,
      between_dataset_se=if(n_valid==3L)sd(bias)/sqrt(3) else NA_real_,
      combined_mcse=if(all(is.finite(mcse)))sqrt(sum(mcse^2))/3 else NA_real_,
      minimum_case_bias=if(n_valid)min(bias,na.rm=TRUE) else NA_real_,
      maximum_case_bias=if(n_valid)max(bias,na.rm=TRUE) else NA_real_,
      maximum_absolute_case_bias=if(n_valid)max(abs(bias),na.rm=TRUE) else NA_real_)
  }))
}

number <- function(x,d=4) ifelse(is.finite(x),formatC(x,format="f",digits=d),"unavailable")
report_case <- function(z) {
  lines <- c(paste0("## ",z$label),"",
    if(z$pilot) "Timing pilot only. This run is excluded from the three-dataset summary and cannot establish parameter recovery or convergence." else
      paste0("Independent dataset at true range ",number(z$true_range),"; two chains, 300 burn-in and 600 retained draws per chain."),"")
  for(i in which(!grepl("^species",z$strata$stratum))) {
    r <- z$strata[i,]
    lines <- c(lines,paste0("- ",r$stratum,": ",r$n," species-site entries; mean probability error ",number(100*r$mean,3),
      " percentage points; MCSE ",number(100*r$mcse,3)," points; R-hat ",number(r$rhat,3),
      "; ESS for the mean ",number(r$ess_mean,0),"; chain 2 minus chain 1 ",number(100*r$chain2_minus_chain1,3)," points.",
      if(nzchar(r$diagnostic_note))paste0(" ",r$diagnostic_note) else ""))
  }
  fp <- z$point_metrics[z$point_metrics$quantity=="field"&z$point_metrics$chain=="combined",]
  qp <- z$point_metrics[z$point_metrics$quantity=="probability"&z$point_metrics$chain=="combined",]
  lines <- c(lines,"",paste0("Posterior mean field slope is ",number(fp$slope),", field SD ratio is ",number(fp$sd_ratio),
    ", and field RMSE is ",number(fp$rmse),". Probability RMSE is ",number(qp$rmse),". These evaluate the posterior mean field; the draw-level field functionals in the CSV describe a different quantity."),"")
  for(i in seq_len(nrow(z$parameters))) {
    r <- z$parameters[i,]
    lines <- c(lines,paste0("- ",r$parameter,": posterior mean ",number(r$mean)," (truth ",number(r$truth),
      "); MCSE ",number(r$mcse),"; R-hat ",number(r$rhat,3),"; ESS for mean ",number(r$ess_mean,0),
      "; chain mean gap ",number(r$chain_gap),if(nzchar(r$diagnostic_note))paste0(". ",r$diagnostic_note) else "."))
  }
  lines <- c(lines,"",paste0("Largest reconstruction discrepancy is ",format(max(z$verification),scientific=TRUE),
    ". The comparison includes saved probability means and chain-specific means. ",
    if(z$fit_psi_output_present) "The fit's own probability output was also checked." else "This binary fit has no separate psi_output, so the saved validation means were the comparison target."),
    "",paste0("Species true mean prevalence ranges from ",number(min(z$species$true_prevalence),3)," to ",number(max(z$species$true_prevalence),3),
    ". Rare-probability subsets do not substitute for rare-species validation."),"")
  if(length(z$warnings)) lines <- c(lines,"Saved fitting warnings:","",paste0("- ",z$warnings),"")
  lines
}

write_report <- function(root,results,pilot=FALSE,analysis_out) {
  out <- analysis_out;dir.create(out,showWarnings=FALSE)
  path <- file.path(out,if(pilot)"timing-pilot-report.md" else "three-dataset-report.md")
  lines <- c(if(pilot)"# Binary timing-pilot analysis" else "# Binary package analysis across three independent datasets","",
    if(pilot) "This report develops and verifies the analysis workflow only. It is excluded from all full-run summaries." else
      "This analysis contains n = 3 independent datasets at three different generating ranges, with one dataset per range. It is not three simulation replicates at each range. Each dataset receives equal weight in the reported mean bias, regardless of its number of observations in a stratum.","",
    "Probability errors are reconstructed for every posterior draw, after applying the logistic transformation, and averaged separately in the predeclared low, medium and high true-probability groups. The five equal-width bins and rare-probability subsets are additional disaggregations. The provisional absolute signed-error target remains five percentage points. No release or general bias pass is inferred from one dataset, a pooled species-site count or a near-zero overall signed mean.","",
    "Monte Carlo SE describes numerical uncertainty from the retained MCMC draws. The between-dataset SE is separately computed as the sample standard deviation of three dataset-level errors divided by sqrt(3). Because ranges differ and n is only three, it includes between-range heterogeneity and cannot quantify uncertainty at any individual range or establish general frequentist calibration.","",
    "Diagnostics use posterior's [rank-normalized R-hat](https://mc-stan.org/posterior/reference/rhat.html), ESS for the posterior mean, bulk ESS, and [Monte Carlo SE for the mean](https://mc-stan.org/posterior/reference/mcse_mean.html). Constant chains produce unavailable diagnostics, not an automatic convergence conclusion. R-hat is suppressed for the single-chain timing pilot.","")
  if(!pilot) {
    pooled <- aggregate_three_cases(results)
    write.csv(pooled,file.path(out,"three-dataset-summary.csv"),row.names=FALSE)
    saveRDS(list(cases=vapply(results,`[[`,character(1),"label"),summary=pooled,
      description="Equal weights; n=3 independently seeded datasets across three ranges, not three replicates per range"),file.path(out,"three-dataset-summary.rds"))
    lines <- c(lines,"## Equal-weight summary across the three datasets","")
    for(i in which(!grepl("^species",pooled$stratum))) {
      r <- pooled[i,]
      lines <- c(lines,paste0("- ",r$stratum,": ",r$datasets," contributing datasets; mean bias ",number(100*r$mean_bias,3),
        " percentage points; between-dataset SE ",number(100*r$between_dataset_se,3)," points; combined MCSE ",
        number(100*r$combined_mcse,3)," points; individual dataset errors range from ",number(100*r$minimum_case_bias,3),
        " to ",number(100*r$maximum_case_bias,3)," points."))
    }
    lines <- c(lines,"","A three-dataset mean is omitted if a subset is absent from any dataset; cases are not silently reweighted. Per-case and per-chain values remain available in the CSV files.","")
  }
  for(z in results) lines <- c(lines,report_case(z))
  lines <- c(lines,"## Reproduction","","Run `Rscript dev/simstudy/summarise_spatial_binary.R --base=INPUT_DIRECTORY --out=NEW_DIRECTORY --pilot` for the timing pilot, or omit `--pilot` for the three full fits. Use explicit `--resume` only to continue an existing analysis output directory. The script reads only completed saved fits and validation objects, records input and script hashes, and never launches fitting or simulation. Detailed CSV files and draw-level traces are stored beside this report in one analysis subdirectory per fit.","")
  writeLines(lines,path)
  path
}

main <- function() {
  args <- commandArgs(TRUE)
  if("--help" %in% args) {
    cat("--base=INPUT_DIRECTORY --out=NEW_DIRECTORY [--pilot] [--resume]\nThe base contains package60-range{4,6,8}-full or package60-timing-pilot. --root is accepted as an alias for --base.\n")
    return(invisible(NULL))
  }
  roots <- sub("^--(root|base)=","",args[grepl("^--(root|base)=",args)])
  outputs <- sub("^--out=","",args[startsWith(args,"--out=")])
  if(length(roots)!=1L||length(outputs)!=1L)stop("Supply one --base=INPUT_DIRECTORY and one --out=NEW_DIRECTORY")
  root <- normalizePath(roots[1],mustWork=TRUE)
  analysis_out <- outputs[1]
  if(dir.exists(analysis_out)&&length(list.files(analysis_out,all.files=TRUE,no..=TRUE))&&!"--resume"%in%args)
    stop("Output directory is not empty; use a new directory or explicit --resume: ",analysis_out)
  dir.create(analysis_out,recursive=TRUE,showWarnings=FALSE)
  analysis_out <- normalizePath(analysis_out,mustWork=TRUE)
  pilot <- "--pilot" %in% args
  labels <- if(pilot)"package60-timing-pilot" else paste0("package60-range",c(4,6,8),"-full")
  results <- lapply(labels,function(label)analyze_case(root,label,pilot,analysis_out))
  ready <- !vapply(results,is.null,logical(1))
  cat("Analyzed:",paste(labels[ready],collapse=", "),"\n")
  if(!all(ready)) {
    cat("Pending completed fit.rds and validation.rds:",paste(labels[!ready],collapse=", "),"\n")
    return(invisible(NULL))
  }
  report <- write_report(root,results,pilot,analysis_out)
  cat("Report:",report,"\n")
  invisible(results)
}
if(sys.nframe()==0L)main()
