finite_mean <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
field_metrics <- function(estimate,truth,link,require_matching_link=FALSE) {
  compatible <- identical(link,"logit")
  if(!compatible && require_matching_link) stop("Native field recovery requires matching logit units")
  constant <- function(x) diff(range(x))<1e-12
  status <- if(!compatible) "incompatible link scale" else if(constant(estimate)) "constant fitted field" else if(constant(truth)) "constant true field" else "available"
  zero <- sqrt(mean(truth^2)); rmse <- sqrt(mean((estimate-truth)^2))
  list(rmse=if(compatible) rmse else NA_real_,zero_rmse=zero,
    rmse_reduction=if(compatible) zero-rmse else NA_real_,
    fractional_rmse_reduction=if(compatible && zero>0) 1-rmse/zero else NA_real_,
    correlation=if(status=="available") cor(estimate,truth) else NA_real_,correlation_status=status,
    attenuation=if(compatible && sum(truth^2)>0) sum(estimate*truth)/sum(truth^2) else NA_real_,
    realised_sd=sd(estimate),truth_realised_sd=sd(truth),native_link=link)
}
probability_metrics <- function(estimate,truth,states) {
  n <- length(estimate)
  if(!n) return(list(n=0L,status="empty",bias=NA_real_,mae=NA_real_,rmse=NA_real_,brier=NA_real_,log_score=NA_real_,clamps=0L))
  stopifnot(length(truth)==n,length(states)==n,all(is.finite(c(estimate,truth,states))))
  q <- pmin(pmax(estimate,1e-12),1-1e-12)
  list(n=n,status="available",bias=mean(estimate-truth),mae=mean(abs(estimate-truth)),rmse=sqrt(mean((estimate-truth)^2)),
    brier=mean((estimate-states)^2),log_score=mean(states*log(q)+(1-states)*log1p(-q)),clamps=sum(q!=estimate))
}
truth_view <- function(truth,scope=c("training","lattice")) {
  scope <- match.arg(scope); land <- truth$landscape
  idx <- if(scope=="training") truth$training_index else truth$lattice_index
  ids <- if(scope=="training") truth$training_site_id else truth$prediction_site_id
  env <- land$environment[idx] %o% land$B
  raw <- land$field[idx,,drop=FALSE]
  terms <- center_terms(land$B0,env,raw,colMeans(land$environment[truth$training_index] %o% land$B),colMeans(land$field[truth$training_index,,drop=FALSE]))
  p <- point_prediction(list(site_id=ids,species_id=truth$species_id),terms,plogis,"logit","generating realised landscape")
  stopifnot(max(abs(p$p_full-land$psi[idx,,drop=FALSE]))<1e-12)
  list(site_id=ids,species_id=truth$species_id,probabilities=p,field=terms$spatial,states=land$z[idx,,drop=FALSE],
    xy=land$points[idx,,drop=FALSE],training_xy=land$points[truth$training_index,,drop=FALSE],
    training_states=land$z[truth$training_index,,drop=FALSE],prevalence=land$prevalence,generating_sd=land$field_sd)
}
score_predictions <- function(prediction_record,truth,training_reference=NULL) {
  p <- prediction_record; t <- truth
  stopifnot(identical(p$site_id,t$site_id),identical(p$species_id,t$species_id))
  validate_prediction(p,list(site_id=t$site_id,species_id=t$species_id))
  nearest <- apply(xy_distances(t$xy,t$training_xy),1,min)
  distance <- cut(nearest,breaks=c(-Inf,.02,.05,.1,Inf),labels=c("<=0.02","(0.02,0.05]","(0.05,0.1]",">0.1"),right=TRUE)
  rows <- bins <- list()
  for(j in seq_along(p$species_id)) {
    truthp <- t$probabilities$p_full[,j]; delta <- t$probabilities$delta_p[,j]
    f <- field_metrics(p$field_mean[,j],t$field[,j],p$link)
    d <- field_metrics(p$delta_p[,j],delta,"logit"); d$native_link <- "probability contrast"
    names(f) <- paste0("field_",names(f)); names(d) <- paste0("delta_",names(d))
    pm <- probability_metrics(p$p_full[,j],truthp,t$states[,j])
    spacezero <- probability_metrics(p$p_space_zero[,j],truthp,t$states[,j])
    intervals <- p$intervals; compatible <- p$link=="logit"
    interval <- list(interval_status=p$uncertainty$status,interval_reason=p$uncertainty$reason,
      probability_inclusion=if(!is.null(intervals)) mean(truthp>=intervals$probability_lower[,j] & truthp<=intervals$probability_upper[,j]) else NA_real_,
      probability_width=if(!is.null(intervals)) mean(intervals$probability_upper[,j]-intervals$probability_lower[,j]) else NA_real_,
      field_inclusion=if(!is.null(intervals)&&compatible) mean(t$field[,j]>=intervals$field_lower[,j] & t$field[,j]<=intervals$field_upper[,j]) else NA_real_,
      field_width_native=if(!is.null(intervals)) mean(intervals$field_upper[,j]-intervals$field_lower[,j]) else NA_real_)
    baseline <- if(is.null(training_reference)) NA_real_ else sqrt(mean((training_reference$p_full[,j]-truthp)^2))
    rows[[j]] <- as.data.frame(c(list(species_id=p$species_id[j],prevalence=t$prevalence[j],training_occurrences=sum(t$training_states[,j]),
      generating_spatial_sd=t$generating_sd,delta_rms=sqrt(mean(p$delta_p[,j]^2)),truth_delta_rms=sqrt(mean(delta^2)),
      historical_field_rmse=if(!is.null(p$historical_field_median_centered) && compatible) sqrt(mean((p$historical_field_median_centered[,j]-t$field[,j])^2)) else NA_real_,
      space_zero_rmse=spacezero$rmse,space_zero_benefit=spacezero$rmse-pm$rmse,env_fit_rmse=baseline,env_fit_benefit=baseline-pm$rmse),f,d,pm,interval),stringsAsFactors=FALSE)
    probability_bin <- cut(truthp,breaks=c(-Inf,.2,.8,Inf),labels=c("<0.2","[0.2,0.8]",">0.8"),right=TRUE)
    probability_bin[truthp==.2] <- "[0.2,0.8]"
    for(type in c("distance","truth_probability")) {
      grouping <- if(type=="distance") distance else probability_bin
      for(bin in levels(grouping)) for(target in c("p_full","p_space_zero","env_fit")) {
        keep <- which(grouping==bin)
        values <- if(target=="env_fit") { if(is.null(training_reference)) NULL else training_reference$p_full[,j] } else p[[target]][,j]
        metrics <- if(is.null(values)) { z <- probability_metrics(numeric(),numeric(),numeric()); z$status <- "unavailable env fit"; z$n <- length(keep); z } else probability_metrics(values[keep],truthp[keep],t$states[keep,j])
        bins[[length(bins)+1L]] <- as.data.frame(c(list(species_id=p$species_id[j],prevalence=t$prevalence[j],stratification=type,bin=bin,target=target),metrics),stringsAsFactors=FALSE)
      }
    }
  }
  list(species=do.call(rbind,rows),strata=do.call(rbind,bins),sites=list(
    observed=attribute_sites(p,t$states),expected=attribute_sites(p,t$probabilities$p_full,expected=TRUE)))
}
equal_community_summary <- function(x,metric,expected_communities=3L) {
  by_community <- vapply(split(x[[metric]],x$community),finite_mean,numeric(1)); available <- by_community[is.finite(by_community)]
  data.frame(metric=metric,mean=finite_mean(available),minimum=if(length(available)) min(available) else NA_real_,
    maximum=if(length(available)) max(available) else NA_real_,communities_available=length(available),communities_expected=expected_communities,
    rows_available=sum(is.finite(x[[metric]])),rows_expected=nrow(x))
}
compare_oracle <- function(scores,oracle_scores) {
  keys <- intersect(c("community","arrangement","scope","species_id"),names(scores))
  merged <- merge(scores,oracle_scores,by=keys,all.x=TRUE,suffixes=c("","_oracle"),sort=FALSE)
  for(metric in intersect(c("field_rmse","field_correlation","delta_rmse","delta_correlation","rmse"),names(scores)))
    merged[[paste0(metric,"_minus_oracle")]] <- merged[[metric]]-merged[[paste0(metric,"_oracle")]]
  merged
}
subset_truth_view <- function(view,species) {
  idx <- match(species,view$species_id); stopifnot(!anyNA(idx)); view$species_id <- species
  for(nm in c("field","states","training_states")) view[[nm]] <- view[[nm]][,idx,drop=FALSE]
  view$prevalence <- view$prevalence[idx]
  p <- view$probabilities; p$species_id <- species
  for(nm in c("p_full","p_space_zero","p_env_zero","p_both_zero","delta_p","field_mean","field_median")) p[[nm]] <- p[[nm]][,idx,drop=FALSE]
  view$probabilities <- p; view
}
load_selected_predictions <- function(run_root,job) {
  directory <- file.path(run_root,"jobs",job$key); path <- file.path(directory,"selection.rds")
  if(!file.exists(path)) return(list(status="not_selected",training=NULL,lattice=NULL))
  selection <- readRDS(path)
  if(is.na(selection$selected_attempt)) return(list(status=selection$status,training=NULL,lattice=NULL,selection=selection))
  selected <- file.path(directory,"attempts",sprintf("%02d",selection$selected_attempt))
  if(md5(file.path(selected,"complete.rds"))!=selection$attempt_hashes[as.character(selection$selected_attempt)]) stop("Selected attempt changed")
  resume_state(selected,selection$identity)
  training <- readRDS(file.path(selected,"training-prediction.rds"))
  pp <- file.path(directory,"prediction/complete.rds")
  lattice <- NULL
  if(file.exists(pp) && resume_state(dirname(pp),selection$identity)$status=="complete") lattice <- readRDS(file.path(dirname(pp),"prediction.rds"))
  list(status=selection$status,training=training,lattice=lattice,selection=selection)
}
score_run <- function(run_root,phase="full") {
  verify_snapshot(run_root); cfg <- readRDS(file.path(run_root,"config.rds"))
  if(!isTRUE(cfg$smoke_only) && !file.exists(file.path(run_root,"runtime-checkpoint.rds"))) stop("Timing checkpoint must precede truth scoring")
  jobs <- make_jobs(run_root)
  if(phase=="pilot") jobs <- subset(jobs,community=="rep01" & arrangement %in% c("spread","clustered"))
  states <- lapply(seq_len(nrow(jobs)),function(i) load_selected_predictions(run_root,as.list(jobs[i,])))
  if(any(vapply(states,function(s) s$status=="not_selected",logical(1)))) stop("Selection must be immutable for every requested cell before scoring")
  manifest <- readRDS(file.path(run_root,"import-manifest.rds")); rows <- strata <- accounting <- list()
  for(i in seq_len(nrow(jobs))) {
    job <- jobs[i,]; state <- states[[i]]
    item <- manifest[manifest$dataset_id==job$dataset_id & manifest$kind=="truth",]
    stopifnot(nrow(item)==1L); truth <- read_checked(file.path(run_root,item$relative_path),item$md5)
    prevalence_map <- setNames(truth$landscape$prevalence,truth$landscape$species)
    if(!is.null(state$training)) {
      selected <- file.path(run_root,"jobs",job$key,"attempts",sprintf("%02d",state$selection$selected_attempt),"fit.rds")
      write_once(native_hyperparameters(readRDS(selected)$fit,truth$landscape$range),file.path(run_root,"scores",job$key,"hyperparameters.rds"))
    }
    accounting[[i]] <- cbind(job,status=state$status,training_supported=!is.null(state$training),lattice_supported=!is.null(state$lattice))
    for(scope in c("training","lattice")) {
      pred <- state[[scope]]
      if(is.null(pred)) next
      tv <- truth_view(truth,scope); if(job$kind=="oracle") tv <- subset_truth_view(tv,job$species)
      env_index <- which(jobs$package==job$package & jobs$dataset_id==job$dataset_id & jobs$arm=="env")
      baseline <- if(length(env_index)==1L) states[[env_index]][[scope]] else NULL
      score <- score_predictions(pred,tv,baseline)
      meta <- job[rep(1,nrow(score$species)),c("key","package","arm","community","arrangement")]; rownames(meta) <- NULL
      rows[[length(rows)+1L]] <- cbind(meta,scope=scope,selected_status=state$status,score$species)
      meta <- job[rep(1,nrow(score$strata)),c("key","package","arm","community","arrangement")]; rownames(meta) <- NULL
      strata[[length(strata)+1L]] <- cbind(meta,scope=scope,score$strata)
      write_once(score,file.path(run_root,"scores",job$key,paste0(scope,".rds")))
    }
  }
  accounting <- do.call(rbind,accounting); species <- if(length(rows)) do.call(rbind,rows) else data.frame(); bins <- if(length(strata)) do.call(rbind,strata) else data.frame()
  species <- complete_species_grid(species,accounting,prevalence_map)
  # Explicit expected denominators accompany every metric row, including failed cells.
  groups <- list()
  for(package in names(cfg$arms)) for(arm in cfg$arms[[package]]) for(arr in cfg$arrangements) for(scope in c("training","lattice")) for(prevalence in c(.05,.25,.75)) {
    sub <- species[species$package==package & species$arm==arm & species$arrangement==arr & species$scope==scope & species$prevalence==prevalence,,drop=FALSE]
    cells <- accounting[accounting$package==package & accounting$arm==arm & accounting$arrangement==arr,,drop=FALSE]
    for(metric in c("field_rmse","field_correlation","field_attenuation","delta_rmse","delta_correlation","delta_rms","rmse","space_zero_benefit","env_fit_benefit","probability_inclusion","probability_width")) {
      if(!metric %in% names(sub)) sub[[metric]] <- rep(NA_real_,nrow(sub))
      if(!"community" %in% names(sub)) sub$community <- character()
      summary <- equal_community_summary(sub,metric,if(phase=="pilot") 1L else 3L)
      summary$rows_expected <- (if(prevalence==.05) 2L else 3L)*(if(phase=="pilot") 1L else 3L)
      groups[[length(groups)+1L]] <- cbind(data.frame(package,arm,arrangement=arr,scope,prevalence,
        cells_available=sum(cells[[paste0(scope,"_supported")]]),cells_expected=if(phase=="pilot") 1L else 3L),summary)
    }
  }
  result <- list(schema_version=1L,accounting=accounting,species=species,strata=bins,groups=do.call(rbind,groups),
    oracle=if(nrow(species)) compare_oracle(subset(species,package!="oracle"),subset(species,package=="oracle")) else data.frame(),
    input_mode=unique(readRDS(file.path(run_root,"comparison-manifest.rds"))$input_mode),
    snapshot_hash=md5(file.path(run_root,"snapshot-record.rds")),phase=phase)
  write_once(result,file.path(run_root,paste0("scores-",phase,".rds"))); result
}
complete_species_grid <- function(species,accounting,prevalence_map) {
  rows <- list()
  for(i in seq_len(nrow(accounting))) for(scope in c("training","lattice")) {
    job <- accounting[i,]; ids <- if(job$kind=="oracle") job$species else names(prevalence_map)
    rows[[length(rows)+1L]] <- data.frame(key=job$key,scope=scope,species_id=ids,package=job$package,arm=job$arm,
      community=job$community,arrangement=job$arrangement,prevalence=unname(prevalence_map[ids]),selected_status=job$status,
      prediction_status=if(job[[paste0(scope,"_supported")]]) "available" else "unavailable")
  }
  grid <- do.call(rbind,rows); grid$.ordering <- seq_len(nrow(grid))
  keys <- c("key","scope","species_id")
  if(nrow(species)) {
    if(anyDuplicated(species[keys])) stop("Duplicate species score identity")
    retained <- c(keys,setdiff(names(species),names(grid)))
    grid <- merge(grid,species[retained],by=keys,all.x=TRUE,sort=FALSE)
  }
  grid <- grid[order(grid$.ordering),]; grid$.ordering <- NULL; rownames(grid) <- NULL
  grid
}
