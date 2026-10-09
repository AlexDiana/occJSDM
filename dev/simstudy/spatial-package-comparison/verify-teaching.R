# Archive-free arithmetic: base R only; never loads a fitting package or raw fit.
verify_teaching_bundle <- function(bundle_path,repo_root=NULL) {
  if(is.character(bundle_path)) {
    checksum <- paste0(bundle_path,".md5")
    if(!file.exists(checksum) || !identical(unname(tools::md5sum(bundle_path)),readLines(checksum,warn=FALSE))) stop("Teaching file checksum mismatch")
    b <- readRDS(bundle_path)
  } else b <- bundle_path
  cfg <- b$config; checks <- list()
  check <- function(ok,label) { if(!isTRUE(ok)) stop("Teaching verification failed: ",label); checks[[length(checks)+1L]] <<- data.frame(check=label,passed=TRUE) }
  check(identical(names(b$landscapes),sprintf("rep%02d",1:3)),"all three communities")
  datasets <- as.vector(outer(cfg$communities,cfg$arrangements,paste,sep="-"))
  check(setequal(names(b$training),datasets) && setequal(names(b$truth_index),datasets) && setequal(names(b$prediction_data),datasets),"all 12 training and prediction inputs")
  a <- b$scores$accounting
  check(sum(a$kind=="comparison")==192L && !anyDuplicated(a$key),"192 comparison cells accounted for")
  for(ds in datasets) {
    tr <- b$training[[ds]]; idx <- b$truth_index[[ds]]; land <- b$landscapes[[idx$community]]
    check(is.matrix(tr$y) && identical(dim(tr$y),c(length(tr$site_id),length(cfg$species))) && all(tr$y %in% 0:1),paste(ds,"binary observations"))
    check(identical(tr$species_id,cfg$species) && identical(colnames(tr$y),cfg$species),paste(ds,"training species order"))
    check(isTRUE(all.equal(unname(tr$y),unname(land$z[idx$training_index,,drop=FALSE]),tolerance=0)),paste(ds,"binary states agree with landscape"))
    check(isTRUE(all.equal(unname(tr$xy),unname(land$points[idx$training_index,,drop=FALSE]),tolerance=0)),paste(ds,"training coordinates"))
    check(all(is.finite(land$psi)) && all(land$psi>=0 & land$psi<=1) && all(land$z %in% 0:1),paste(ds,"complete landscape probabilities and binary states"))
  }
  near <- function(x,y,label) check(isTRUE(all.equal(unname(x),unname(y),tolerance=1e-10,check.attributes=FALSE)),label)
  for(i in seq_len(nrow(a))) {
    job <- a[i,]; state <- b$predictions[[job$key]]
    idx <- b$truth_index[[job$dataset_id]]; land <- b$landscapes[[idx$community]]
    species <- if(job$kind=="oracle") job$species else cfg$species
    sj <- match(species,land$species)
    for(scope in c("training","lattice")) {
      p <- state[[scope]]; declared <- job[[paste0(scope,"_supported")]]
      check(identical(!is.null(p),declared),paste(job$key,scope,"prediction support"))
      if(is.null(p)) next
      site <- if(scope=="training") idx$training_index else idx$lattice_index
      ids <- if(scope=="training") idx$training_site_id else idx$prediction_site_id
      check(identical(p$species_id,species),paste(job$key,scope,"species order"))
      check(identical(p$site_id,ids),paste(job$key,scope,"site order"))
      cols <- c("p_full","p_space_zero","p_env_zero","p_both_zero")
      check(all(vapply(p[cols],function(x)is.matrix(x)&&identical(dim(x),c(length(site),length(species)))&&all(is.finite(x))&&all(x>=0 & x<=1),logical(1))),paste(job$key,scope,"probabilities"))
      near(p$delta_p,p$p_full-p$p_space_zero,paste(job$key,scope,"probability contrast"))
      check(identical(is.null(p$intervals),p$uncertainty$status=="unavailable"),paste(job$key,scope,"interval support"))
      truth <- land$psi[site,sj,drop=FALSE]; states <- land$z[site,sj,drop=FALSE]
      true_field <- sweep(land$field[site,sj,drop=FALSE],2,colMeans(land$field[idx$training_index,sj,drop=FALSE]),"-")
      raw_environment <- land$environment[site] %o% land$B[sj]
      true_space_zero <- plogis(sweep(raw_environment,2,land$B0[sj]+colMeans(land$field[idx$training_index,sj,drop=FALSE]),"+"))
      delta <- truth-true_space_zero
      s <- b$scores$species[b$scores$species$key==job$key & b$scores$species$scope==scope,,drop=FALSE]
      check(nrow(s)==length(species)&&!anyDuplicated(s$species_id),paste(job$key,scope,"species score accounting"))
      s <- s[match(species,s$species_id),,drop=FALSE]
      q <- pmin(pmax(p$p_full,1e-12),1-1e-12)
      metrics <- list(bias=colMeans(p$p_full-truth),mae=colMeans(abs(p$p_full-truth)),rmse=sqrt(colMeans((p$p_full-truth)^2)),
        brier=colMeans((p$p_full-states)^2),log_score=colMeans(states*log(q)+(1-states)*log1p(-q)),clamps=colSums(q!=p$p_full),
        delta_rmse=sqrt(colMeans((p$delta_p-delta)^2)),delta_rms=sqrt(colMeans(p$delta_p^2)),truth_delta_rms=sqrt(colMeans(delta^2)),
        space_zero_rmse=sqrt(colMeans((p$p_space_zero-truth)^2)))
      metrics$space_zero_benefit <- metrics$space_zero_rmse-metrics$rmse
      for(nm in names(metrics)) near(s[[nm]],metrics[[nm]],paste(job$key,scope,nm))
      if(p$link=="logit") {
        near(s$field_rmse,sqrt(colMeans((p$field_mean-true_field)^2)),paste(job$key,scope,"field_rmse"))
        near(s$field_attenuation,colSums(p$field_mean*true_field)/colSums(true_field^2),paste(job$key,scope,"field attenuation"))
      } else check(all(is.na(s$field_rmse)),paste(job$key,scope,"unavailable native-link field recovery"))
      if(!is.null(p$intervals)) {
        v <- p$intervals
        near(s$probability_width,colMeans(v$probability_upper-v$probability_lower),paste(job$key,scope,"probability interval width"))
        near(s$probability_inclusion,colMeans(truth>=v$probability_lower & truth<=v$probability_upper),paste(job$key,scope,"probability interval inclusion"))
      }
    }
  }
  if(!is.null(b$scores$groups)) for(i in seq_len(nrow(b$scores$groups))) {
    g <- b$scores$groups[i,]; s <- b$scores$species
    keep <- s$package==g$package & s$arm==g$arm & s$arrangement==g$arrangement & s$scope==g$scope & s$prevalence==g$prevalence
    s <- s[keep,,drop=FALSE]
    values <- split(s[[g$metric]],s$community)
    community <- vapply(values,function(x)if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_,numeric(1))
    available <- community[is.finite(community)]
    near(g$mean,if(length(available)) mean(available) else NA_real_,paste("group mean",i))
    near(g$communities_available,length(available),paste("group community denominator",i))
    near(g$rows_available,sum(is.finite(s[[g$metric]])),paste("group available rows",i))
  }
  do.call(rbind,checks)
}
