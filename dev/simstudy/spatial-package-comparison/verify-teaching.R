# Archive-free arithmetic: base R only; never loads a fitting package or raw fit.
verify_teaching_bundle <- function(bundle_path,repo_root=NULL) {
  if(is.character(bundle_path)) {
    checksum <- paste0(bundle_path,".md5")
    if(!file.exists(checksum) || !identical(unname(tools::md5sum(bundle_path)),readLines(checksum,warn=FALSE))) stop("Teaching file checksum mismatch")
    b <- readRDS(bundle_path)
  } else b <- bundle_path
  cfg <- b$config; checks <- new.env(parent=emptyenv())
  check <- function(ok,label) {
    if(!isTRUE(ok)) stop("Teaching verification failed: ",label)
    # Keep one compact receipt per stratum audit, with the actual assertion count.
    label <- sub("stratum.*$","stratum arithmetic and accounting",label)
    checks[[label]] <- if(is.null(checks[[label]])) 1L else checks[[label]]+1L
  }
  check(identical(names(b$landscapes),sprintf("rep%02d",1:3)),"all three communities")
  datasets <- as.vector(outer(cfg$communities,cfg$arrangements,paste,sep="-"))
  check(setequal(names(b$training),datasets) && setequal(names(b$truth_index),datasets) && setequal(names(b$prediction_data),datasets),"all 12 training and prediction inputs")
  a <- b$scores$accounting
  check(sum(a$kind=="comparison")==192L && !anyDuplicated(a$key),"192 comparison cells accounted for")
  if(identical(b$scores$phase,"full")) {
    check(sum(a$kind=="oracle")==96L,"96 oracle cells accounted for")
    s <- b$scores$species
    check(nrow(s)==2L*(192L*length(cfg$species)+96L) && !anyDuplicated(s[c("key","scope","species_id")]),"full species denominator")
    for(i in seq_len(nrow(a))) for(scope in c("training","lattice")) {
      job <- a[i,]; rows <- s[s$key==job$key & s$scope==scope,,drop=FALSE]
      species <- if(job$kind=="oracle") job$species else cfg$species
      check(setequal(rows$species_id,species),paste(job$key,scope,"species denominator"))
      check(all(rows$prediction_status==if(job[[paste0(scope,"_supported")]]) "available" else "unavailable"),paste(job$key,scope,"species availability denominator"))
    }
  }
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
      env_key <- a$key[a$package==job$package & a$arm=="env" & a$dataset_id==job$dataset_id]
      baseline <- if(length(env_key)==1L) b$predictions[[env_key]][[scope]] else NULL
      metrics$env_fit_rmse <- if(is.null(baseline)) rep(NA_real_,length(species)) else sqrt(colMeans((baseline$p_full-truth)^2))
      metrics$env_fit_benefit <- metrics$env_fit_rmse-metrics$rmse
      metrics$training_occurrences <- colSums(land$z[idx$training_index,sj,drop=FALSE])
      metrics$generating_spatial_sd <- rep(land$field_sd,length(species))
      for(nm in names(metrics)) near(s[[nm]],metrics[[nm]],paste(job$key,scope,nm))
      verify_field <- function(estimate,truth,prefix,compatible) {
        zero <- sqrt(colMeans(truth^2)); error <- sqrt(colMeans((estimate-truth)^2))
        constant <- function(x)diff(range(x))<1e-12
        ec <- apply(estimate,2,constant); tc <- apply(truth,2,constant)
        correlation <- vapply(seq_len(ncol(estimate)),function(j)if(!compatible||ec[j]||tc[j]) NA_real_ else cor(estimate[,j],truth[,j]),numeric(1))
        attenuation <- if(compatible) colSums(estimate*truth)/colSums(truth^2) else rep(NA_real_,length(species))
        attenuation[colSums(truth^2)==0] <- NA_real_
        fraction <- if(compatible) 1-error/zero else rep(NA_real_,length(species)); fraction[zero==0] <- NA_real_
        metrics <- list(rmse=if(compatible) error else rep(NA_real_,length(species)),zero_rmse=zero,
          rmse_reduction=if(compatible) zero-error else rep(NA_real_,length(species)),fractional_rmse_reduction=fraction,
          correlation=correlation,attenuation=attenuation,realised_sd=apply(estimate,2,sd),truth_realised_sd=apply(truth,2,sd))
        for(nm in names(metrics)) near(s[[paste0(prefix,nm)]],metrics[[nm]],paste(job$key,scope,prefix,nm))
        status <- if(!compatible) rep("incompatible link scale",length(species)) else ifelse(ec,"constant fitted field",ifelse(tc,"constant true field","available"))
        check(identical(s[[paste0(prefix,"correlation_status")]],status),paste(job$key,scope,prefix,"correlation status"))
      }
      verify_field(p$field_mean,true_field,"field_",p$link=="logit")
      verify_field(p$delta_p,delta,"delta_",TRUE)
      if(p$link=="logit") {
        near(s$field_rmse,sqrt(colMeans((p$field_mean-true_field)^2)),paste(job$key,scope,"field_rmse"))
        near(s$field_attenuation,colSums(p$field_mean*true_field)/colSums(true_field^2),paste(job$key,scope,"field attenuation"))
      } else check(all(is.na(s$field_rmse)),paste(job$key,scope,"unavailable native-link field recovery"))
      if(!is.null(p$intervals)) {
        v <- p$intervals
        near(s$probability_width,colMeans(v$probability_upper-v$probability_lower),paste(job$key,scope,"probability interval width"))
        near(s$probability_inclusion,colMeans(truth>=v$probability_lower & truth<=v$probability_upper),paste(job$key,scope,"probability interval inclusion"))
        near(s$field_width_native,colMeans(v$field_upper-v$field_lower),paste(job$key,scope,"field interval width"))
        near(s$field_inclusion,if(p$link=="logit") colMeans(true_field>=v$field_lower & true_field<=v$field_upper) else rep(NA_real_,length(species)),paste(job$key,scope,"field interval inclusion"))
      }
      if(!is.null(b$site_attribution[[job$key]][[scope]])) {
        saved <- b$site_attribution[[job$key]][[scope]]
        for(type in c("observed","expected")) {
          weights <- if(type=="observed") states else truth
          logs <- lapply(p[c("p_both_zero","p_space_zero","p_env_zero","p_full")],function(q) {
            q <- pmin(pmax(q,1e-12),1-1e-12); weights*log(q)+(1-weights)*log1p(-q)
          })
          names(logs) <- c("L0","LE","LS","LES")
          E <- (logs$LE-logs$L0+logs$LES-logs$LS)/2
          S <- (logs$LS-logs$L0+logs$LES-logs$LE)/2
          table <- saved[[type]]
          check(identical(table$site_id,rep(ids,length(species))) && identical(table$species_id,rep(species,each=length(ids))),paste(job$key,scope,type,"attribution IDs"))
          for(nm in names(logs)) near(table[[nm]],as.vector(logs[[nm]]),paste(job$key,scope,type,"attribution",nm))
          near(table$environment,as.vector(E),paste(job$key,scope,type,"attribution environment"))
          near(table$space,as.vector(S),paste(job$key,scope,type,"attribution space"))
          near(table$total_gain,as.vector(logs$LES-logs$L0),paste(job$key,scope,type,"attribution total"))
          near(table$environment+table$space,table$total_gain,paste(job$key,scope,type,"attribution identity"))
          clamps <- Reduce(`+`,lapply(p[cols],function(q)q<1e-12|q>1-1e-12))
          near(table$clamped,as.vector(clamps),paste(job$key,scope,type,"attribution clamps"))
          check(all(table$expected==(type=="expected")),paste(job$key,scope,type,"attribution weights"))
        }
      }
      if(!is.null(b$scores$strata)) {
        bins <- b$scores$strata[b$scores$strata$key==job$key & b$scores$strata$scope==scope,,drop=FALSE]
        check(nrow(bins)==length(species)*7L*3L && !anyDuplicated(bins[c("species_id","stratification","bin","target")]),paste(job$key,scope,"stratum accounting"))
        xy <- land$points[site,,drop=FALSE]; tx <- land$points[idx$training_index,,drop=FALSE]
        nearest <- vapply(seq_len(nrow(xy)),function(k)min(sqrt(rowSums(sweep(tx,2,xy[k,],"-")^2))),numeric(1))
        distance <- cut(nearest,c(-Inf,.02,.05,.1,Inf),labels=c("<=0.02","(0.02,0.05]","(0.05,0.1]",">0.1"))
        for(k in seq_len(nrow(bins))) {
          row <- bins[k,]; j <- match(row$species_id,species)
          group <- if(row$stratification=="distance") distance else {
            z <- cut(truth[,j],c(-Inf,.2,.8,Inf),labels=c("<0.2","[0.2,0.8]",">0.8"));z[truth[,j]==.2] <- "[0.2,0.8]";z
          }
          keep <- which(group==row$bin)
          values <- if(row$target=="env_fit") { if(!is.null(baseline)) baseline$p_full[,j] else NULL } else p[[row$target]][,j]
          near(row$n,length(keep),paste(job$key,scope,"stratum count",k))
          if(is.null(values)||!length(keep)) {
            check(row$status==if(is.null(values)) "unavailable env fit" else "empty",paste(job$key,scope,"stratum availability",k))
            check(all(is.na(unlist(row[c("bias","mae","rmse","brier","log_score")]))),paste(job$key,scope,"stratum unavailable metrics",k))
          } else {
            estimate <- values[keep]; t <- truth[keep,j]; y <- states[keep,j]; q <- pmin(pmax(estimate,1e-12),1-1e-12)
            metrics <- list(bias=mean(estimate-t),mae=mean(abs(estimate-t)),rmse=sqrt(mean((estimate-t)^2)),
              brier=mean((estimate-y)^2),log_score=mean(y*log(q)+(1-y)*log1p(-q)),clamps=sum(q!=estimate))
            for(nm in names(metrics)) near(row[[nm]],metrics[[nm]],paste(job$key,scope,"stratum",k,nm))
          }
        }
      }
    }
  }
  if(!is.null(b$scores$oracle) && nrow(b$scores$oracle)) {
    o <- b$scores$oracle; s <- b$scores$species
    keys <- c("community","arrangement","scope","species_id")
    identity <- function(x)do.call(paste,c(x[keys],sep="/"))
    reference <- s[s$package=="oracle",,drop=FALSE]
    stopifnot(!anyDuplicated(identity(reference)))
    ri <- match(identity(o),identity(reference)); mi <- match(paste(o$key,o$scope,o$species_id),paste(s$key,s$scope,s$species_id))
    check(!anyNA(ri)&&!anyNA(mi),"oracle comparison identities")
    check(nrow(o)==sum(s$package!="oracle"),"oracle comparison denominator")
    for(metric in c("field_rmse","field_correlation","delta_rmse","delta_correlation","rmse"))
      near(o[[paste0(metric,"_minus_oracle")]],s[[metric]][mi]-reference[[metric]][ri],paste("oracle gap",metric))
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
    near(g$minimum,if(length(available)) min(available) else NA_real_,paste("group minimum",i))
    near(g$maximum,if(length(available)) max(available) else NA_real_,paste("group maximum",i))
    near(g$communities_expected,3L,paste("group expected communities",i))
    near(g$rows_expected,3L*if(g$prevalence==.05) 2L else 3L*3L,paste("group expected species rows",i))
    cells <- a[a$package==g$package & a$arm==g$arm & a$arrangement==g$arrangement,,drop=FALSE]
    near(g$cells_expected,3L,paste("group expected cells",i))
    near(g$cells_available,sum(cells[[paste0(g$scope,"_supported")]]),paste("group available cells",i))
  }
  labels <- ls(checks,all.names=TRUE)
  data.frame(check=labels,passed=TRUE,assertions=vapply(labels,function(label)checks[[label]],integer(1)))
}
