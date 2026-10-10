# Reporting audit only: ordinary compact data, no fitting packages or archive.
# Run after the export's hash/process completion receipt has passed.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L, file.exists(args[1]), !file.exists(args[2]))
b <- readRDS(args[1])
stopifnot(all(b$verification$compact$passed))
checks <- list(); distances <- list(); metadata <- list()
check <- function(ok,label) {
  if(!isTRUE(ok)) stop(label)
  checks[[length(checks)+1L]] <<- data.frame(check=label,passed=TRUE)
}
near <- function(x,y) length(x)==length(y) && all(is.finite(c(x,y))) && all(abs(x-y)<1e-10)
for(id in names(b$training)) {
  training <- b$training[[id]]
  xy <- training$xy; scale <- apply(xy,2,sd)
  h <- b$hyperparameters[[paste0("occJSDM-gp-",id)]]
  ell <- h[h$parameter=="native_length_scale",c("lower","median","upper")]
  check(nrow(ell)==1L,paste(id,"GP length-scale row"))
  for(k in 1:2) {
    row <- h[h$parameter==paste0("half_correlation_distance_axis_",k),]
    check(near(unlist(row[c("lower","median","upper")]),
               sqrt(2*log(2))*scale[k]*unlist(ell)),paste(id,"GP raw-axis conversion",k))
    check(identical(row$units,"raw unit-square coordinate distance"),paste(id,"GP distance units",k))
  }
  unequal <- abs(scale[1]-scale[2])>=1e-12
  check(!unequal || !any(h$range_recovery_supported),paste(id,"anisotropic scaling excludes scalar range recovery"))
  check(grepl("distinct from realised",h$definition[h$parameter=="spatial_prior_sd"]),paste(id,"coefficient prior SD label"))
  hh <- b$hyperparameters[[paste0("Hmsc-spatial2-",id)]]
  check(nrow(hh)==2L && all(hh$lower<=hh$median & hh$median<=hh$upper),paste(id,"Hmsc interval ordering"))
  check(hh$median[1]<=hh$median[2] && all(grepl("ordered within each draw",hh$definition)),paste(id,"Hmsc factor ordering definition"))
  check(all(grepl("zero denotes spatial independence",hh$definition)) && !any(hh$range_recovery_supported),paste(id,"Hmsc range semantics"))
  hg <- b$hyperparameters[[paste0("gllvm-spatial2-",id)]]
  check(nrow(hg)==2L && all(is.finite(hg$median) & hg$median>=0),paste(id,"gllvm finite point ranges"))
  check(all(is.na(hg$lower) & is.na(hg$upper)) && !any(hg$range_recovery_supported),paste(id,"gllvm interval and kernel limits"))
  for(package in c("occJSDM","Hmsc","gllvm")) {
    hp <- b$hyperparameters[[paste0(package,"-",if(package=="occJSDM") "gp" else "spatial2","-",id)]]
    hp$key <- paste0(package,"-",if(package=="occJSDM") "gp" else "spatial2","-",id)
    hp$selected_status <- b$selections[[hp$key[1]]]$status
    metadata[[length(metadata)+1L]] <- hp
  }
  ti <- b$truth_index[[id]]; land <- b$landscapes[[ti$community]]
  lattice <- land$points[ti$lattice_index,,drop=FALSE]
  nearest <- vapply(seq_len(nrow(lattice)),function(i) sqrt(min(rowSums((xy-matrix(lattice[i,],nrow(xy),2,byrow=TRUE))^2))),numeric(1))
  saved <- b$basis_summaries[[id]]$distance
  check(near(unname(quantile(nearest)),unname(saved$nearest_training_quantiles)),paste(id,"nearest-distance quantiles"))
  check(near(max(nearest),saved$maximum_nearest_training_distance),paste(id,"nearest-distance maximum"))
  check(identical(b$basis_summaries[[id]]$fitted_range_available,FALSE),paste(id,"Moran graph is not fitted range"))
  bins <- table(cut(nearest,c(-Inf,.02,.05,.1,Inf),labels=c("<=0.02","(0.02,0.05]","(0.05,0.1]",">0.1")))
  check(sum(bins)==1600L,paste(id,"distance bins retain all lattice sites"))
  distances[[id]] <- data.frame(dataset=id,bin=names(bins),sites=as.integer(bins),sites_expected=1600L,
    nearest_median=median(nearest),nearest_maximum=max(nearest),
    fraction_within_true_half_correlation=mean(nearest<=sqrt(2*log(2))*land$range))
}
check(length(intersect(loadedNamespaces(),c("occJSDM","sjSDM","Hmsc","gllvm")))==0L,"archive-free audit")
result <- list(passed=TRUE,checks=do.call(rbind,checks),native=do.call(rbind,metadata),
  distances=do.call(rbind,distances),input_md5=unname(tools::md5sum(args[1])),
  generated=Sys.time(),interpretation="Axis-specific base-kernel distances are descriptive. Exponential factor ranges, coefficient prior SD and Moran graph distances are not recovery estimates of the isotropic generating squared-exponential range.")
saveRDS(result,args[2],compress="xz")
cat("Native/distance reporting audit passed",nrow(result$checks),"checks\n")
