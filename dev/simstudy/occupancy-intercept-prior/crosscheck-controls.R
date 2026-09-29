#!/usr/bin/env Rscript
# Cross-check of the scorer's control outcomes against the archived control
# results (R10: the archives are a cross-check only, never the source).
#
#   Rscript crosscheck-controls.R --repo=REPO --study=STUDY --phase=A1|A2
#     [--archives=DIR] [--summary=DIR] [--out=DIR]
#
# For every selected control fit, compares summarise.R's tables (SUMMARY
# defaults to STUDY/summary) and analysis.R's cell records with the archived
# quantities that are defined the same way:
#   A1  the pr11-current-20260927 result file of the selected fit and the
#       published current-main-recheck/results/community-scores.csv ("current"
#       version): per-cell posterior-mean probability and truth; signed error,
#       mean absolute error and cell count in each band of the original 100
#       sites and, for 300-site fits, of all sites; the rare group, recomputed
#       from the archived per-cell estimates and truth. The A1 archives record
#       no interval and no B0 summary, so coverage and B0 are not cross-checked.
#   A2  the spatial-amplitude robust-v1 result file of the selected fit and the
#       published robust-v1/summary-binary-final/groups.csv: per-cell estimate,
#       interval and containment; signed error, MAE and coverage in each band
#       and prevalence group; the 1% plus 5% group from the archived cells; B0
#       bias, absolute bias and coverage.
# Writes OUT/control-validation-<phase>.csv (one row per control) and
# OUT/control-validation-summary-<phase>.csv (largest difference per
# quantity); OUT defaults to the study directory's results/. Exits 1 if any
# difference exceeds 1e-10 (points, probability or logit) or any containment,
# cell count or group set differs.

TOL <- 1e-10

mx <- function(x) if(!length(x)) NA_real_ else max(x)

crosscheck_a1 <- function(key,fit_label,rec,be,archives,pub) {
  r <- readRDS(file.path(archives,sub('-fit[.]rds$','-result.rds',fit_label)))
  b <- be[be$key==key,,drop=FALSE];p <- pub[pub$key==key & pub$version=='current',,drop=FALSE]
  n <- nrow(r$truth);cells_ok <- TRUE;de <- dm <- pe <- pm <- numeric()
  for(k in seq_len(nrow(r$groups))) {
    g <- r$groups[k,];scope <- if(g$scope=='original100') 'primary' else 'allsites'
    if(scope=='allsites' && n<=100L) next
    m <- b[b$scope==scope & b$group==g$group,,drop=FALSE]
    if(nrow(m)!=1L || m$cells!=g$cells) {cells_ok <- FALSE;next}
    de <- c(de,abs(m$signed_error-100*g$bias));dm <- c(dm,abs(m$mean_abs_cell_error-100*g$mae))
    metric <- if(scope=='primary') 'occupancy_original_sites' else 'occupancy_all_sites'
    q <- p[p$metric==metric & p$group==if(g$group=='middle') 'medium' else g$group,,drop=FALSE]
    if(nrow(q)!=1L || q$n_elements!=m$cells) {cells_ok <- FALSE;next}
    pe <- c(pe,abs(m$signed_error-100*q$bias));pm <- c(pm,abs(m$mean_abs_cell_error-100*q$mae))
  }
  rare <- which(colMeans(r$truth)<.2);rm <- b[b$group=='rare_below_20pct',,drop=FALSE];rare_d <- NA_real_
  if(length(rare)) {
    idx <- unlist(lapply(rare,function(s) (s-1L)*n+seq_len(n)));e <- as.vector(r$estimate)[idx]-as.vector(r$truth)[idx]
    if(nrow(rm)!=1L || rm$cells!=length(idx)) cells_ok <- FALSE else
      rare_d <- max(abs(rm$signed_error-100*mean(e)),abs(rm$mean_abs_cell_error-100*mean(abs(e))))
  } else if(nrow(rm)) cells_ok <- FALSE
  data.frame(phase='A1',key=key,cells=nrow(rec$cells),archive=basename(sub('-fit[.]rds$','-result.rds',fit_label)),
    cells_and_groups_equal=cells_ok,groups_compared=length(de),
    estimate=max(abs(rec$cells$estimate-as.vector(r$estimate))),truth=max(abs(rec$cells$truth-as.vector(r$truth))),
    signed_error=mx(de),mean_abs_cell_error=mx(dm),published_signed_error=mx(pe),published_mean_abs_cell_error=mx(pm),
    rare_group=rare_d,rare_species=length(rare),stringsAsFactors=FALSE)
}

crosscheck_a2 <- function(key,rec,be,cv,b0,archives,fits,pub) {
  row <- fits[fits$prior=='inverse_gamma' & fits$key==key,,drop=FALSE]
  if(nrow(row)!=1L || !identical(row$fit_md5,rec$fit_md5) || !identical(unname(tools::md5sum(row$result_file)),row$result_md5))
    stop('No archived robust-v1 result of the selected control fit ',key)
  r <- readRDS(row$result_file)
  el <- r$elements[r$elements$metric=='occupancy',,drop=FALSE];stopifnot(identical(el$element,seq_len(nrow(rec$cells))))
  b <- be[be$key==key,,drop=FALSE];c <- cv[cv$key==key,,drop=FALSE];p <- pub[pub$key==key & pub$prior=='inverse_gamma',,drop=FALSE]
  cells_ok <- TRUE;de <- dm <- dc <- pe <- pm <- pc <- numeric()
  og <- r$groups[r$groups$metric=='occupancy',,drop=FALSE]
  for(k in seq_len(nrow(og))) {
    g <- og[k,];name <- if(g$group=='medium') 'middle' else g$group
    m <- b[b$scope=='primary' & b$group==name,,drop=FALSE];mc <- c[c$scope=='primary' & c$group==name,,drop=FALSE]
    q <- p[p$metric=='occupancy' & p$group==g$group,,drop=FALSE]
    if(nrow(m)!=1L || nrow(mc)!=1L || nrow(q)!=1L || m$cells!=g$n || q$n!=g$n) {cells_ok <- FALSE;next}
    de <- c(de,abs(m$signed_error-100*g$bias));dm <- c(dm,abs(m$mean_abs_cell_error-100*g$mae));dc <- c(dc,abs(mc$coverage-g$coverage))
    pe <- c(pe,abs(m$signed_error-100*q$bias));pm <- c(pm,abs(m$mean_abs_cell_error-100*q$mae));pc <- c(pc,abs(mc$coverage-q$coverage))
  }
  if(!setequal(b$group[b$scope=='primary'],c(ifelse(og$group=='medium','middle',og$group),'rare_1_5pct'))) cells_ok <- FALSE
  sp <- which(abs(r$species$target-.01)<1e-10 | abs(r$species$target-.05)<1e-10);n <- nrow(rec$cells)/nrow(r$species)
  idx <- unlist(lapply(sp,function(s) (s-1L)*n+seq_len(n)))
  rm <- b[b$group=='rare_1_5pct',,drop=FALSE];rc <- c[c$group=='rare_1_5pct',,drop=FALSE]
  rare_d <- if(nrow(rm)==1L && rm$cells==length(idx)) max(abs(rm$signed_error-100*mean(el$bias[idx])),
    abs(rm$mean_abs_cell_error-100*mean(abs(el$bias[idx]))),abs(rc$coverage-mean(el$covered[idx]))) else {cells_ok <- FALSE;NA_real_}
  ie <- r$elements[r$elements$metric=='intercept',,drop=FALSE];ig <- r$groups[r$groups$metric=='intercept' & r$groups$group=='all',]
  bb <- b0[b0$key==key,,drop=FALSE]
  b0_d <- max(abs(bb$b0_bias-mean(ie$bias)),abs(bb$b0_abs_bias-mean(abs(ie$bias))),abs(bb$b0_coverage-mean(ie$covered)),
    abs(bb$b0_bias-ig$bias),abs(bb$b0_coverage-ig$coverage),abs(rec$b0_species$bias-ie$bias))
  data.frame(phase='A2',key=key,cells=nrow(rec$cells),archive=basename(row$result_file),cells_and_groups_equal=cells_ok,
    groups_compared=length(de),estimate=max(abs(rec$cells$estimate-el$estimate)),truth=max(abs(rec$cells$truth-el$truth)),
    lower=max(abs(rec$cells$lower-el$lower)),upper=max(abs(rec$cells$upper-el$upper)),containment_differences=sum(rec$cells$covered!=el$covered),
    signed_error=mx(de),mean_abs_cell_error=mx(dm),coverage=mx(dc),published_signed_error=mx(pe),published_mean_abs_cell_error=mx(pm),
    published_coverage=mx(pc),rare_group=rare_d,b0=b0_d,stringsAsFactors=FALSE)
}

crosscheck_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','archives','phase','summary','out'),required=c('repo','study','phase'))
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  summary <- o$summary %||% file.path(study,'summary');out <- o$out %||% file.path(repo,ANALYSIS_REL,'results')
  phase <- o$phase;if(!phase %in% c('A1','A2')) stop('--phase must be A1 or A2')
  read <- function(name) {x <- utils::read.csv(file.path(summary,phase,name),stringsAsFactors=FALSE,
    colClasses=c(key='character',community='character',fit='character',fit_md5='character'));x[x$sd==1 & x$kind=='selected',,drop=FALSE]}
  be <- read('band-error.csv');cv <- read('coverage.csv');b0 <- read('b0.csv')
  ctl <- control_items(phase,study,archives,repo)
  if(!setequal(unique(be$key),ctl$key)) stop('The summary does not hold every control of phase ',phase)
  rec <- function(i) {r <- readRDS(score_path(summary,phase,1,ctl$key[i],ctl$schedule[i]))
    if(!identical(r$fit_md5,ctl$expected_md5[i])) stop('Record of another fit: ',ctl$key[i]);r}
  rows <- if(phase=='A1') {
    pub <- utils::read.csv(file.path(repo,'dev/simstudy/current-main-recheck/results/community-scores.csv'),stringsAsFactors=FALSE)
    lapply(seq_len(nrow(ctl)),function(i) crosscheck_a1(ctl$key[i],ctl$fit_label[i],rec(i),be,archives,pub))
  } else {
    final <- file.path(archive_dir(archives,'amplitude'),'robust-v1/summary-binary-final')
    fits <- utils::read.csv(file.path(final,'fits.csv'),stringsAsFactors=FALSE)
    pub <- utils::read.csv(file.path(final,'groups.csv'),stringsAsFactors=FALSE)
    lapply(seq_len(nrow(ctl)),function(i) crosscheck_a2(ctl$key[i],rec(i),be,cv,b0,archives,fits,pub))
  }
  x <- do.call(rbind,rows)
  q <- setdiff(names(x)[vapply(x,is.numeric,logical(1))],c('cells','groups_compared','rare_species','containment_differences'))
  s <- data.frame(phase=phase,quantity=q,controls=nrow(x),compared=vapply(q,function(k) sum(!is.na(x[[k]])),integer(1)),
    max_abs_difference=vapply(q,function(k) mx(x[[k]][!is.na(x[[k]])]),numeric(1)),stringsAsFactors=FALSE)
  s <- rbind(s,data.frame(phase=phase,quantity=c('cells_and_groups_unequal',if(phase=='A2') 'containment_differences'),controls=nrow(x),
    compared=nrow(x),max_abs_difference=c(sum(!x$cells_and_groups_equal),if(phase=='A2') sum(x$containment_differences)),stringsAsFactors=FALSE))
  rownames(s) <- NULL
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(x,file.path(out,paste0('control-validation-',phase,'.csv')),row.names=FALSE)
  utils::write.csv(s,file.path(out,paste0('control-validation-summary-',phase,'.csv')),row.names=FALSE)
  print(s,row.names=FALSE)
  ok <- all(x$cells_and_groups_equal) && all(s$max_abs_difference[s$quantity %in% q]<=TOL,na.rm=TRUE) &&
    (phase!='A2' || sum(x$containment_differences)==0L)
  cat(if(ok) 'All' else 'NOT all',nrow(x),'phase',phase,'controls reproduce the archived results.\n')
  if(ok) 0L else 1L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R')) source(file.path(here,f))
  status <- tryCatch(crosscheck_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
