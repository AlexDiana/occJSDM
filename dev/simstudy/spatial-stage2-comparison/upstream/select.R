# Selection reads diagnostics and training predictions only.
select_attempts <- function(attempt_records,config=study_config(),kind=c("point","bayesian")) {
  kind <- match.arg(kind); a <- attempt_records
  ids <- if(length(a)) as.integer(vapply(a,`[[`,numeric(1),"attempt_id")) else integer()
  if(anyDuplicated(ids)) stop("Duplicate attempt identity")
  if(any(vapply(a,function(x) identical(x$status,"resource_interrupted"),logical(1))))
    return(list(status="budget_amendment_required",selected_attempt=NA_integer_,next_attempts=integer(),attempts=ids))
  returned <- vapply(a,function(x) identical(x$status,"complete"),logical(1))
  acceptable <- returned & vapply(a,function(x) isTRUE(x$acceptable),logical(1))
  result <- function(status,selected=NA_integer_,next_ids=integer(),reason="")
    list(schema_version=1L,status=status,selected_attempt=selected,next_attempts=next_ids,attempts=ids,reason=reason)
  if(kind=="bayesian") {
    if(!length(a)) return(result("pending",next_ids=1L))
    if(any(acceptable)) return(result("selected",selected=max(ids[acceptable]),reason="MCMC diagnostics passed"))
    if(!2L %in% ids) return(result("pending",next_ids=2L,reason="Initial diagnostics or fit failed; prespecified extended schedule"))
    if(!any(returned)) return(result("failed",reason="No complete Bayesian fit"))
    return(result("unresolved_diagnostics",selected=max(ids[returned]),reason="Finite MCMC budget exhausted; diagnostics remain failed"))
  }
  missing <- setdiff(1:config$diagnostics$initial_starts,ids)
  if(length(missing)) return(result("pending",next_ids=missing))
  if(!any(acceptable)) return(result("failed",reason="No acceptable solution among attempted starts"))
  eligible <- which(acceptable); obj <- vapply(a[eligible],`[[`,numeric(1),"objective")
  best <- eligible[order(obj,ids[eligible])[1]]
  agree <- vapply(setdiff(eligible,best),function(j)
    abs(a[[j]]$objective-a[[best]]$objective)<=config$diagnostics$objective_agreement &&
    max(abs(a[[j]]$training_probability-a[[best]]$training_probability))<=config$diagnostics$probability_agreement,logical(1))
  if(any(agree)) return(result("selected",ids[best],reason="Best acceptable objective replicated within objective and training-probability tolerances"))
  missing <- setdiff(1:config$diagnostics$max_starts,ids)
  if(length(missing)) return(result("pending",next_ids=missing,reason="Best acceptable solution not replicated"))
  result("unstable",ids[best],reason="Best acceptable objective selected; agreement unresolved after six starts")
}
runtime_gate <- function(remaining_hours,limit_hours) {
  if(length(limit_hours)!=1L || !is.finite(limit_hours) || limit_hours<=0) stop("Runtime limit is unset")
  list(passed=is.finite(remaining_hours)&&remaining_hours<=limit_hours,
       remaining_hours=remaining_hours,limit_hours=limit_hours)
}
