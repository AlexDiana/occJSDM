attribute_sites <- function(prediction_record,outcomes_or_probabilities,expected=FALSE) {
  p <- prediction_record; y <- outcomes_or_probabilities
  stopifnot(identical(dim(y),dim(p$p_full)),all(is.finite(y)),all(y>=0 & y<=1))
  if(!expected && any(!y %in% 0:1)) stop("Observed attribution requires binary states")
  matrices <- p[c("p_both_zero","p_space_zero","p_env_zero","p_full")]
  if(any(vapply(matrices,function(x) any(!is.finite(x)|x<0|x>1),logical(1)))) stop("Invalid probability")
  clamps <- Reduce(`+`,lapply(matrices,function(x) (x<1e-12)|(x>1-1e-12)))
  ll <- lapply(matrices,function(x) { x <- pmin(pmax(x,1e-12),1-1e-12); y*log(x)+(1-y)*log1p(-x) })
  names(ll) <- c("L0","LE","LS","LES")
  E <- .5*((ll$LE-ll$L0)+(ll$LES-ll$LS)); S <- .5*((ll$LS-ll$L0)+(ll$LES-ll$LE))
  total <- ll$LES-ll$L0
  if(max(abs(E+S-total))>1e-10) stop("Attribution identity failed")
  data.frame(site_id=rep(p$site_id,length(p$species_id)),species_id=rep(p$species_id,each=length(p$site_id)),
    expected=expected,L0=as.vector(ll$L0),LE=as.vector(ll$LE),LS=as.vector(ll$LS),LES=as.vector(ll$LES),
    environment=as.vector(E),space=as.vector(S),total_gain=as.vector(total),clamped=as.integer(clamps))
}
