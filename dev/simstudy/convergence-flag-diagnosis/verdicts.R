# Verdicts of the Task 3 diagnostic fits (README.md, "Decision rules" and "How
# the rules are measured", as amended by AMENDMENT-1.md, rulings R13 to R18).
# Pure functions of tables already computed by anatomy.R and modes.R; they
# read no fit. Source modes.R first (region_pattern(), region_note()).
#
# Interfaces
#
#   extended_verdict(regions, anatomy, quantities = SPECIES6_QUANTITIES)
#     `regions`: anchored_regions() of the extended run's primary (anchored)
#     assignment; `anatomy`: chain_anatomy() rows for species 6 (loadings and
#     other quantities are ignored). One row: verdict ('Slow mixing',
#     'Separated modes', 'one mode only: near-truth', 'one mode only: mirror'
#     or 'Mixed'), the region pattern and R18 note, chain counts per region,
#     the largest all-chain Rhat and per-chain split-Rhat over the eleven
#     species-6 quantities (with where they occur), and whether each is within
#     its limit.
#       Slow mixing: every chain visits both regions and rhat_all_chains is at
#         most 1.01 for all eleven quantities.
#       Separated modes: every chain stays in one region, each region holds a
#         chain, and split_rhat_within_chain is at most 1.05 for every chain
#         and all eleven quantities.
#       one mode only (R15): all chains stay in the same region, named.
#       Mixed: anything else, including any chain in an unknown region (R18).
#
#   explains_verdict(separation, extended = '', quantities = SPECIES6_QUANTITIES)
#     `separation`: chain_separation()$separation rows for species 6 of a
#     variant's chains. A quantity agrees if its label is 'agrees' and its
#     all-chain Rhat is at most 1.05 (R14); a quantity fixed by design is left
#     out. The variant explains the split if every other quantity agrees.
#     Returns list(quantities = one row per quantity, verdict = one row),
#     the verdict flagged uninformative when `extended` (the extended-run
#     verdict) is one mode only.
#
#   occupancy_change(extended, variant)
#     Named lists (by species name) of iterations x chains matrices of
#     mean_psi_original_sites. Per species of `variant`: the pooled posterior
#     means, change (variant minus extended), each run's posterior::mcse_mean
#     and their combination in quadrature, with isolation_checks() applied.
#   isolation_checks(changes)
#     Adds below_change_limit (|change| < 0.01), within_2_mcse (|change| at
#     most 2 combined MCSEs) and pass (either) (R13).
#   isolates_verdict(changes, extended = '')
#     One row: 'isolates' if every species passes, else 'does not isolate',
#     with counts, failing species and the uninformative flag.

if(!exists('region_pattern',mode='function') || !exists('region_note',mode='function'))
  stop('Source modes.R before verdicts.R')

# The eleven non-loading quantities of chain_anatomy() (README.md, ruling R6).
SPECIES6_QUANTITIES <- c('B0','theta0','beta_theta_intercept','beta_theta_slope',
  'p_primer1','p_primer2','q_primer1','q_primer2','B_slope1','B_slope2','mean_psi_original_sites')
SLOW_MIXING_RHAT <- 1.01
SEPARATED_SPLIT_RHAT <- 1.05
AGREES_RHAT <- 1.05
ISOLATES_CHANGE <- .01
ISOLATES_MCSE <- 2

one_mode_only <- function(verdict) startsWith(as.character(verdict),'one mode only')
join_names <- function(x) paste(x,collapse=';')
# TRUE only when every value is present and within the limit.
all_within <- function(x,limit) length(x)>0L && !anyNA(x) && all(x<=limit)

need_quantities <- function(have,quantities) {
  absent <- setdiff(quantities,unique(have))
  if(length(absent)) stop('Table lacks species-6 quantities: ',paste(absent,collapse=', '))
}

extended_verdict <- function(regions,anatomy,quantities=SPECIES6_QUANTITIES) {
  need_quantities(anatomy$quantity,quantities)
  a <- anatomy[anatomy$quantity %in% quantities,,drop=FALSE]
  chains <- sort(unique(regions$chain))
  for(q in quantities) {
    ch <- a$chain[a$quantity==q]
    if(anyDuplicated(ch)) stop('Anatomy is not one row per chain for ',q)
    if(!identical(sort(as.integer(ch)),as.integer(chains))) stop('Anatomy chains for ',q,' differ from the region chains')
  }
  pattern <- region_pattern(regions)
  rhat <- tapply(a$rhat_all_chains,a$quantity,function(x) unique(x)[1])[quantities]
  rhat_ok <- all_within(rhat,SLOW_MIXING_RHAT)
  split_ok <- all_within(a$split_rhat_within_chain,SEPARATED_SPLIT_RHAT)
  top <- a[which.max(a$split_rhat_within_chain),,drop=FALSE]
  verdict <- if(identical(pattern,'every chain visits both') && rhat_ok) 'Slow mixing' else
    if(identical(pattern,'each chain in one region') && split_ok) 'Separated modes' else
    if(one_mode_only(pattern)) pattern else 'Mixed'
  count <- function(r) sum(regions$region==r)
  data.frame(verdict=verdict,pattern=pattern,note=region_note(regions),n_chains=length(chains),
    n_near_truth=count('near-truth'),n_mirror=count('mirror'),n_both=count('both'),n_unknown=count('unknown'),
    max_rhat_all_chains=max(rhat,na.rm=TRUE),max_rhat_quantity=names(rhat)[which.max(rhat)],
    rhat_ok=rhat_ok,max_split_rhat_within_chain=top$split_rhat_within_chain,
    max_split_chain=as.integer(top$chain),max_split_quantity=top$quantity,split_rhat_ok=split_ok,
    stringsAsFactors=FALSE)
}

explains_verdict <- function(separation,extended='',quantities=SPECIES6_QUANTITIES) {
  need_quantities(separation$quantity,quantities)
  s <- separation[match(quantities,separation$quantity),,drop=FALSE]
  fixed <- s$fixed %in% TRUE
  agrees <- s$label=='agrees' & !is.na(s$rhat) & s$rhat<=AGREES_RHAT
  detail <- data.frame(quantity=s$quantity,label=s$label,separation=s$separation,
    max_split_rhat_within_chain=s$max_split_rhat_within_chain,rhat=s$rhat,fixed=fixed,
    considered=!fixed,agrees=agrees & !fixed,stringsAsFactors=FALSE)
  considered <- detail[detail$considered,,drop=FALSE]
  explains <- nrow(considered)>0L && all(considered$agrees)
  verdict <- data.frame(verdict=if(explains) 'explains' else 'does not explain',explains=explains,
    n_considered=nrow(considered),n_agree=sum(considered$agrees),
    excluded_fixed=join_names(detail$quantity[detail$fixed]),
    not_agreeing=join_names(considered$quantity[!considered$agrees]),
    max_rhat=if(any(!is.na(considered$rhat))) max(considered$rhat,na.rm=TRUE) else NA_real_,
    max_separation=if(any(is.finite(considered$separation))) max(considered$separation[is.finite(considered$separation)]) else NA_real_,
    uninformative=one_mode_only(extended),extended_verdict=as.character(extended),stringsAsFactors=FALSE)
  list(quantities=detail,verdict=verdict)
}

isolation_checks <- function(changes) {
  changes$below_change_limit <- abs(changes$change)<ISOLATES_CHANGE
  changes$within_2_mcse <- abs(changes$change)<=ISOLATES_MCSE*changes$combined_mcse
  changes$pass <- changes$below_change_limit | changes$within_2_mcse
  changes
}

occupancy_change <- function(extended,variant) {
  species <- names(variant)
  if(is.null(species) || any(!nzchar(species))) stop('Variant draws must be named by species')
  absent <- setdiff(species,names(extended))
  if(length(absent)) stop('Species not in the extended run: ',paste(absent,collapse=', '))
  row <- function(s) {
    e <- as.matrix(extended[[s]]);v <- as.matrix(variant[[s]])
    se_e <- posterior::mcse_mean(e);se_v <- posterior::mcse_mean(v)
    data.frame(species_name=s,extended_chains=ncol(e),variant_chains=ncol(v),
      extended_mean=mean(e),variant_mean=mean(v),change=mean(v)-mean(e),
      extended_mcse=se_e,variant_mcse=se_v,combined_mcse=sqrt(se_e^2+se_v^2),stringsAsFactors=FALSE)
  }
  out <- do.call(rbind,lapply(species,row));rownames(out) <- NULL
  isolation_checks(out)
}

isolates_verdict <- function(changes,extended='') {
  if(!nrow(changes)) stop('No species to compare')
  pass <- changes$pass %in% TRUE
  isolates <- all(pass)
  data.frame(verdict=if(isolates) 'isolates' else 'does not isolate',isolates=isolates,
    n_species=nrow(changes),n_pass=sum(pass),n_by_change=sum(changes$below_change_limit %in% TRUE),
    n_by_mcse_only=sum(pass & !(changes$below_change_limit %in% TRUE)),
    failing=join_names(changes$species_name[!pass]),max_abs_change=max(abs(changes$change)),
    uninformative=one_mode_only(extended),extended_verdict=as.character(extended),stringsAsFactors=FALSE)
}
