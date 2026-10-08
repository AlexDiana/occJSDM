# Frozen protocol. Scientific changes need a dated spec amendment first.
study_config <- function() {
  list(schema_version=1L, input_modes=c("reconstructed", "restored"),
    arms=list(occJSDM=c("env","trend","mem","gp"),
              sjSDM=c("env","trend","mem","dnn"),
              Hmsc=c("env","trend","mem","spatial2"),
              gllvm=c("env","trend","mem","spatial2")),
    communities=sprintf("rep%02d",1:3), arrangements=c("spread","pairs","clustered","grid"),
    species=sprintf("species%02d",1:8),
    occ_revision="9af859795974839ee6038d00b6d8422ae77b047f",
    sj_revision="d2ca508853a6e39df493f87c21d9c0136cfe652b",
    generator_md5="80a8f43c3654d33a5551049273670402",
    bundle_md5="e2d6374046b60dadde6063dd6ed17ad1", tolerance=1e-12,
    workers=4L, threads=1L, remaining_wall_hours=48, attempt_hours=24,
    mcmc=list(initial=c(nchain=4L,nburn=3000L,niter=4000L,nthin=1L),
              extended=c(nchain=4L,nburn=8000L,niter=16000L,nthin=1L)),
    diagnostics=list(rhat=1.01,ess=400,gradient=0.01,objective_agreement=0.1,
                     probability_agreement=0.01,initial_starts=3L,max_starts=6L),
    sjsdm=list(hidden=c(30L,30L),activation="selu",dropout=0,bias=TRUE,
      lambda=0.002,alpha=0.2,env_lambda=0,weight_decay=0.0001,
      epochs=c(3000L,1000L),lr=c(0.002,0.0002),sampling=2000L,
      link="probit",dtype="float64",device="cpu",diag=TRUE),
    gllvm=list(maxit=12000L,method="VA",family="binomial",link="logit"),
    limits=c(comparison=768L,oracle=192L,total=960L), prediction_block=25L,
    audit_lattice_cells=as.integer(outer(c(1,11,20,30,40),
                                          c(1,5,10,14,18,23,27,31,36,40),
                                          function(x,y) x+40*(y-1))))
}

check_input_mode <- function(mode) {
  if (length(mode)!=1L || !mode %in% c("reconstructed","restored"))
    stop("An explicit input mode (reconstructed or restored) is required")
  mode
}
comparison_manifest <- function(config=study_config(), input_mode) {
  check_input_mode(input_mode)
  rows <- list()
  for (p in names(config$arms)) for (a in config$arms[[p]])
    for (comm in config$communities) for (arr in config$arrangements)
      rows[[length(rows)+1L]] <- data.frame(package=p,arm=a,community=comm,arrangement=arr)
  m <- do.call(rbind,rows); m$job_id <- seq_len(nrow(m))
  m$dataset_id <- paste(m$community,m$arrangement,sep="-")
  m$key <- paste(m$package,m$arm,m$dataset_id,sep="-")
  m$input_mode <- input_mode; m$new_fit <- TRUE
  m$fit_seed <- attempt_seed(m$job_id,1L,"fit")
  m$prediction_seed <- attempt_seed(m$job_id,1L,"prediction")
  m$diagnostic_seed <- attempt_seed(m$job_id,1L,"diagnostic")
  stopifnot(nrow(m)==192L,!anyDuplicated(m$key),sum(m$arm=="gp")==12L)
  m
}
oracle_manifest <- function(config=study_config(), input_mode) {
  check_input_mode(input_mode); rows <- list()
  for (comm in config$communities) for (arr in config$arrangements) for (s in config$species)
    rows[[length(rows)+1L]] <- data.frame(community=comm,arrangement=arr,species=s)
  m <- do.call(rbind,rows); m$job_id <- seq_len(nrow(m)); m$input_mode <- input_mode
  m$key <- paste(m$community,m$arrangement,m$species,sep="-")
  m$fit_seed <- attempt_seed(m$job_id,1L,"oracle")
  m$prediction_seed <- attempt_seed(m$job_id,1L,"oracle_prediction")
  stopifnot(nrow(m)==96L,!anyDuplicated(m$key)); m
}
attempt_seed <- function(job_id,attempt_id,kind) {
  stopifnot(all(attempt_id %in% 1:6),all(job_id>=1),all(job_id==as.integer(job_id)))
  base <- c(fit=27080000L,prediction=27180000L,diagnostic=27280000L,
            oracle=27380000L,oracle_prediction=27480000L)[kind]
  stopifnot(length(base)==1L,!is.na(base)); as.integer(base+100L*job_id+attempt_id)
}
