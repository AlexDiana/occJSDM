args <- commandArgs(trailingOnly=TRUE); code <- normalizePath(args[1]); repo <- normalizePath(file.path(code,"../../..")); root <- args[2]
for(f in c("config.R","archive.R","basis.R","adapter-contract.R","runner.R")) source(file.path(code,f))
dir.create(root,recursive=TRUE,showWarnings=FALSE)
set.seed(504); n <- 20L; S <- 8L
tr <- list(schema_version=1L,dataset_id="rep01-spread",site_id=sprintf("smoke%03d",1:n),species_id=sprintf("species%02d",1:S),
  environment=matrix(rnorm(n),n,1,dimnames=list(NULL,"environment")),xy=matrix(runif(2*n),n,2),input_hashes=c(synthetic="CLI fixture"))
tr$y <- matrix(rbinom(n*S,1,.4),n,S,dimnames=list(NULL,tr$species_id))
pred <- tr; pred$y <- NULL; pred$site_id <- sprintf("new%03d",1:n); pred$xy <- matrix(runif(2*n),n,2)
land <- list(points=rbind(tr$xy,pred$xy),environment=c(tr$environment,pred$environment),
  B0=rep(-1,S),B=rep(.4,S),range=.1,field_sd=1,species=tr$species_id,prevalence=c(.05,.05,.25,.25,.25,.75,.75,.75),field=matrix(rnorm(2*n*S),2*n,S))
land$psi <- plogis(matrix(rep(land$B0,each=2*n),2*n,S)+land$environment %o% land$B+land$field)
land$z <- matrix(rbinom(2*n*S,1,land$psi),2*n,S)
tr$y <- land$z[1:n,,drop=FALSE]; colnames(tr$y) <- tr$species_id
truth <- list(landscape=land,training_index=1:n,lattice_index=n+1:n,training_site_id=tr$site_id,prediction_site_id=pred$site_id,species_id=tr$species_id)
truth_path <- file.path(root,"truth/rep01-spread.rds"); write_once(truth,truth_path)
write_once(data.frame(dataset_id="rep01-spread",kind="truth",relative_path="truth/rep01-spread.rds",md5=md5(truth_path)),file.path(root,"import-manifest.rds"))
for(kind in c("training","prediction")) write_once(if(kind=="training") tr else pred,file.path(root,kind,"rep01-spread.rds"))
cfg <- study_config(); cfg$mcmc$initial <- cfg$mcmc$extended <- c(nchain=2L,nburn=10L,niter=20L,nthin=1L)
cfg$sjsdm$epochs <- c(2L,2L); cfg$sjsdm$sampling <- 16L; cfg$gllvm$maxit <- 100L; cfg$smoke_only <- TRUE
write_once(cfg,file.path(root,"config.rds"))
m <- comparison_manifest(study_config(),"reconstructed")
m <- m[m$community=="rep01" & m$arrangement=="spread" & ((m$package=="occJSDM" & m$arm=="gp") | (m$package=="sjSDM" & m$arm=="trend") | (m$package %in% c("Hmsc","gllvm") & m$arm=="spatial2")),]
write_once(m,file.path(root,"comparison-manifest.rds"))
o <- oracle_manifest(study_config(),"reconstructed"); o <- o[o$community=="rep01"&o$arrangement=="spread"&o$species=="species01",]
write_once(o,file.path(root,"oracle-manifest.rds"))
for(arm in unique(m$arm)) write_once(fit_basis(tr,arm),file.path(root,"bases",paste0("rep01-spread-",arm,".rds")))
minimal <- list(landscape=list(B0=rep(-1,S),B=rep(.4,S),range=.1),surveys=list(spread=list(truth=list(xy=tr$xy,environment=as.vector(tr$environment),z=tr$y))))
write_once(minimal,file.path(root,"oracle-inputs/rep01-spread.rds"))
main <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-package-comparison-20261008"
dir.create(file.path(root,"runtime")); stopifnot(file.copy(file.path(main,"runtime/package-records.rds"),file.path(root,"runtime")))
dir.create(file.path(root,"sources")); stopifnot(file.copy(file.path(main,"sources/acceptance"),file.path(root,"sources"),recursive=TRUE))
freeze_sources(code,repo,root)
write_once(list(smoke_only=TRUE),file.path(root,"preflight.rds"))
cat("Synthetic CLI root",root,"\n")
