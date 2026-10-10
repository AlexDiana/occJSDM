library(testthat)
code <- normalizePath("..")
for(f in c("config.R","archive.R","basis.R","adapter-contract.R","score.R","attribution.R")) source(file.path(code,f))
if(file.exists("../verify-teaching.R")) source("../verify-teaching.R")
fixture_bundle <- function() {
  cfg <- study_config(); m <- comparison_manifest(cfg,"reconstructed")
  landscapes <- training <- prediction_data <- truth_index <- list()
  for(comm in cfg$communities) {
    S <- length(cfg$species); xy <- cbind(seq(0,1,length.out=6),seq(0,1,length.out=6))
    land <- list(points=xy,environment=seq(-1,1,length.out=6),B0=rep(-1,S),B=rep(.4,S),species=cfg$species,
      field=outer(seq(-1,1,length.out=6),rep(1,S)),field_sd=1,range=.03,prevalence=c(.05,.05,.25,.25,.25,.75,.75,.75))
    land$psi <- plogis(-1+land$environment %o% land$B+land$field)
    land$z <- matrix(rep(c(0,0,0,1,0,1),S),6,S,dimnames=list(NULL,cfg$species)); landscapes[[comm]] <- land
    for(arr in cfg$arrangements) {
      ds <- paste(comm,arr,sep="-"); ids <- paste0(comm,"-",1:6)
      tr <- list(schema_version=1L,dataset_id=ds,site_id=ids[1:3],species_id=cfg$species,y=land$z[1:3,,drop=FALSE],environment=matrix(land$environment[1:3],ncol=1),xy=xy[1:3,,drop=FALSE],input_hashes=c(synthetic="test"))
      training[[ds]] <- tr
      pd <- tr; pd$y <- NULL; pd$site_id <- ids[4:6]; pd$environment <- matrix(land$environment[4:6],ncol=1); pd$xy <- xy[4:6,,drop=FALSE]; prediction_data[[ds]] <- pd
      truth_index[[ds]] <- list(community=comm,species_id=cfg$species,training_index=1:3,lattice_index=4:6,training_site_id=ids[1:3],prediction_site_id=ids[4:6])
    }
  }
  accounting <- transform(m,kind="comparison",species=NA_character_,status="failed",training_supported=FALSE,lattice_supported=FALSE)
  key <- m$key[which(m$package=="occJSDM" & m$arm=="gp")[1]]; row <- match(key,m$key); ds <- m$dataset_id[row]
  tv <- truth_index[[ds]]; tv$landscape <- landscapes[[tv$community]]
  view <- truth_view(tv,"training"); p <- view$probabilities
  accounting[row,c("training_supported","lattice_supported")] <- c(TRUE,FALSE); accounting$status[row] <- "selected"
  s <- score_predictions(p,view)$species
  s$key <- key; s$scope <- "training"; s$community <- m$community[row]; s$arrangement <- m$arrangement[row]; s$package <- "occJSDM"; s$arm <- "gp"
  list(schema_version=1L,config=cfg,landscapes=landscapes,training=training,prediction_data=prediction_data,truth_index=truth_index,
    predictions=setNames(list(list(training=p,lattice=NULL)),key),scores=list(accounting=accounting,species=s),input_mode="reconstructed")
}
test_that("archive-free verification protects all communities and binary observations", {
  b <- fixture_bundle(); expect_true(all(verify_teaching_bundle(b)$passed))
  bad <- b; bad$landscapes$rep03 <- NULL
  expect_error(verify_teaching_bundle(bad),"three communities")
  bad <- b; bad$training[[1]]$y <- NULL
  expect_error(verify_teaching_bundle(bad),"binary observations")
})
test_that("species order and declared support cannot drift during export", {
  b <- fixture_bundle(); bad <- b; bad$predictions[[1]]$training$species_id <- rev(b$config$species)
  expect_error(verify_teaching_bundle(bad),"species order")
  bad <- b; bad$scores$accounting$lattice_supported[1] <- TRUE
  expect_error(verify_teaching_bundle(bad),"support")
  bad <- b; bad$predictions[[1]]$training$uncertainty$status <- "available"
  expect_error(verify_teaching_bundle(bad),"interval")
})
test_that("independent compact arithmetic detects inconsistent recovery summaries", {
  b <- fixture_bundle(); b$scores$species$rmse[1] <- 10
  expect_error(verify_teaching_bundle(b),"rmse")
})
test_that("native species dimnames preserve exact correlation status checks", {
  b <- fixture_bundle(); p <- b$predictions[[1]]$training
  for(nm in c("p_full","p_space_zero","p_env_zero","p_both_zero","delta_p","field_mean","field_median")) {
    colnames(p[[nm]]) <- p$species_id
  }
  b$predictions[[1]]$training <- p
  expect_true(all(verify_teaching_bundle(b)$passed))
  bad <- b; bad$scores$species$delta_correlation_status[1] <- "constant fitted field"
  expect_error(verify_teaching_bundle(bad),"correlation status")
})
test_that("near-zero recovery gains use absolute arithmetic tolerance", {
  b <- fixture_bundle(); key <- names(b$predictions)[1]
  ds <- b$scores$accounting$dataset_id[match(key,b$scores$accounting$key)]
  tv <- b$truth_index[[ds]]; tv$landscape <- b$landscapes[[tv$community]]
  view <- truth_view(tv,"training"); p <- b$predictions[[key]]$training
  p$field_mean <- p$field_mean*1e-8
  b$predictions[[key]]$training <- p
  s <- score_predictions(p,view)$species
  for(nm in names(s)) b$scores$species[[nm]] <- s[[nm]]
  b$scores$species$field_rmse_reduction[1] <- b$scores$species$field_rmse_reduction[1]+1e-15
  b$scores$species$field_fractional_rmse_reduction[1] <- b$scores$species$field_fractional_rmse_reduction[1]+1e-15
  expect_true(all(verify_teaching_bundle(b)$passed))
  bad <- b; bad$scores$species$field_rmse_reduction[1] <- bad$scores$species$field_rmse_reduction[1]+1e-8
  expect_error(verify_teaching_bundle(bad),"rmse_reduction")
  bad <- b; bad$scores$species$field_rmse_reduction[1] <- NA_real_
  expect_error(verify_teaching_bundle(bad),"rmse_reduction")
})
test_that("group denominators retain two rare and three other species per community", {
  b <- fixture_bundle(); rows <- lapply(c(.05,.25,.75),function(prevalence) {
    s <- b$scores$species[b$scores$species$prevalence==prevalence,,drop=FALSE]
    summary <- equal_community_summary(s,"rmse",3L)
    summary$rows_expected <- 3L*(if(prevalence==.05) 2L else 3L)
    cbind(data.frame(package="occJSDM",arm="gp",arrangement="spread",scope="training",prevalence,
      cells_available=1L,cells_expected=3L),summary)
  })
  b$scores$groups <- do.call(rbind,rows)
  expect_true(all(verify_teaching_bundle(b)$passed))
  bad <- b; bad$scores$groups$rows_expected[2] <- 27L
  expect_error(verify_teaching_bundle(bad),"expected species rows")
})
test_that("compact verification checks per-site attribution and distance-bin arithmetic", {
  b <- fixture_bundle(); key <- names(b$predictions)[1]; ds <- b$scores$accounting$dataset_id[match(key,b$scores$accounting$key)]
  tv <- b$truth_index[[ds]]; tv$landscape <- b$landscapes[[tv$community]]
  score <- score_predictions(b$predictions[[key]]$training,truth_view(tv,"training"))
  b$site_attribution <- setNames(list(list(training=score$sites)),key)
  b$scores$strata <- transform(score$strata,key=key,scope="training")
  expect_true(all(verify_teaching_bundle(b)$passed))
  bad <- b; bad$site_attribution[[key]]$training$expected$space[1] <- .01
  expect_error(verify_teaching_bundle(bad),"attribution")
  bad <- b; bad$scores$strata$rmse[1] <- 5
  expect_error(verify_teaching_bundle(bad),"stratum")
})
if(file.exists("../export-teaching.R")) source("../export-teaching.R")
test_that("export refuses archives without independent full-run verification", {
  root <- tempfile(); dir.create(root)
  expect_error(build_teaching_bundle(root,tempdir()),"independent full-run verification")
})
test_that("Moran reporting separates graph distance from basis amplitude", {
  set.seed(441); xy <- matrix(runif(80),40,2)
  training <- list(dataset_id="fixture",site_id=paste0("t",1:40),xy=xy,environment=matrix(rnorm(40),ncol=1))
  basis <- fit_basis(training,"mem")
  prediction <- list(site_id=paste0("n",1:3),xy=rbind(c(0,0),c(.5,.5),c(1,1)))
  x <- moran_reporting_summary(basis,prediction)
  expect_equal(x$basis$training_sd,apply(basis$training_matrix,2,sd),ignore_attr=TRUE)
  expect_equal(x$graph$k,8L); expect_equal(x$graph$components,1L)
  expect_equal(x$distance$maximum_nearest_training_distance,max(apply(xy_distances(prediction$xy,xy),1,min)))
  expect_false(x$fitted_range_available)
})
test_that("saved bundle checksum is independent of the reader's R serializer", {
  b <- fixture_bundle(); path <- tempfile(fileext=".rds")
  save_teaching_bundle(b,path)
  expect_true(file.exists(paste0(path,".md5")))
  expect_true(all(verify_teaching_bundle(path)$passed))
  con <- file(path,"ab"); writeBin(as.raw(1),con); close(con)
  expect_error(verify_teaching_bundle(path),"file checksum")
})
source("../runner.R")
source("export-fixture.R")
test_that("a compact export retains full inputs, all cells and independent arithmetic", {
  fixture <- fixture_bundle(); root <- fixture_export_archive(fixture,code,normalizePath(file.path(code,"../../..")))
  on.exit(unlink(root,recursive=TRUE))
  before <- loadedNamespaces()
  bundle <- build_teaching_bundle(root,normalizePath(file.path(code,"../../..")))
  expect_length(bundle$landscapes,3L); expect_length(bundle$training,12L)
  expect_equal(nrow(bundle$scores$accounting),288L); expect_length(bundle$predictions,288L)
  expect_true(all(bundle$verification$compact$passed))
  bundle$scores$phase <- "full"
  bundle$scores$species <- complete_species_grid(bundle$scores$species,bundle$scores$accounting,
    setNames(bundle$landscapes$rep01$prevalence,bundle$config$species))
  expect_true(all(verify_teaching_bundle(bundle)$passed))
  bad <- bundle; bad$scores$species <- bad$scores$species[-nrow(bad$scores$species),]
  expect_error(verify_teaching_bundle(bad),"species denominator")
  expect_length(intersect(setdiff(loadedNamespaces(),before),c("occJSDM","sjSDM","Hmsc","gllvm")),0L)
  target <- file.path(root,"teaching.rds"); save_teaching_bundle(bundle,target)
  expect_true(all(verify_teaching_bundle(target)$passed))
})
test_that("lesson sources teach the spatial comparison without fitting while knitting", {
  repo <- normalizePath(file.path(code,"../../.."))
  expect_true(all(verify_spatial_lesson_sources(repo)$passed))
})
test_that("export CLI creates a verified teaching artifact without running a fit", {
  repo <- normalizePath(file.path(code,"../../.."))
  root <- fixture_export_archive(fixture_bundle(),code,repo)
  on.exit(unlink(root,recursive=TRUE))
  path <- file.path(root,"cli-teaching.rds")
  output <- system2(file.path(R.home(),"bin/exec/R"),c("--vanilla","--slave","-f",shQuote(file.path(code,"run.R")),"--args",
    "--mode=export","--input-mode=reconstructed",shQuote(paste0("--run=",root)),shQuote(paste0("--repo=",repo)),shQuote(paste0("--bundle=",path))),stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(output,"status"))) cat(paste(output,collapse="\n"),"\n")
  expect_null(attr(output,"status"))
  expect_true(file.exists(path))
  expect_true(all(verify_teaching_bundle(path)$passed))
})
