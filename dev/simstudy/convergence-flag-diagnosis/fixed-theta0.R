# Research-only fitting clone for variant (b) of Task 3 (PLAN.md and README.md
# in this directory): one species' field false-positive probability theta0 is
# held at a fixed value. Production code (R/, src/) is not changed; run.R
# sources this file and hashes it into every fit.
#
# make_fixed_theta0_fitter(species, value)
#   A copy of the installed occJSDM::runOccJSDM whose enclosing environment is
#   a child of the occJSDM namespace holding one binding, sample_theta0, which
#   shadows the package's own: the clone below. runOccJSDM's body, formals and
#   byte code are unchanged, and every other name it uses (sampler functions,
#   data preparation, output) resolves in the namespace exactly as before. The
#   fit therefore differs from runOccJSDM() only in the theta0 update, which
#   runOccJSDM calls once per iteration (R/runOccJSDM.R line 1258 at 707540a,
#   `theta0 <- sample_theta0(z, w, idx_z_w, a_theta0, b_theta0)`).
#
# sample_theta0_fixed
#   The package's sample_theta0() (R/mcmcfun.R lines 120-138 at 707540a) with
#   one added statement, marked below. Every species' Beta draw is taken
#   exactly as in the package, so R's random-number stream, and with it every
#   other species' value in that update, is unchanged; the fixed species' draw
#   is then replaced by the fixed value. Its deparsed text differs from the
#   package function's by that one statement (clone_matches_package()), and both
#   texts and their md5s are archived in results/fixed-theta0/.
#   Not changed: runOccJSDM's starting value (theta0 = 0.05 for every species
#   before the first iteration, R/runOccJSDM.R line 1115), so each chain's first
#   z and w updates use 0.05 for the fixed species; every later update and every
#   retained draw uses the fixed value.

# The added statement as deparse() prints it (two lines), without indentation.
FIXED_THETA0_ADDED_LINES <- c('if (s == fixed_species)','theta0[s] <- fixed_value')

sample_theta0_fixed <- function(z, w, idx_z, a_theta0, b_theta0){

  S <- ncol(z)

  z_all <- z[idx_z,]

  theta0 <- rep(NA, S)

  for (s in 1:S) {

    z0_1 <- sum(z_all[,s] == 0 & w[,s] == 1)
    z0_0 <- sum(z_all[,s] == 0 & w[,s] == 0)

    theta0[s] <- rbeta(1, a_theta0 + z0_1, b_theta0 + z0_0)

    # Research-only change (variant b): hold one species' theta0 fixed.
    if (s == fixed_species) theta0[s] <- fixed_value

  }

  theta0
}

package_function <- function(name) get(name,envir=asNamespace('occJSDM'),inherits=FALSE)

# TRUE when the clone's deparsed text is the package function's with exactly
# the added statement's two lines inserted, and the formals are identical.
clone_matches_package <- function(package_fun=package_function('sample_theta0')) {
  clone <- deparse(sample_theta0_fixed);pkg <- deparse(package_fun)
  at <- which(trimws(clone)==FIXED_THETA0_ADDED_LINES[1])
  if(length(at)!=1L || at==length(clone)) return(FALSE)
  added <- c(at,at+1L)
  identical(trimws(clone[added]),FIXED_THETA0_ADDED_LINES) && identical(clone[-added],pkg) &&
    identical(formals(sample_theta0_fixed),formals(package_fun))
}

make_fixed_sample_theta0 <- function(species,value) {
  if(length(species)!=1L || !is.numeric(species) || species!=round(species) || species<1)
    stop('species must be one positive whole number')
  if(length(value)!=1L || !is.finite(value) || value<=0 || value>=1) stop('value must be one number in (0, 1)')
  env <- new.env(parent=asNamespace('occJSDM'))
  env$fixed_species <- as.integer(species);env$fixed_value <- as.numeric(value)
  f <- sample_theta0_fixed;environment(f) <- env
  f
}

# runOccJSDM with the theta0 update replaced by `update`; nothing else changes.
make_theta0_override_fitter <- function(update,target=package_function('runOccJSDM')) {
  env <- new.env(parent=asNamespace('occJSDM'))
  env$sample_theta0 <- update
  f <- target;environment(f) <- env
  f
}

make_fixed_theta0_fitter <- function(species,value,target=package_function('runOccJSDM')) {
  if(!clone_matches_package()) stop('sample_theta0_fixed no longer matches the installed sample_theta0 plus one line')
  make_theta0_override_fitter(make_fixed_sample_theta0(species,value),target)
}

# ---- Archive of the clone ---------------------------------------------------------

text_md5 <- function(lines) {
  f <- tempfile();on.exit(unlink(f),add=TRUE)
  writeLines(lines,f);unname(tools::md5sum(f))
}
FIXED_THETA0_TEXT_FILES <- c(package='sample_theta0-package.txt',clone='sample_theta0-clone.txt')

fixed_theta0_texts <- function(package_fun=package_function('sample_theta0'))
  list(package=deparse(package_fun),clone=deparse(sample_theta0_fixed))

write_fixed_theta0_archive <- function(dir,revision,package_fun=package_function('sample_theta0')) {
  if(!clone_matches_package(package_fun)) stop('The clone does not match the package function plus one line')
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  texts <- fixed_theta0_texts(package_fun)
  for(k in names(FIXED_THETA0_TEXT_FILES)) writeLines(texts[[k]],file.path(dir,FIXED_THETA0_TEXT_FILES[[k]]))
  utils::write.csv(data.frame(item=c('package_sample_theta0','clone_sample_theta0_fixed','added_lines'),
    file=c(FIXED_THETA0_TEXT_FILES[['package']],FIXED_THETA0_TEXT_FILES[['clone']],''),
    md5=c(text_md5(texts$package),text_md5(texts$clone),text_md5(FIXED_THETA0_ADDED_LINES)),
    revision=revision,stringsAsFactors=FALSE),file.path(dir,'hashes.csv'),row.names=FALSE)
  invisible(dir)
}

# The archived md5s, after checking that the installed package function and
# the clone are exactly the archived texts. Refuses otherwise.
check_fixed_theta0_archive <- function(dir,package_fun=package_function('sample_theta0')) {
  f <- file.path(dir,'hashes.csv')
  if(!file.exists(f)) stop('Missing fixed-theta0 archive: ',f)
  h <- utils::read.csv(f,colClasses='character')
  texts <- fixed_theta0_texts(package_fun)
  want <- stats::setNames(h$md5,h$item)
  found <- c(package_sample_theta0=text_md5(texts$package),clone_sample_theta0_fixed=text_md5(texts$clone),
    added_lines=text_md5(FIXED_THETA0_ADDED_LINES))
  if(!identical(sort(names(want)),sort(names(found))) || !identical(unname(want[names(found)]),unname(found)))
    stop('The fixed-theta0 clone or the installed sample_theta0 differs from the archive in ',dir)
  for(k in names(FIXED_THETA0_TEXT_FILES)) if(!identical(readLines(file.path(dir,FIXED_THETA0_TEXT_FILES[[k]])),texts[[k]]))
    stop('Archived ',FIXED_THETA0_TEXT_FILES[[k]],' differs from the current function text')
  if(!clone_matches_package(package_fun)) stop('The clone does not match the package function plus one line')
  found
}
