# Reproduce algebra checks from site-variation-partitioning.md without fitting
# or importing sjSDM/Python. Run from the occJSDM root:
# Rscript dev/simstudy/site-variation-partitioning-audit.R /path/to/s-jSDM
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)

read_functions <- function(path, names, environment) {
  for (expression in parse(path)) {
    if (is.call(expression) && as.character(expression[[1]]) %in% c("<-", "=") &&
        is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% names) {
      eval(expression, envir = environment)
    }
  }
  stopifnot(all(vapply(names, exists, logical(1), envir = environment,
                      inherits = FALSE)))
}

occ <- new.env(parent = globalenv())
occ$logistic <- stats::plogis
read_functions("R/jsdmfun.R",
               c("fast_sd", "SD", "computeVariancePartitioning_stdevs"), occ)
x <- matrix(seq(-2, 2, length.out = 10), ncol = 1)
zero <- matrix(0, 10, 1)
cancel <- occ$computeVariancePartitioning_stdevs(x, -x, zero, "binary")
constant <- occ$computeVariancePartitioning_stdevs(zero, zero, zero, "binary")
intercept <- occ$computeVariancePartitioning_stdevs(
  matrix(-3, 10, 1), matrix(seq(2, 4, length.out = 10), 10, 1), zero, "binary"
)
stopifnot(all.equal(cancel$Environmental, .5),
          all.equal(cancel$Spatial, .5), cancel$Total == 0,
          all(is.nan(unlist(constant[1:3]))), constant$Total == 0,
          intercept$Environmental > .35, intercept$Environmental < .45)
print(rbind(cancellation = cancel, constant = constant, intercept = intercept))

sj <- new.env(parent = globalenv())
read_functions(file.path(args[1], "sjSDM/R/anova.R"), "get_shared_anova", sj)
shared <- sj$get_shared_anova(
  list(A = .3, B = 0, S = .3, AB = .3, AS = .3, BS = .3, Full = .3), "3"
)$equal
stopifnot(abs(sum(unlist(shared[1:3])) - .1) < 1e-7,
          shared$F_B == 0, shared$F_S == 0)
print(shared)
cat("Current sjSDM equal allocation sums to 0.1 for a 0.3 pure E:S overlap.\n")

# The proposed order-neutral allocation: average contributions across all six
# orders in which three components can be added. This is a research check,
# not an exported occJSDM function or a scored model fit.
orders <- rbind(c("E", "S", "B"), c("E", "B", "S"),
                c("S", "E", "B"), c("S", "B", "E"),
                c("B", "E", "S"), c("B", "S", "E"))
allocate <- function(scores) {
  contribution <- c(E = 0, S = 0, B = 0)
  for (i in seq_len(nrow(orders))) {
    previous <- "none"
    included <- character()
    for (component in orders[i, ]) {
      included <- sort(c(included, component))
      current <- paste(included, collapse = "")
      contribution[component] <- contribution[component] +
        scores[current] - scores[previous]
      previous <- current
    }
  }
  contribution / nrow(orders)
}
overlap <- c(none = 0, E = .3, S = .3, B = 0,
             ES = .3, BE = .3, BS = .3, BES = .3)
stopifnot(isTRUE(all.equal(unname(allocate(overlap)), c(.15, .15, 0))))
identical_models <- setNames(rep(0, 8), names(overlap))
stopifnot(all(allocate(identical_models) == 0))
negative <- c(none = 0, E = .2, S = -.1, B = 0,
              ES = .1, BE = .2, BS = -.1, BES = .1)
stopifnot(isTRUE(all.equal(unname(allocate(negative)), c(.2, -.1, 0))))
set.seed(22)
for (i in 1:100) {
  scores <- setNames(rnorm(8), names(overlap))
  stopifnot(abs(sum(allocate(scores)) - (scores["BES"] - scores["none"])) < 1e-12)
}
cat("Proposed allocation conserves total improvement, splits shared effects symmetrically, and retains negative contributions.\n")
