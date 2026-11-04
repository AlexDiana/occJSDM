# A version of simulateOccJSDMData() whose true occupancy intercepts are drawn
# around `b0_mean` instead of 0, for testing b0_prior = "community" on a
# community whose average intercept is not the fixed prior's 0. It rewrites
# the one line `B0 <- rnorm(S, sd = 1)` in simulateData() and errors unless it
# finds that line exactly once, so a change to the simulator cannot silently
# leave the intercepts unshifted. The simulation-study harness has the same
# function as simstudy_shifted_simulator().
shifted_b0_simulator <- function(b0_mean) {
  stopifnot(is.numeric(b0_mean), length(b0_mean) == 1L, is.finite(b0_mean))
  ns <- asNamespace("occJSDM")
  target <- quote(B0 <- rnorm(S, sd = 1))
  replacement <- bquote(B0 <- rnorm(S, sd = 1) + .(b0_mean))
  hits <- 0L
  swap <- function(e) {
    if (identical(e, target)) {
      hits <<- hits + 1L
      return(replacement)
    }
    # Only calls can contain the target. Testing e[[i]] in place, rather than
    # binding it, keeps empty arguments such as the one in x[, 1] harmless.
    for (i in seq_along(e)) if (is.call(e[[i]])) e[[i]] <- swap(e[[i]])
    e
  }
  inner <- utils::removeSource(get("simulateData", envir = ns))
  body(inner) <- swap(body(inner))
  if (hits != 1L) {
    stop("shifted_b0_simulator(): expected the line `B0 <- rnorm(S, ",
         "sd = 1)` once in occJSDM's simulateData(), found it ", hits,
         " times. Has the simulator changed?")
  }
  environment(inner) <- ns
  env <- new.env(parent = ns)
  assign("simulateData", inner, envir = env)
  outer <- get("simulateOccJSDMData", envir = ns)
  environment(outer) <- env
  outer
}
