# Prior-invariance run of the scale draws plus both scale moves (sigma_emb, then the axis scales),
# with a constant likelihood. Measures the Monte Carlo error of the quartile shares over seeds, and
# what two wrong acceptance ratios in the sigma_emb move do to them. usage: Rscript <this> <lib>
.libPaths(c(commandArgs(TRUE)[1], .libPaths())); suppressMessages(library(occJSDM))
seed_both <- function(s) { set.seed(s); occJSDM:::setOccJSDMSeed(as.integer(s)) }
run <- function(sigma_move, seed, N = 50000L, burn = 1000L) {
  S <- 3L; p <- 3L; gt <- 2L; g_b <- 2L; n <- 4L
  seed_both(seed)
  Tr_b <- cbind(1, rnorm(S)); X <- matrix(rnorm(n * p), n, p); zero <- matrix(0, n, S)
  cp <- occJSDM:::read_coef_prior(list(sigma_b_scale = 0.8, sigma_emb_scale = 1.3, axis_scale = 0.5),
    list(embeddingVars = c("e1", "e2")), S, gt,
    data.frame(column = c("m1", "e1", "e2"), variable = c("m1", "e1", "e2"), kind = "numeric"))
  start <- occJSDM:::coef_scales_start(cp, p)
  sigma_b <- cp$sigma_b_scale; sigma_emb <- start$sigma_emb; lambda <- start$lambda
  G <- matrix(rnorm(g_b * p), g_b, p); C <- matrix(rnorm(gt * p), gt, p); A <- matrix(rnorm(S * gt), S, gt)
  keep <- matrix(NA_real_, N + burn, 4L); acc <- 0L
  ok <- tryCatch({ for (i in seq_len(N + burn)) {
    sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
    centre <- Tr_b %*% G + A %*% C
    B <- t(centre + matrix(rnorm(S * p), S, p) * rep(sds$b, each = S))
    sc <- occJSDM:::sample_coef_scales(t(B) - centre, G, C, cp, sigma_b, sigma_emb, lambda)
    sigma_b <- sc$sigma_b; sigma_emb <- sc$sigma_emb; lambda <- sc$lambda
    eta <- zero
    if (!is.null(sigma_move)) { sg <- sigma_move(B, t(B) - centre, sigma_emb, X, zero, zero, zero, cp, 0.5)
      B <- sg$B; sigma_emb <- sg$sigma_emb; eta <- sg$eta; acc <- acc + sg$accepted }
    mv <- occJSDM:::axis_scale_move(B, G, C, lambda, X, zero, zero, eta, cp, 0.5)
    B <- mv$B; G <- mv$G; C <- mv$C; lambda <- mv$lambda
    sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
    gc <- occJSDM:::sample_GC(B, Tr_b, A, sds$b, sds$g, sds$c); G <- gc$G; C <- gc$C
    A <- occJSDM:::sample_A(B, C, Tr_b, G, sds$b)$A
    keep[i, ] <- c(sigma_b / cp$sigma_b_scale, sigma_emb / cp$sigma_emb_scale, lambda[2:3] / cp$axis_scale)
  }; TRUE }, error = function(e) conditionMessage(e))
  if (!isTRUE(ok)) return(list(error = ok))
  keep <- keep[-seq_len(burn), ]
  q <- tan(pi / 2 * c(.25, .5, .75))
  share <- vapply(q, function(x) colMeans(keep < x), numeric(4)) - rep(c(.25, .5, .75), each = 4)
  list(err = share, acc = acc / (N + burn), ess = apply(log(keep), 2, function(v) coda::effectiveSize(coda::mcmc(v))))
}
real <- occJSDM:::sigma_emb_move
mutate <- function(from, to) { txt <- paste(deparse(real), collapse = "\n"); stopifnot(grepl(from, txt, fixed = TRUE))
  f <- eval(parse(text = sub(from, to, txt, fixed = TRUE))); environment(f) <- asNamespace("occJSDM"); f }
no_jacobian <- mutate("+ delta -", "-")
wrong_scale <- mutate("(proposed/coef_prior$sigma_emb_scale)^2", "(proposed/(2 * coef_prior$sigma_emb_scale))^2")
show <- function(label, r) { if (!is.null(r$error)) { cat(label, "ERROR:", r$error, "\n"); return(invisible()) }
  cat(sprintf("%-22s acc %.2f  max|err| sigma_b %.3f sigma_emb %.3f lambda %.3f %.3f   ESS(log) %s\n", label, r$acc,
      max(abs(r$err[1, ])), max(abs(r$err[2, ])), max(abs(r$err[3, ])), max(abs(r$err[4, ])), paste(round(r$ess), collapse = " "))) }
errs <- list()
for (s in 1:8) { r <- run(real, s); show(paste("both moves seed", s), r); errs[[s]] <- r$err }
a <- simplify2array(errs); cat("sd over seeds of each share error (rows sigma_b, sigma_emb, lambda_e1, lambda_e2):\n"); print(round(apply(a, c(1, 2), sd), 4))
cat("largest mean error:", round(max(abs(apply(a, c(1, 2), mean))), 4), "\nlargest |error| per row:\n"); print(round(apply(abs(a), 1, max), 4))
for (s in 1:2) show(paste("axis move only seed", s), run(NULL, s))
for (s in 1:2) show(paste("no Jacobian seed", s), run(no_jacobian, s))
for (s in 1:2) show(paste("wrong scale seed", s), run(wrong_scale, s))
