# Parallel MCMC chains: core budget, per-chain seeds, the single-chain
# function, the PSOCK driver and the combination helpers.
# Design: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.

# Holds the prepared setup list while chains run, in the parent and in every
# worker. Package level, so that a worker can reach it by name after the task
# function arrives as a reference to the namespace.
chain_env <- new.env(parent = emptyenv())

# `available` is an argument so that the warning can be tested without
# depending on the operating system: parallel::detectCores() is documented to
# return NA when the count is unknown, and then no warning is possible.
resolve_ncores <- function(MCMCparams, available = parallel::detectCores()) {
  ncores <- get_param(MCMCparams, "ncores", 1)
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || ncores != floor(ncores)) {
    stop("MCMCparams$ncores must be a single whole number of at least 1", call. = FALSE)
  }
  ncores <- as.integer(ncores)
  if (!is.na(available) && ncores > available) {
    warning(sprintf("MCMCparams$ncores = %d exceeds the %d cores detected on this machine",
                    ncores, available), call. = FALSE)
  }
  ncores
}

# ncores is the total core budget: workers = min(nchain, ncores), each with
# floor(ncores / workers) RcppParallel threads. ncores = 1 is today's fit and
# leaves the session's thread setting alone ("ambient").
plan_chain_budget <- function(nchain, ncores) {
  workers <- as.integer(min(nchain, ncores))
  list(ncores = as.integer(ncores),
       workers = workers,
       threads_per_worker = max(1L, as.integer(floor(ncores / workers))),
       thread_setting = if (ncores == 1L) "ambient" else "budget")
}

ambient_rcpp_threads <- function() {
  env <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "")
  if (nzchar(env)) as.integer(env) else as.integer(RcppParallel::defaultNumThreads())
}

# Merge two Welford accumulators (running mean, M2 and mean likelihood per
# observation, each over n draws) exactly, with the standard pairwise formula.
# Chains are merged in chain order by merge_waic_parts(), so the result cannot
# depend on which worker produced which chain.
merge_waic_pair <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  n <- a$n + b$n
  delta <- b$mean_log - a$mean_log
  list(M2 = a$M2 + b$M2 + delta^2 * a$n * b$n / n,
       mean_log = a$mean_log + delta * b$n / n,
       mean_lik = (a$n * a$mean_lik + b$n * b$mean_lik) / n,
       n = n)
}

merge_waic_parts <- function(parts) {
  Reduce(merge_waic_pair, parts, accumulate = FALSE)
}

# Per-chain sums of the retained latent states, added in chain order. A name
# that is NULL or absent on either side comes out NULL, so a shape mismatch
# between chains cannot produce a silently partial sum. NULL entries are kept
# in the result (`[<-` with list(NULL), not `[[<-` with NULL), so the names
# stay the same however many parts are folded in.
add_latent_sums <- function(total, part) {
  if (is.null(part)) return(total)
  if (is.null(total)) return(part)
  nms <- union(names(total), names(part))
  out <- stats::setNames(vector("list", length(nms)), nms)
  for (nm in nms) {
    if (!is.null(total[[nm]]) && !is.null(part[[nm]])) out[[nm]] <- total[[nm]] + part[[nm]]
  }
  out
}

# The fields the chain body reads from runOccJSDM()'s frame. Computed on
# 13 October 2026 with codetools::findGlobals() on the loop body, plus the
# names that analysis misses because the body also assigns them: list_data,
# and the states z and w, which runOccJSDM() sets from y before the loop for
# the binary, continuous and occupancy models (R/runOccJSDM.R:601-604) and
# which the body only assigns for the other models; plus Tr and N3, read
# inside calls. A PSOCK worker enforces completeness: a missing field is an error.
chain_setup_fields <- function() c(
  "model", "jsdmModel", "threshold",
  "y", "y_NA", "z", "w", "M", "K", "P", "sumM", "sumP", "sumK", "maxP", "primerId_p",
  "idx_w_k", "idx_p_k", "idx_z_w", "list_idx",
  "n", "N", "N3", "S", "d", "g", "gt", "ps", "ncov_theta", "ncov_psi",
  "X_theta", "X_psi", "Tr",
  "list_Xs", "list_SoRSummaries",
  "list_data", "list_priors", "spatial_sd_prior",
  "a_theta0", "b_theta0", "a_p", "b_p", "a_q", "b_q", "b_betatheta", "B_betatheta",
  "nburn", "niter", "nthin", "summarisedLatentPresences")

set_chain_setup <- function(setup) {
  assign("setup", setup, envir = chain_env)
  invisible(NULL)
}

clear_chain_setup <- function() {
  if (exists("setup", envir = chain_env, inherits = FALSE)) rm("setup", envir = chain_env)
  invisible(NULL)
}

# Two seeds per chain, drawn in chain order so that chain k's seeds do not
# depend on how many chains follow it. Distinct by construction except for an
# astronomically unlikely collision, which is redrawn.
draw_chain_seeds <- function(nchain) {
  draw <- function() {
    seeds <- matrix(NA_integer_, nrow = 2L, ncol = nchain, dimnames = list(c("r", "cpp"), NULL))
    for (k in seq_len(nchain)) seeds[, k] <- sample.int(.Machine$integer.max, 2L)
    seeds
  }
  seeds <- draw()
  while (anyDuplicated(as.vector(seeds)) > 0L) seeds <- draw()
  seeds
}

save_session_rng <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) invisible(sample.int(1L))
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}

# Reset first, then restore: set.seed() with the pinned kinds clears R's
# cached Box-Muller deviate the same way in every path, and the saved state
# carries the user's kinds back with it. Registered with on.exit() in
# runOccJSDM() so it runs on every exit. Spec: "After the chains".
restore_session_rng <- function(saved) {
  set.seed(0L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  assign(".Random.seed", saved, envir = globalenv())
  invisible(NULL)
}

seed_chain_rng <- function(seed_r, seed_cpp) {
  set.seed(seed_r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  setOccJSDMSeed(seed_cpp)
  invisible(NULL)
}

chain_failure <- function(k, message) {
  structure(list(chain = k, message = message), class = "occJSDM_chain_failure")
}

stop_if_chain_failed <- function(res) {
  if (inherits(res, "occJSDM_chain_failure")) {
    stop(sprintf("Chain %d failed: %s", res$chain, res$message), call. = FALSE)
  }
  res
}

# The task a worker (or the in-process loop) runs for one chain. Package
# level, so it crosses the socket as a namespace reference; setup is read
# from chain_env, which the bootstrap filled.
run_chain_task <- function(k, seeds) {
  tryCatch(run_one_chain(get("setup", envir = chain_env, inherits = FALSE),
                         seeds["r", k], seeds["cpp", k], k),
           error = function(e) chain_failure(k, conditionMessage(e)))
}

# With handle = NULL this is the in-process path; the cluster branch is
# exercised from Task 4 onward.
run_chain_batch <- function(handle, batch, seeds, task = run_chain_task) {
  if (is.null(handle)) return(lapply(batch, task, seeds = seeds))
  parallel::clusterApplyLB(handle$cl, batch, task, seeds = seeds)
}

drop_execution <- function(fit) {
  fit$infos$execution <- NULL
  fit
}

# One chain, from its starting values to its last retained draw. Everything
# the body reads arrives in setup and is unpacked into this frame, so the
# moved code runs unchanged. The chain's randomness is fixed by its two seeds
# alone: it does not continue any other chain's stream, in process or not.
run_one_chain <- function(setup, seed_r, seed_cpp, chain) {
  started <- proc.time()[["elapsed"]]
  # Unpack setup, one field per line in the order of chain_setup_fields(), from
  # which this block was generated: keep the two in step when a field changes.
  # Explicit assignments rather than list2env(), so that R CMD check's code
  # analysis sees every name defined here. A NULL field assigns NULL.
  model <- setup$model
  jsdmModel <- setup$jsdmModel
  threshold <- setup$threshold
  y <- setup$y
  y_NA <- setup$y_NA
  z <- setup$z
  w <- setup$w
  M <- setup$M
  K <- setup$K
  P <- setup$P
  sumM <- setup$sumM
  sumP <- setup$sumP
  sumK <- setup$sumK
  maxP <- setup$maxP
  primerId_p <- setup$primerId_p
  idx_w_k <- setup$idx_w_k
  idx_p_k <- setup$idx_p_k
  idx_z_w <- setup$idx_z_w
  list_idx <- setup$list_idx
  n <- setup$n
  N <- setup$N
  N3 <- setup$N3
  S <- setup$S
  d <- setup$d
  g <- setup$g
  gt <- setup$gt
  ps <- setup$ps
  ncov_theta <- setup$ncov_theta
  ncov_psi <- setup$ncov_psi
  X_theta <- setup$X_theta
  X_psi <- setup$X_psi
  Tr <- setup$Tr
  list_Xs <- setup$list_Xs
  list_SoRSummaries <- setup$list_SoRSummaries
  list_data <- setup$list_data
  list_priors <- setup$list_priors
  spatial_sd_prior <- setup$spatial_sd_prior
  a_theta0 <- setup$a_theta0
  b_theta0 <- setup$b_theta0
  a_p <- setup$a_p
  b_p <- setup$b_p
  a_q <- setup$a_q
  b_q <- setup$b_q
  b_betatheta <- setup$b_betatheta
  B_betatheta <- setup$B_betatheta
  nburn <- setup$nburn
  niter <- setup$niter
  nthin <- setup$nthin
  summarisedLatentPresences <- setup$summarisedLatentPresences
  seed_chain_rng(seed_r, seed_cpp)

  # ---- moved from runOccJSDM(): chain output, starting values, iteration loop ----
  # (lines 1055 to 1366 of R/runOccJSDM.R before this change, edited per the plan)
    # chain output
    {
      if(model %in% c("occupancy","two_stage")){
        beta_theta_output_chain <- array(NA, dim = c(ncov_theta, S, niter))
        theta0_output_chain <- array(NA, dim = c(S, niter))
        if(!summarisedLatentPresences){
          z_output_chain <- array(NA, dim = c(n, S, niter))
          psi_output_chain <- array(NA, dim = c(n, S, niter))
        }
      }

      if(model == "two_stage"){
        p_output_chain <- array(NA, dim = c(maxP, S, niter))
        q_output_chain <- array(NA, dim = c(maxP, S, niter))
      }

      # jsdm params
      {
        B0_output_chain <- array(NA, dim = c(S, niter))
        B_output_chain <- array(NA, dim = c(ncov_psi, S, niter))
        G_output_chain <- array(NA, dim = c(g, ncov_psi, niter))
        A_output_chain <- array(NA, dim = c(S, gt, niter))
        C_output_chain <- array(NA, dim = c(gt, ncov_psi, niter))
        Bs_output_chain <- array(NA, dim = c(ps, S, niter))
        Gs_output_chain <- array(NA, dim = c(g, ps, niter))
        As_output_chain <- array(NA, dim = c(S, gt, niter))
        Cs_output_chain <- array(NA, dim = c(gt, ps, niter))
        U_output_chain <- array(NA, dim = c(n, d, niter))
        L_output_chain <- array(NA, dim = c(d, S, niter))
        sigmab_output_chain <- rep(NA, niter)
        sigmabs_output_chain <- rep(NA, niter)
        sigmah_output_chain <- rep(NA, niter)
        idx_ls_output_chain <- rep(NA, niter)
        tau_output_chain <- matrix(NA, S, niter)
        varPart_output_chain <-  array(NA, dim = c(S, 4, niter))
      }

    }

    # Per-chain sums of the retained latent states and per-chain WAIC
    # accumulators. The parent combines them in chain order.
    z_sum <- psi_sum <- theta_sum <- w_sum <- NULL
    if (model %in% c("occupancy", "two_stage")) {
      if (summarisedLatentPresences) {
        z_sum <- matrix(0, n, S)
        psi_sum <- matrix(0, n, S)
      }
      theta_sum <- matrix(0, N, S)
    }
    if (model == "two_stage") w_sum <- matrix(0, N, S)
    list_waic_jsdm <- create_waic_quantities(n * S)
    list_waic_w <- if (model %in% c("occupancy", "two_stage")) create_waic_quantities(N * S) else NULL
    list_waic_y <- if (model == "two_stage") create_waic_quantities(N3 * S) else NULL
    currentWAICiter <- 1

    # starting values
    {
      if(model == "two_stage"){

        w <- matrix(NA, N, S)

        for (s in 1:S) {
          for (i in 1:n) {
            for (m in 1:M[i]) {
              idx_im1 <- sumK[sumP[sumM[i] + m] + 1] + 1
              idx_im2 <- sumK[sumP[sumM[i] + m] + P[sumM[i] + m]] +
                K[sumP[sumM[i] + m] + P[sumM[i] + m]]
              y_subset <- y[idx_im1:idx_im2,s]
              y_subset <- y_subset[!is.na(y_subset)]
              if(length(y_subset) > 0){
                w[sumM[i] + m,s] <- as.numeric(any(y_subset > 0))
              } else {
                w[sumM[i] + m,s] <- 1
              }
            }
          }
        }

        c_imk <- array(as.numeric(y > 0), dim = dim(y), dimnames = dimnames(y))

        p <- matrix(.9, maxP, S)
        q <- matrix(.05, maxP, S)

        y_pos <- (y > 0)

      }

      if(model %in% c("occupancy","two_stage")){

        z <- matrix(NA, n, S)

        for (s in 1:S) {
          for (i in 1:n) {
            idx_i1 <- sumM[i] + 1
            idx_i2 <- sumM[i] + M[i]

            z[i,s] <- as.numeric(any(w[idx_i1:idx_i2, s] > 0))
          }
        }

        beta_theta <- matrix(0, ncov_theta, S)
        theta <- computeTheta(X_theta, beta_theta)
        theta0 <- rep(.05, S)
      }

      # jsdmParams
      {
        B0 <- rep(0, S)
        B <- matrix(0, ncov_psi, S)
        Bs <- matrix(0, ps, S)
        L <- matrix(1, d, S)
        diag(L) <- 1
        L[lower.tri(L)] <- 0
        G <- matrix(0, g, ncov_psi)
        C <- matrix(1, gt, ncov_psi)
        diag(C) <- 1
        C[lower.tri(C)] <- 0
        A <- matrix(0, S, gt)
        Gs <- matrix(0, g, ps)
        Cs <- matrix(1, gt, ps)
        diag(Cs) <- 1
        Cs[lower.tri(Cs)] <- 0
        As <- matrix(0, S, gt)
        U <- matrix(0, n, d)
        sigma_b <- 1
        sigma_bs <- if (identical(spatial_sd_prior$type,"half_cauchy"))
          spatial_sd_prior$scale else .001
        sigma_h <- 1
        idx_ls <- 3 # dim(list_SoRSummaries$Ks_all)[3][5]
        tau <- rep(1, S)

        Bt <- t(B) - computeBtcoef(G, Tr, A, C, matrix(0, S, ncov_psi))
        Bst <- t(Bs) - computeBtcoef(Gs, Tr, As, Cs, matrix(0, S, ps))

        Ks <- matrix(list_SoRSummaries$Ks_all[,,idx_ls],nrow=n)

        list_jSDMparams <- list(
          "B0" = B0,
          "B" = B,
          "G" = G,
          "A" = A,
          "C" = C,
          "Bt" = Bt,
          "Bs" = Bs,
          "Gs" = Gs,
          "As" = As,
          "Cs" = Cs,
          "Bst" = Bst,
          "U" = U,
          "L" = L,
          "sigma_b" = sigma_b,
          "sigma_bs" = sigma_bs,
          "sigma_h" = sigma_h,
          "idx_ls" = idx_ls,
          "tau" = tau
        )

        list_psiCoef <- computePsiCoef(
          X_psi, Ks, list_Xs$Xs_centers, Tr,
          B0, G, A, C, Bt,
          Gs, As, Cs, Bst,
          U, L)

      }

      psi <- logistic(list_psiCoef$eta)

    }

    for (iter in 1:(nburn + niter * nthin)) {

      if(iter > nburn){
        currentIter <- (iter - nburn) / nthin
        if(currentIter %% 100 == 0){
          message(paste0("Chain ", chain, " - Iteration ",currentIter))
        }
      } else {
        if(iter %% 100 == 0){
          message(paste0("Chain ", chain, " - Burn in Iteration ",iter))
        }
      }

      # sample z
      if(model %in% c("occupancy","two_stage")){
        z <- sample_z_cpp_parallel(w, psi, theta, theta0, M, sumM)
      }

      # sample jsdm coef
      {
        list_data$z <- z

        list_jSDMparams <- update_jSDMcoef(
          list_data,
          list_jSDMparams,
          list_priors,
          list_Xs,
          list_SoRSummaries,
          model = jsdmModel,
          computeVarPart = (iter > nburn & (iter - nburn) %% nthin == 0)
        )
        psi <- list_jSDMparams$psi
        eta <- list_jSDMparams$eta
      }

      # sample w
      if(model == "two_stage"){
        if(threshold > 0){

          w <- sample_w_cim_cipp_parallel(y, y_NA, theta, theta0, p, q,
                                 M, K, sumP, sumM, sumK, P, primerId_p, maxP, z)

        }
      }

      # sample cimk
      if(model == "two_stage"){
        if (threshold > 0) {

          w_all <- w[idx_w_k,,drop=F]

          # faster way to assign c_imk to 1 if logy1 > 0 for w_all = 1 and to 2
          # if log1 > 0 when w_all = 0
          c_imk <- y_pos * (2 - (w_all == 1))

        }
      }

      # sample theta
      if(model %in% c("occupancy","two_stage")){
        # beta_theta <- sample_betatheta_cpp(w, z, beta_theta, idx_z_w, X_theta,
        #                                    b_betatheta, B_betatheta)
        beta_theta <- sample_betatheta_cpp_parallel(w, z, beta_theta, idx_z_w, X_theta,
                                           b_betatheta, B_betatheta)
        theta <- computeTheta(X_theta, beta_theta)
      }

      # sample pq
      if(model == "two_stage"){
        list_pq <- sample_pq_cpp_parallel(c_imk, y_NA, w, idx_p_k, idx_w_k, maxP, a_p, b_p, a_q, b_q)
        p <- list_pq$p
        q <- list_pq$q
      }

      # sample theta0
      if(model %in% c("occupancy","two_stage")){
        theta0 <- sample_theta0(z, w, idx_z_w, a_theta0, b_theta0)
      }

      {

        if(iter > nburn & (iter - nburn) %% nthin == 0){
          currentIter <- (iter - nburn) / nthin

          if(model %in% c("occupancy","two_stage")){
            beta_theta_output_chain[,,currentIter] <- beta_theta
            theta0_output_chain[,currentIter] <- theta0
            if(summarisedLatentPresences){
              z_sum <- z_sum + z
              psi_sum <- psi_sum + psi
            } else {
              z_output_chain[,,currentIter] <- z
              psi_output_chain[,,currentIter] <- psi
            }
            theta_sum <- theta_sum + theta
          }

          if(model == "two_stage"){
            p_output_chain[,,currentIter] <- p
            q_output_chain[,,currentIter] <- q
            w_sum <- w_sum + w
          }

          # save jsdm params
          {
            B0_output_chain[,currentIter] <- list_jSDMparams$B0
            B_output_chain[,,currentIter] <- list_jSDMparams$B
            G_output_chain[,,currentIter] <- list_jSDMparams$G
            A_output_chain[,,currentIter] <- list_jSDMparams$A
            C_output_chain[,,currentIter] <- list_jSDMparams$C
            Bs_output_chain[,,currentIter] <- list_jSDMparams$Bs
            Gs_output_chain[,,currentIter] <- list_jSDMparams$Gs
            As_output_chain[,,currentIter] <- list_jSDMparams$As
            Cs_output_chain[,,currentIter] <- list_jSDMparams$Cs
            U_output_chain[,,currentIter] <- list_jSDMparams$U
            L_output_chain[,,currentIter] <- list_jSDMparams$L
            sigmab_output_chain[currentIter] <- list_jSDMparams$sigma_b
            sigmabs_output_chain[currentIter] <- list_jSDMparams$sigma_bs
            sigmah_output_chain[currentIter] <- list_jSDMparams$sigma_h
            idx_ls_output_chain[currentIter] <- list_jSDMparams$idx_ls
            if(model == "continuous") tau_output_chain[,currentIter] <- list_jSDMparams$tau

            varPart_output_chain[,,currentIter] <-
              as.matrix(list_jSDMparams$variancePartitioning)
          }

          # update WIAC
          {

            loglik_jsdm <- computeModelLoglikJSDM_cpp(z, eta, jsdmModel, tau)
            list_waic_jsdm <- update_waic_summary(loglik_jsdm, list_waic_jsdm, currentWAICiter)

            if(model %in% c("occupancy","two_stage")){
              logliks_firststage <- computeModelLoglikFirstStage_cpp(w, z, theta, theta0,
                                                                 list_idx$idx_z_w)
              list_waic_w <- update_waic_summary(logliks_firststage, list_waic_w, currentWAICiter)
            }

            if(model == "two_stage"){
              logliks_secondstage <- computeModelLoglikSecondStage_cpp(y, w, p, q,
                                                                   list_idx$idx_w_k,
                                                                   list_idx$idx_p_k)
              list_waic_y <- update_waic_summary(logliks_secondstage, list_waic_y, currentWAICiter)
            }

            currentWAICiter <- currentWAICiter + 1

          }

        }

      }

    }
  # ---- end of moved code ----

  if (!model %in% c("occupancy", "two_stage")) {
    beta_theta_output_chain <- NULL; theta0_output_chain <- NULL
    z_output_chain <- NULL; psi_output_chain <- NULL
  }
  if (model != "two_stage") { p_output_chain <- NULL; q_output_chain <- NULL }
  if (model %in% c("occupancy", "two_stage") && summarisedLatentPresences) {
    z_output_chain <- NULL; psi_output_chain <- NULL
  }

  # The number of WAIC accumulator updates that actually happened. The parent
  # checks it against niter for every part before merging by it.
  n_waic <- currentWAICiter - 1

  list(
    chain = chain, seed_r = seed_r, seed_cpp = seed_cpp,
    elapsed = proc.time()[["elapsed"]] - started,
    draws = list(
      beta_theta = beta_theta_output_chain, theta0 = theta0_output_chain,
      z = z_output_chain, psi = psi_output_chain,
      p = p_output_chain, q = q_output_chain,
      B0 = B0_output_chain, B = B_output_chain, G = G_output_chain, A = A_output_chain,
      C = C_output_chain, Bs = Bs_output_chain, Gs = Gs_output_chain, As = As_output_chain,
      Cs = Cs_output_chain, U = U_output_chain, L = L_output_chain,
      sigmab = sigmab_output_chain, sigmabs = sigmabs_output_chain, sigmah = sigmah_output_chain,
      idx_ls = idx_ls_output_chain, varPart = varPart_output_chain, tau = tau_output_chain),
    latent_sums = list(z = z_sum, psi = psi_sum, theta = theta_sum, w = w_sum),
    waic = list(jsdm = c(list_waic_jsdm, list(n = n_waic)),
                w = if (is.null(list_waic_w)) NULL else c(list_waic_w, list(n = n_waic)),
                y = if (is.null(list_waic_y)) NULL else c(list_waic_y, list(n = n_waic))))
}
