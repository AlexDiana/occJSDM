# The implementation of computeDiagnostics() shipped up to occJSDM 0.3.0,
# kept verbatim as a test oracle.
#
# It took R-hat and effective sample size for a whole parameter block from
# coda::gelman.diag() and coda::effectiveSize(), which is correct but forms
# parameters-by-parameters matrices, so its cost grew with the square of the
# block size. The replacement computes the same statistics parameter by
# parameter. Its contract is that everything a user sees (the printed summary,
# the warnings, the skipped blocks) is unchanged, and the direct way to assert
# that is to run both on the same input and compare the console output.
#
# Do not edit this function to make a test pass: it is the reference. Delete
# it, and the tests that use it, only if the printed summary is changed on
# purpose.
legacy_computeDiagnostics <- function(results_output){

  nchain <- dim(results_output$jsdm_output$B0_output)[3]
  niter <- dim(results_output$jsdm_output$B0_output)[2]

  results_output_all <- c(
    results_output[names(results_output) != "jsdm_output"],
    results_output$jsdm_output)

  for (param_name in names(results_output_all)) {

    if(!(param_name %in% c("z_output","w_output","psi_output",
                           "theta_output","idx_ls_output","varPart_output"))){

      arr <- results_output_all[[param_name]]
      dims <- dim(arr)

      if (is.null(dims) || length(dims) < 2 || any(dims== 0 ) || any(is.na(arr)) ) {
        next
      }

      ndim <- length(dims)

      if (ndim == 2) {
        n_params <- 1
      } else {
        n_params <- prod(dims[1:(ndim - 2)])
      }

      # Reshape uniformly to 3D: (n_params, niter, nchain)
      arr_3d <- array(arr, dim = c(n_params, niter, nchain))

      arr_flat <- matrix(arr_3d, nrow = n_params, ncol = niter * nchain)
      param_vars <- apply(arr_flat, 1, var, na.rm = TRUE)

      # Keep only parameters with variance greater than a tiny threshold
      is_variable <- !is.na(param_vars) & param_vars > 1e-12
      valid_params <- sum(is_variable)

      if (valid_params == 0) {
        message(sprintf("\nParameter Set: %s", param_name))
        message("--------------------------------------------------")
        message("  [Skipped: All parameters in this set are fixed/constant]")
        next
      }

      # Subset to remove the fixed parameters
      arr_3d <- arr_3d[is_variable, , , drop = FALSE]

      # Convert to coda mcmc.list
      chain_list <- lapply(1:nchain, function(c) {
        # Use matrix() to prevent R from dropping dimensions on single parameters
        chain_slice <- matrix(arr_3d[, , c], nrow = valid_params, ncol = niter)
        coda::mcmc(t(chain_slice)) # Transpose so rows = iterations, cols = parameters
      })

      mcmc_obj <- coda::mcmc.list(chain_list)

      # Diagnostics
      ess <- coda::effectiveSize(mcmc_obj)

      if (nchain > 1) {
        rhat <- tryCatch({
          diag_res <- coda::gelman.diag(mcmc_obj, autoburnin = FALSE, multivariate = FALSE)
          diag_res$psrf[, 1]
        }, error = function(e) rep(NA, n_params))
      } else {
        rhat <- rep(NA, n_params)
      }

      # Print Summary
      param_name_to_print <- gsub("_output","",param_name)
      message(sprintf("\nParameter Set: %s (Flattened N = %d)", param_name_to_print, n_params))
      message("--------------------------------------------------")

      # Use suppressWarnings for min/max to prevent warnings if all values are NA (e.g., fixed params)
      max_rhat <- suppressWarnings(max(rhat, na.rm = TRUE))
      mean_rhat <- mean(rhat, na.rm = TRUE)
      min_ess <- suppressWarnings(min(ess, na.rm = TRUE))
      mean_ess <- mean(ess, na.rm = TRUE)

      if (nchain > 1) {
        max_rhat <- suppressWarnings(max(rhat, na.rm = TRUE))
        mean_rhat <- mean(rhat, na.rm = TRUE)
        message(sprintf("  R-hat | Max:  %7.3f | Mean: %7.3f", max_rhat, mean_rhat))

        if (is.finite(max_rhat) && max_rhat > 1.1) {
          warning(sprintf("Convergence issue in '%s'. Max R-hat > 1.1", param_name), call. = FALSE)
        }
      } else {
        message("  R-hat | [Requires multiple chains to compute]")
      }

      message(sprintf("  ESS   | Min:  %7.1f | Mean: %7.1f", min_ess, mean_ess))

      if (is.finite(min_ess) && min_ess < 50) {
        warning(sprintf("Low ESS in '%s'. Min ESS < 50", param_name_to_print), call. = FALSE)
      }

    }

  }

  message("\n==================================================")

  invisible(NULL)

}
