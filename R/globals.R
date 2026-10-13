# Columns evaluated inside ggplot2 aesthetics, dplyr data masks or tidyselect.
# Each name and its source functions were reviewed after dead-code cleanup;
# see dev/maintenance/post-beta-maintenance-20261007/data-masked-columns.csv.
# Undefined symbols in retained experimental helpers are intentionally excluded.
utils::globalVariables(c(
  "2.5%",
  "97.5%",
  "Biotic",
  "Env",
  "Max",
  "Min",
  "Output",
  "Primer",
  "Sample",
  "Site",
  "Spatial",
  "Species",
  "chain",
  "high",
  "iter",
  "label",
  "label1",
  "label2",
  "low",
  "lower",
  "mean_lower",
  "med",
  "name",
  "p1",
  "p2",
  "primer_group",
  "primer_parity",
  "q1",
  "q2",
  "r",
  "rownum",
  "sample_group",
  "sample_parity",
  "shade_level",
  "site_color_idx",
  "site_group",
  "species",
  "upper",
  "value",
  "x",
  "y"
))

# Fields of the setup list that run_one_chain() in R/chains.R unpacks into its
# frame with list2env(); the moved chain body reads them as plain names.
# Keep in step with chain_setup_fields().
utils::globalVariables(c(
  "model", "jsdmModel", "threshold",
  "y", "y_NA", "z", "w", "M", "K", "P", "sumM", "sumP", "sumK", "maxP", "primerId_p",
  "idx_w_k", "idx_p_k", "idx_z_w", "list_idx",
  "n", "N", "N3", "S", "d", "g", "gt", "ps", "ncov_theta", "ncov_psi",
  "X_theta", "X_psi", "Tr",
  "list_Xs", "list_SoRSummaries",
  "list_data", "list_priors", "spatial_sd_prior",
  "a_theta0", "b_theta0", "a_p", "b_p", "a_q", "b_q", "b_betatheta", "B_betatheta",
  "nburn", "niter", "nthin", "summarisedLatentPresences"
))
