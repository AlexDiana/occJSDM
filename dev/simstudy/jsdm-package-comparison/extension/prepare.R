args <- commandArgs(trailingOnly = TRUE)
code <- normalizePath(args[1]); root <- args[2]
source(file.path(code, "study.R"))
stopifnot(!file.exists(file.path(root, "manifest.csv")))
for (d in c("inputs", "truth", "jobs", "logs")) dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
rows <- list()
for (replicate in 1:10) for (scenario in c("baseline", "rare", "correlated", "curved", "traits")) {
  community <- simulate_community(replicate, scenario)
  dataset <- sprintf("%s-r%02d", scenario, replicate)
  saveRDS(community, file.path(root, "truth", paste0(dataset, ".rds")))
  species_counts <- if (scenario == "traits") c(10, 30) else 10
  responses <- if (scenario == "curved") c("linear", "quadratic") else "linear"
  trait_arms <- if (scenario == "traits") c(FALSE, TRUE) else FALSE
  packages <- if (scenario == "traits") c("occJSDM", "gllvm", "Hmsc") else c("occJSDM", "gllvm", "sjSDM", "Hmsc")
  for (n_sites in c(100, 300)) for (n_species in species_counts) for (response in responses) for (use_traits in trait_arms) {
    id <- sprintf("%s-n%d-s%d-%s-traits%d", dataset, n_sites, n_species, response, use_traits)
    input <- make_input(community, n_sites, n_species, response, use_traits)
    path <- file.path(root, "inputs", paste0(id, ".rds"))
    saveRDS(input, path)
    for (package in packages) rows[[length(rows) + 1L]] <- data.frame(
      job = paste(id, package, sep = "-"), dataset, input = id,
      replicate, scenario, n_sites, n_species, response, use_traits, package,
      input_md5 = unname(tools::md5sum(path)))
  }
}
manifest <- do.call(rbind, rows)
manifest$job_number <- seq_len(nrow(manifest))
write.csv(manifest, file.path(root, "manifest.csv"), row.names = FALSE)
cat("Prepared", nrow(manifest), "package/model combinations, before any fitting.\n")
print(with(manifest, table(scenario, package)))
