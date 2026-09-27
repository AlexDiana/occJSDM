# Independently check exported identities, probabilities, and example truth.
# Rscript dev/simstudy/vignette-lesson/verify_latent_tables.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
source("dev/simstudy/vignette-lesson/helpers.R")

lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
lesson <- readRDS(lesson_path)
bundle <- readRDS("vignettes/teaching-data/latent-presence-lesson.rds")
stopifnot(identical(bundle$lesson_md5, unname(tools::md5sum(lesson_path))),
          identical(bundle$source_hashes, lesson_source_hashes()),
          identical(bundle$exporter_md5, unname(tools::md5sum(
            "dev/simstudy/vignette-lesson/summarise_latent_tables.R"))))
fit_path <- file.path(archive, bundle$fit_manifest$file)
stopifnot(identical(unname(tools::md5sum(fit_path)), bundle$fit_manifest$md5))
fit <- readRDS(fit_path)$fit
x <- bundle$tables
stopifnot(nrow(x) == 24000L)

sample_ids <- unique(fit$infos$data_info[c("Site", "Sample")])
site_index <- match(as.character(x$Site), as.character(fit$infos$siteNames))
sample_key <- function(x) paste(x$Site, x$Sample, sep = "/")
sample_index <- match(sample_key(x), sample_key(sample_ids))
primer_index <- match(as.character(x$Primer), as.character(fit$infos$primerNames))
species_index <- match(x$species, fit$infos$speciesNames)
stopifnot(!anyNA(c(site_index, sample_index, primer_index, species_index)))

equal <- function(a, b) stopifnot(isTRUE(all.equal(unname(a), unname(b),
                                                  tolerance = 1e-12)))
out <- fit$results_output
equal(x$CondOccProb, out$z_output[cbind(site_index, species_index)])
equal(x$CondSampleProb, out$w_output[cbind(sample_index, species_index)])
equal(x$PredOccProb, out$psi_output[cbind(site_index, species_index)])
equal(x$CollectionProb, out$theta_output[cbind(sample_index, species_index)])
p_mean <- apply(out$p_output, c(1, 2), mean)
equal(x$DetectionProb, p_mean[cbind(primer_index, species_index)])

# Look up observed reads directly from the original input. Do not use the
# exporter's join or rely on its row order.
info <- lesson$input$sim$data_list$info
pcr <- ave(seq_len(nrow(info)), interaction(info$Site, info$Sample, info$Primer),
           FUN = seq_along)
observation_key <- function(site, sample, primer, pcr) {
  paste(site, sample, primer, pcr, sep = "/")
}
row_index <- match(observation_key(x$Site, x$Sample, x$Primer, x$PCR),
                   observation_key(info$Site, info$Sample, info$Primer, pcr))
equal(x$OTU, lesson$input$sim$data_list$OTU[cbind(row_index, species_index)])

# The student joins on identities, so permuting either table must not change
# which known state/read is attached to an exported record.
keys <- c("species", "Site", "Sample", "Primer", "PCR")
observed_truth <- lesson$observations[, c(keys, "reads", "z", "w", "source")]
set.seed(20260922)
joined <- dplyr::left_join(x[sample.int(nrow(x)), ],
                          observed_truth[sample.int(nrow(observed_truth)), ],
                          by = keys, relationship = "one-to-one")
equal(joined$OTU, joined$reads)
truth <- lesson$input$sim$true_params
zi <- cbind(match(as.character(joined$Site), rownames(truth$z_true)),
            match(joined$species, colnames(truth$z_true)))
wi <- cbind(match(as.character(joined$Sample), rownames(truth$w_true)),
            match(joined$species, colnames(truth$w_true)))
equal(joined$z, truth$z_true[zi])
equal(joined$w, truth$w_true[wi])

duplicated_truth <- rbind(observed_truth, observed_truth[1, ])
duplicate_failure <- tryCatch({
  dplyr::left_join(x, duplicated_truth, by = keys, relationship = "one-to-one")
  FALSE
}, error = function(e) TRUE)
stopifnot(duplicate_failure)

example <- joined[joined$species == "OTU_1" & joined$Site %in% c(2, 3), ]
stopifnot(nrow(example) == 48L, all(example$z == 1),
          any(example$source == "Laboratory false positive"),
          any(example$source == "True detection"))
tbl <- occJSDM::plotLatentPresences(example, species_name = "OTU_1")
stopifnot(inherits(tbl, "gt_tbl"), nrow(tbl$`_data`) == 48L)
cat("Verified all 24,000 native records against archived probabilities and input reads.\n")
cat("Truth joins survive row permutations and reject duplicate identities; native table builds.\n")
