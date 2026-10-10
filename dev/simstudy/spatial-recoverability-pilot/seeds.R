# Seed namespace for the spatial recoverability pilot. Labels are hashed by a
# position-weighted byte sum so that every random stream is reproducible and
# distinct from the sweep's label*1000+replicate scheme.
pilot_label <- function(...) paste(c(...), collapse = "-")

pilot_seed <- function(label) {
  stopifnot(is.character(label), length(label) == 1L, nzchar(label))
  b <- as.numeric(utf8ToInt(paste0("recoverability-pilot:", label)))
  as.integer((sum(b * seq_along(b)) * 7919 + sum(b) + 1) %% 2147483646 + 1)
}
