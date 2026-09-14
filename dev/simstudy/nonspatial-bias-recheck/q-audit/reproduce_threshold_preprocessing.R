# Execute the installed entry point's actual preprocessing, then return
# immediately before any covariate construction or MCMC. Package unchanged.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
.libPaths(c(file.path(root, "library"), .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
capture_preprocessing <- get("runOccJSDM", asNamespace("occJSDM"))
original_body <- body(capture_preprocessing)
parts <- as.list(original_body)
block <- which(vapply(parts, function(x) any(grepl("y <- OTU", deparse(x), fixed = TRUE)), logical(1)))
stopifnot(length(block) == 1L)
body(capture_preprocessing) <- as.call(c(parts[seq_len(block)],
  list(quote(return(list(original = OTU, actual = y, model = model))))))
dat <- list(info = data.frame(Site = rep(1:2, each = 4),
                             Sample = rep(1:4, each = 2), Primer = 1),
            OTU = cbind(sp1 = c(0,1,2,3,4,5,NA,100),
                        sp2 = c(100,NA,5,4,3,2,1,0)))
rows <- lapply(c(1,2,3), function(threshold) {
  result <- suppressMessages(capture_preprocessing(dat, threshold = threshold))
  stopifnot(result$model == "two_stage")
  expected <- (result$original >= threshold) * 1
  data.frame(threshold = threshold, species = rep(colnames(dat$OTU), each = nrow(dat$OTU)),
             row = rep(seq_len(nrow(dat$OTU)), ncol(dat$OTU)),
             original = as.vector(result$original), actual = as.vector(result$actual),
             expected = as.vector(expected))
})
rows <- do.call(rbind, rows)
stopifnot(all(rows$actual[rows$threshold == 1] == rows$expected[rows$threshold == 1], na.rm = TRUE),
          all(rows$actual[rows$threshold > 1] == 0, na.rm = TRUE),
          any(rows$actual[rows$threshold > 1] != rows$expected[rows$threshold > 1], na.rm = TRUE),
          identical(body(get("runOccJSDM", asNamespace("occJSDM"))), original_body))
write.csv(rows, file.path(root, "q-audit/threshold-preprocessing-reproduction.csv"), row.names = FALSE)
print(rows)
cat("Confirmed in the installed entry point: threshold one is correct; thresholds two and three erase every observed detection. No MCMC or package modifications.\n")
