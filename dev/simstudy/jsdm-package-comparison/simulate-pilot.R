# New independent community for Lesson 4. All coefficients are fixed before fitting.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1])
stopifnot(!file.exists(file.path(root, "inputs", "training.rds")))

parameters <- data.frame(
  species = sprintf("species_%02d", 1:10),
  intercept = c(-1.8, -1.3, -0.9, -0.5, -0.2, 0.2, 0.5, 0.9, 1.3, 1.8),
  environment_1 = c(1.2, 0.8, -1.0, -0.7, 1.0, -1.1, 0.6, -0.8, 1.0, -1.2),
  environment_2 = c(-0.6, 0.9, 0.7, -1.0, 0.4, 0.8, -0.9, 0.6, -0.7, 0.5),
  hidden_1 = c(0.9, 0.7, -0.8, -0.6, 0.1, 0.2, 0.8, -0.9, 0.5, -0.6),
  hidden_2 = c(0.2, -0.3, 0.2, -0.4, 0.9, -0.8, 0.6, -0.5, -0.8, 0.7)
)
set.seed(26092201)

raw_x <- matrix(rnorm(400 * 2), nrow = 400, ncol = 2,
                dimnames = list(c(sprintf("train_%03d", 1:100),
                                  sprintf("test_%03d", 1:300)),
                                c("environment_1", "environment_2")))
hidden <- matrix(rnorm(400 * 2), nrow = 400, ncol = 2)
coefficients <- t(as.matrix(parameters[, 2:4]))
loadings <- t(as.matrix(parameters[, 5:6]))
linear_predictor <- cbind(1, raw_x) %*% coefficients + hidden %*% loadings
probability <- plogis(linear_predictor)
colnames(probability) <- parameters$species
occurrence <- matrix(rbinom(length(probability), size = 1, prob = probability),
                     nrow = 400, dimnames = dimnames(probability))

training_rows <- 1:100
test_rows <- 101:400
centre <- colMeans(raw_x[training_rows, ])
spread <- apply(raw_x[training_rows, ], 2, sd)
x <- sweep(sweep(raw_x, 2, centre, "-"), 2, spread, "/")
scaled_coefficients <- coefficients
scaled_coefficients[1, ] <- coefficients[1, ] + centre %*% coefficients[-1, ]
scaled_coefficients[-1, ] <- coefficients[-1, ] * spread
stopifnot(max(abs(cbind(1, x) %*% scaled_coefficients -
                  cbind(1, raw_x) %*% coefficients)) < 1e-12)

training <- list(x = as.data.frame(x[training_rows, ]),
                 y = occurrence[training_rows, ], centre = centre, spread = spread)
saveRDS(training, file.path(root, "inputs", "training.rds"))
saveRDS(as.data.frame(x[test_rows, ]), file.path(root, "inputs", "test-x.rds"))
saveRDS(list(parameters = parameters, raw_x = raw_x, hidden = hidden,
             conditional_probability = probability, occurrence = occurrence,
             scaled_coefficients = scaled_coefficients, loadings = loadings,
             seed = 26092201, training_rows = training_rows, test_rows = test_rows),
        file.path(root, "truth", "truth.rds"))
write.csv(parameters, file.path(root, "truth", "generating-parameters.csv"), row.names = FALSE)
write.table(tools::md5sum(list.files(file.path(root, "inputs"), full.names = TRUE)),
            file.path(root, "inputs", "hashes.tsv"), sep = "\t", col.names = FALSE)
cat("Saved one fixed community. Training presences per species:\n")
print(colSums(training$y))
