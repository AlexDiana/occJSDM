# Simulate the Quickstart's example survey: a two-stage eDNA river survey
# modelled on Ji et al. (2025, Ecology Letters), the Gaoligongshan study.
# 100 sites, 3 water samples per site, 2 primers, 3 PCR replicates per primer.
#
# simulateOccJSDMData() draws its own anonymous covariates and coefficients,
# so this script simulates the same model directly in order to plant named
# covariates, a two-level factor and chosen effects:
#   occupancy:  logit(psi[i, s]) = B0[s] + X[i, ] B[, s] + U[i, ] L[, s]
#   collection: w ~ Bernoulli(theta[i, s]) if z = 1, Bernoulli(theta0[s]) if z = 0
#   PCR:        y ~ Bernoulli(p[primer, s]) if w = 1, Bernoulli(q[primer, s]) if w = 0
# Read counts follow simulateOccJSDMData()'s defaults (R/simulateData.R).
#
# Run from the repository root: Rscript dev/release/2026-11-01-quickstart-data/simulate-quickstartdata.R
# Writes data/quickstartdata.rda and dev/release/2026-11-01-quickstart-data/truth.rds.

set.seed(20261101)
logistic <- function(x) 1 / (1 + exp(-x))
std <- function(x) (x - mean(x)) / sd(x)

n <- 100   # sites
M <- 3     # water samples per site
P <- 2     # primers
K <- 3     # PCR replicates per primer

# Species: vertebrates named Class_Order_number, so richness can be grouped.
groups <- c(Mammalia_Carnivora = 3, Mammalia_Artiodactyla = 3, Mammalia_Rodentia = 3,
            Mammalia_Eulipotyphla = 2, Aves_Passeriformes = 3, Amphibia_Anura = 3,
            Actinopterygii_Cypriniformes = 3)
speciesNames <- unlist(lapply(names(groups), function(g) paste0(g, "_", seq_len(groups[[g]]))))
S <- length(speciesNames)
class <- sub("_.*", "", speciesNames)

# Traits: body mass (kg) by group, and IUCN Red List status (0/1), more
# likely for large mammals.
log_mass_mean <- c(Mammalia_Carnivora = log(8), Mammalia_Artiodactyla = log(60),
                   Mammalia_Rodentia = log(0.3), Mammalia_Eulipotyphla = log(0.02),
                   Aves_Passeriformes = log(0.03), Amphibia_Anura = log(0.05),
                   Actinopterygii_Cypriniformes = log(0.2))
group_of <- sub("_[0-9]+$", "", speciesNames)
body_mass <- signif(exp(rnorm(S, log_mass_mean[group_of], 0.5)), 3)
red_listed <- c(1, 0, 1,  1, 1, 0,  0, 0, 0,  0, 0,  0, 0, 1,  1, 0, 0,  0, 0, 0)
traits <- cbind(body_mass = body_mass, red_listed = red_listed)
rownames(traits) <- speciesNames

# Sites: elevation, Human Modification Index and protection status. Protected
# streams tend to be higher and less modified, but not perfectly so.
protected01 <- rbinom(n, 1, 0.45)
elevation <- round(pmin(pmax(rnorm(n, 1800 + 600 * protected01, 900), 200), 4800))
hmi <- round(logistic(-0.3 - 0.8 * protected01 - 0.2 * std(elevation) + rnorm(n, 0, 0.7)), 3)
protected <- factor(ifelse(protected01 == 1, "protected", "not_protected"),
                    levels = c("not_protected", "protected"))

# Occupancy coefficients on the scale runOccJSDM() fits: numeric covariates
# standardised, protected as a 0/1 indicator.
X <- cbind(elevation = std(elevation), hmi = std(hmi), protected = protected01)
z_mass <- std(log(body_mass))
B0 <- rnorm(S, -0.4, 0.9)
B <- rbind(
  # Fish and frogs favour low elevations; the others vary.
  elevation = ifelse(class %in% c("Actinopterygii", "Amphibia"), -1.0, 0.2) + rnorm(S, 0, 0.4),
  # Planted trait-by-environment effect: larger species decline more steeply
  # with human modification.
  hmi = -0.4 - 0.7 * z_mass + rnorm(S, 0, 0.25),
  # Planted trait-by-protected effect: Red Listed species favour protected sites.
  protected = 0.1 + 1.4 * red_listed + rnorm(S, 0, 0.3)
)
d <- 2
U <- matrix(rnorm(n * d), n, d)
L <- matrix(rnorm(d * S, 0, 0.7), d, S)
psi <- logistic(matrix(B0, n, S, byrow = TRUE) + X %*% B + U %*% L)
z <- matrix(rbinom(n * S, 1, psi), n, S)

# Field collection: 3 samples per site, water volume filtered (litres) raises
# collection probability. Two species carry noticeable field contamination.
N <- n * M
site_of_sample <- rep(seq_len(n), each = M)
water_volume <- round(exp(rnorm(N, log(2), 0.4)), 2)
beta_theta <- rbind(intercept = rnorm(S, qlogis(0.55), 0.5),
                    water_volume = 0.6 + rnorm(S, 0, 0.15))
theta <- logistic(cbind(1, std(water_volume)) %*% beta_theta)
theta0 <- rep(0.01, S)
contaminated <- c("Mammalia_Rodentia_2", "Aves_Passeriformes_2")
theta0[match(contaminated, speciesNames)] <- 0.12
z_sample <- z[site_of_sample, ]
w <- matrix(rbinom(N * S, 1, ifelse(z_sample == 1, theta, matrix(theta0, N, S, byrow = TRUE))), N, S)

# PCR: primer 1 is better for mammals and birds, primer 2 for fish and frogs,
# with species-level variation, so neither primer is best for every species.
primer1_better <- class %in% c("Mammalia", "Aves")
p <- rbind(primer1 = logistic(qlogis(ifelse(primer1_better, 0.8, 0.45)) + rnorm(S, 0, 0.4)),
           primer2 = logistic(qlogis(ifelse(primer1_better, 0.5, 0.8)) + rnorm(S, 0, 0.4)))
q <- matrix(0.01, P, S)

# One row per PCR: by sample, then primer, then replicate.
sample_of_pcr <- rep(seq_len(N), each = P * K)
primer_of_pcr <- rep(rep(seq_len(P), each = K), times = N)
w_pcr <- w[sample_of_pcr, ]
p_pcr <- p[primer_of_pcr, ]
q_pcr <- q[primer_of_pcr, ]
state <- matrix(0L, nrow(w_pcr), S)  # 0 no reads, 1 true detection, 2 false positive
state[w_pcr == 1] <- rbinom(sum(w_pcr == 1), 1, p_pcr[w_pcr == 1])
state[w_pcr == 0] <- 2L * rbinom(sum(w_pcr == 0), 1, q_pcr[w_pcr == 0])
logy <- matrix(0, nrow(state), S)
logy[state == 1] <- rnorm(sum(state == 1), 5, 1)
logy[state == 2] <- rnorm(sum(state == 2), 1.5, 1)
OTU <- pmax(round(exp(logy) - 1), 0)
colnames(OTU) <- speciesNames

info <- data.frame(
  Site = site_of_sample[sample_of_pcr],
  Sample = sample_of_pcr,
  Primer = primer_of_pcr,
  elevation = elevation[site_of_sample[sample_of_pcr]],
  hmi = hmi[site_of_sample[sample_of_pcr]],
  protected = protected[site_of_sample[sample_of_pcr]],
  water_volume = water_volume[sample_of_pcr]
)

quickstartdata <- list(info = info, OTU = OTU, traits = traits)
save(quickstartdata, file = "data/quickstartdata.rda", compress = "xz")

dimnames(z) <- dimnames(psi) <- list(as.character(seq_len(n)), speciesNames)
colnames(B) <- colnames(L) <- names(B0) <- names(theta0) <- colnames(beta_theta) <- colnames(p) <- speciesNames
saveRDS(list(B0 = B0, B = B, U = U, L = L, psi = psi, z = z, beta_theta = beta_theta,
             theta0 = theta0, p = p, q = q, w = w, contaminated = contaminated),
        "dev/release/2026-11-01-quickstart-data/truth.rds")
