# Projection of the stage 2 run time from the Accelerate per-iteration costs in results/probe-accelerate.txt (SHELVED.md 9.3).
# Pricing as the plan and SHELVED.md 4.1: four chains of 7,000 iterations per job, 20 jobs at 900 sites and 20 at 300 per arm,
# the 300-site jobs at their measured 300-site cost, 4.43 times for the all-extension upper bound, gllvm excluded, 8 workers.
# env is priced at its reference-BLAS cost from SHELVED.md 4.2 (not re-measured). Usage:  Rscript probe-accelerate-projection.R
s <- c(gp900 = 1.461, gp600 = 0.667, gp450 = 0.470, gp300 = 0.273, gp_n300 = 0.140, h2_900 = 0.307, h2_300 = 0.022, h6_900 = 0.545, h6_300 = 0.043, env = 0.0075)
h <- unname(s * 28000 / 3600); names(h) <- names(s)
arm <- c(gp900 = 20*h[["gp900"]] + 20*h[["gp_n300"]], gp600 = 20*h[["gp600"]] + 20*h[["gp_n300"]], gp450 = 20*h[["gp450"]] + 20*h[["gp_n300"]],
         gp300 = 20*h[["gp300"]] + 20*h[["gp_n300"]], h2 = 20*h[["h2_900"]] + 20*h[["h2_300"]], h6 = 20*h[["h6_900"]] + 20*h[["h6_300"]], env = 40*h[["env"]])
cat("worker-hours per arm (40 jobs), initial schedule:\n"); print(round(arm, 1))
fixed <- arm[["h2"]] + arm[["h6"]] + arm[["env"]]
for (g in c("gp900", "gp600", "gp450", "gp300")) { tot <- arm[[g]] + fixed
  cat(sprintf("ladder %s: initial %.0f worker-hours = %.0f wall-hours at 8 workers; all-extension x4.43 = %.0f worker-hours = %.0f wall-hours; one third extend = %.0f wall-hours\n",
      g, tot, tot/8, tot*4.43, tot*4.43/8, tot*2.13/8)) }
cat("speed-ups vs SHELVED.md: gp300", round(7.9/0.273, 1), " gp900 (contended baseline)", round(89.3/1.461, 1), " Hmsc2", round(1.65/0.307, 1), " Hmsc6", round(26/0.545, 1), "\n")
