args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
out <- file.path(code, "stability-results")
dir.create(out, showWarnings = FALSE)
profiles <- read.csv(file.path(root, "stability/checks/scale-profiles.csv"))
p <- ggplot(profiles, aes(scale, -improvement, colour = factor(species))) +
  geom_hline(yintercept = 0, colour = "grey65") + geom_line(linewidth = .8) +
  geom_point(size = 2) +
  scale_x_log10(breaks = c(.5,1,2,4,8,16), labels = c("half","1x","2x","4x","8x","16x")) +
  scale_colour_manual(values = c("#007A87", "#CB6E17"), labels = c("Species 2", "Species 7")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x = "Multiplier on one species' intercept, environmental slopes and factor loadings",
       y = "Loss of fit relative to the original estimate\n(total log-likelihood units)", colour = NULL,
       title = "Very different effect sizes can fit nearly equally well",
       subtitle = "Other species' parameters are held fixed; no true values are used",
       caption = "Zero is the original fit. Positive values mean a worse fit.\nSpecies 2 changes by only 0.084 at 16x; this is weak scale information, not proof of an infinite optimum.") +
  theme(plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "effect-scale-profile.png"), p, width = 10, height = 5.7, dpi = 160)

parameters <- readRDS(file.path(root, "stability/fits/sjSDM-low-rate-start-1.rds"))$parameters
x <- data.frame(environment_1 = seq(-2,2,length.out=81), environment_2 = 0)
curves <- lapply(c(.5,1,4,16), function(scale) {
  changed <- parameters
  changed$beta[,2] <- changed$beta[,2]*scale
  changed$loading[,2] <- changed$loading[,2]*scale
  pred <- point_marginal(changed,x)[,2]
  data.frame(environment=x$environment_1, probability=pred, scale=factor(scale,levels=c(.5,1,4,16)))
}) |> bind_rows()
p <- ggplot(curves, aes(environment, probability, colour=scale, linetype=scale)) +
  geom_line(linewidth=.9) +
  scale_y_continuous(labels=scales::label_percent(), limits=c(0,1)) +
  scale_colour_manual(values=c("#777777","#007A87","#CB6E17","#705A93"),
                      labels=c("Half","Original","4x effects","16x effects")) +
  scale_linetype_manual(values=c("dotted","solid","longdash","twodash"),
                       labels=c("Half","Original","4x effects","16x effects")) +
  theme_bw(base_size=12) + theme(legend.position="bottom",plot.caption=element_text(hjust=0)) +
  labs(x="Environmental predictor 1 (training standard deviations)",
       y="Predicted probability at a new site",colour=NULL,linetype=NULL,
       title="Large changes in coefficients need not mean large changes in predictions",
       subtitle="Species 2; environmental predictor 2 held at its training mean",
       caption="Environmental effects and hidden-factor effects are scaled together.\nThese are diagnostic fitted curves, not comparisons with simulated truth.")
ggsave(file.path(out,"effect-scale-predictions.png"),p,width=10,height=5.7,dpi=160)
write.csv(curves,file.path(out,"effect-scale-predictions.csv"),row.names=FALSE)
for (f in c("initial-gradients.csv","low-rate-summary.csv","scale-profiles.csv")) {
  file.copy(file.path(root,"stability/checks",f),out,overwrite=TRUE)
}
refs <- lapply(1:3,function(i){
  r<-readRDS(file.path(root,"stability/fits",sprintf("reference-start-%d.rds",i)))$summary
  data.frame(start=i,convergence=r$convergence,loglik_81=r$loglik_81,
    loglik_161=r$loglik_161,loglik_241=r$loglik_241,
    maximum_residual_sd=r$maximum_residual_sd,
    accepted=FALSE,reason="Likelihood fails finer-grid verification")
}) |> bind_rows()
write.csv(refs,file.path(out,"rejected-reference-checks.csv"),row.names=FALSE)
