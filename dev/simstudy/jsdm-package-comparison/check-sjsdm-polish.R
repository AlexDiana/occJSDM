args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
source(file.path(code, "sjsdm-gradient.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
evaluate <- make_joint_gradient(input$x, input$y, 121L)
rows <- predictions <- list()
for (i in 1:3) {
  label <- sprintf("sjSDM-default-penalty-polish-start-%d", i)
  base <- file.path(root, "stability/fits", label)
  f <- readRDS(paste0(base, ".rds"))
  r <- readRDS(paste0(base, "-record.rds"))
  stopifnot(isTRUE(r$ok), identical(f$input,input), length(f$history)==1000,
            r$settings$weight_decay==.0001, r$settings$seed==26092640L+i)
  p <- f$parameters
  theta <- c(p$beta,p$loading)
  stopifnot(max(abs((cbind(1,as.matrix(input$x))%*%p$beta)*.999999+.0000005-f$native_raw))<1e-10)
  a <- evaluate(theta)
  b <- joint_integration(p,input$x,input$y,241L,clamp=TRUE)
  stopifnot(abs(a$value+b$loglik)<.001)
  penalty <- 100*.0001*sum(theta^2)/2
  gradient <- a$gradient+100*.0001*theta
  predictions[[i]] <- point_marginal(p,fixed_grid())
  trajectory <- lapply(c(0,250,500,750,1000),function(epoch){
    checkpoint <- readRDS(sprintf("%s-epoch-%04d.rds",base,epoch))
    t <- c(checkpoint$parameters$beta,checkpoint$parameters$loading)
    data.frame(start=i,epoch=epoch,penalised_score=-checkpoint$diagnostic$value-100*.0001*sum(t^2)/2)
  })
  trajectory <- do.call(rbind,trajectory)
  rows[[i]] <- data.frame(start=i,loglik=-a$value,penalty=penalty,
    penalised_score=-a$value-penalty,max_penalised_gradient=max(abs(gradient)),
    difference_121_241=abs(a$value+b$loglik),
    late_improvement=tail(trajectory$penalised_score,1)-tail(trajectory$penalised_score,2)[1],
    maximum_residual_sd=max(sqrt(colSums(p$loading^2))))
  write.csv(trajectory,file.path(root,"stability/checks",paste0(label,"-trajectory.csv")),row.names=FALSE)
  saveRDS(list(parameters=p,gradient=gradient,refined=b),
          file.path(root,"stability/checks",paste0(label,".rds")))
}
rows <- do.call(rbind,rows)
spread <- max(apply(simplify2array(predictions),1:2,function(p) diff(range(p))))
summary <- list(rows=rows,objective_spread=diff(range(rows$penalised_score)),
  prediction_spread_pp=100*spread,passed=diff(range(rows$penalised_score))<=.1&&spread<=.01)
saveRDS(summary,file.path(root,"stability/checks/penalty-polish-summary.rds"))
write.csv(rows,file.path(root,"stability/checks/penalty-polish-summary.csv"),row.names=FALSE)
print(summary)
