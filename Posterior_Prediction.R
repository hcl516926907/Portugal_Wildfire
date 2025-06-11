dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire'


library(INLA)
library(doParallel)

print(detectCores())


inla.setOption(scale.model.default = TRUE)

load(file.path(dir.out, 'dataset_perf_evaluate.RData'))



# load(file=file.path(dir.out, 'Model_egp_bym2.RData'))
# load(file=file.path(dir.out, 'Model_weibull_bym2.RData'))
# load(file=file.path(dir.out, 'Model_gamma_bym2.RData'))

load(file=file.path(dir.out, 'Model_gamma_bym2_h1.RData'))
# load(file=file.path(dir.out, 'Model_egp_bym2_h1.RData'))
# load(file=file.path(dir.out, 'Model_weibull_bym2_h1.RData'))
# load(file=file.path(dir.out, 'Model_egp_bym2_h1_noxgb.RData'))

result <- res
print(summary(result))



####################################################################
# Generate posterior samples of effects and hyperparameters
####################################################################

n.samples <- 1000
samples = inla.posterior.sample(n.samples, result = result, seed=1234)
apply(sapply(samples, function(x) x$hyperpar),1,summary)


n1 <- nrow(data.fit)

####################################################################
# Define the sample function of posterior predictive distribution on
# three eGP, Gamma and Weibull likelhoods
####################################################################

# weillbull  
post.pred.weibull.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
  pred.cnt.mat <- c()
  pred.z.mat <- c()
  pred.ba.mat <- c()
  param.cnt <- list()
  param.z <- list()
  param.ba <- list()
  hyper.ba <- list()

  for (i in 1:n.samples){
    eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
    lambda.pois <- exp(eta.pois)


    eta.z <- samples[[i]]$latent[idx.pred.z,1]
    p <- exp(eta.z)/(1 + exp(eta.z))


    eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
    alpha.ba <- samples[[i]]$hyperpar[1]

    lambda.ba  <-  exp(eta.ba)
    scale.ba <- lambda.ba^(-1/alpha.ba)

    pred.z <- rbinom(length(idx.pred.pois), size=1, prob=p)
    pred.cnt <- rpois(length(idx.pred.z), lambda.pois )
    pred.ba <- rweibull(length(idx.pred.ba), shape = alpha.ba, scale = scale.ba)

    zero.idx <- which(pred.z==0)
    pred.cnt[zero.idx] <- 0
    pred.ba[zero.idx] <- 0


    pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
    pred.z.mat <- cbind(pred.z.mat, pred.z)
    pred.ba.mat <- cbind(pred.ba.mat, pred.ba)

    param.cnt[[i]] <- lambda.pois
    param.z[[i]] <- p
    param.ba[[i]] <- scale.ba
    hyper.ba[[i]] <- alpha.ba
  }


  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:length(idx.pred.pois)){
    pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
    pred.z[[i]] <- as.numeric(pred.z.mat[i,])
    pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
  }

  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z
              # 'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 
              # 'hyper.ba'=hyper.ba
              ))
}

#egpd
post.pred.egp.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
  pred.cnt.mat <- c()
  pred.z.mat <- c()
  pred.ba.mat <- c()
  param.cnt <- list()
  param.z <- list()
  param.ba <- list()
  hyper.xi <- list()
  hyper.kappa <- list()
  
  regp <- function(n, eta, kappa, xi,alpha=0.5){
    a <- ((1-alpha^(1/kappa))^(-xi) -1)
    # x <- rnorm(n, sd = 0.3)
    # eta <- 0.9 + 1.1 * x
    q <- exp(eta)
    sigma <- xi * q / a
    y <- - (sigma / xi) * (1- (1-runif(n)^(1/kappa))^(-xi))
    return(y)
  }
  
  
  for (i in 1:n.samples){
    eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
    lambda.pois <- exp(eta.pois)
    
    
    eta.z <- samples[[i]]$latent[idx.pred.z,1]
    p <- exp(eta.z)/(1 + exp(eta.z))
    
    
    eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
    xi <- samples[[i]]$hyperpar[1]
    # xi <- -0.284
    kappa <- samples[[i]]$hyperpar[2]
    
    pred.z <- rbinom(length(idx.pred.pois), size=1, prob=p)
    pred.cnt <- rpois(length(idx.pred.z), lambda.pois )
    pred.ba <- regp(length(idx.pred.ba), eta = eta.ba, kappa=kappa, xi=xi)
    
    zero.idx <- which(pred.z==0)
    pred.cnt[zero.idx] <- 0
    pred.ba[zero.idx] <- 0
    
    
    pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
    pred.z.mat <- cbind(pred.z.mat, pred.z)
    pred.ba.mat <- cbind(pred.ba.mat, pred.ba)
    
    param.cnt[[i]] <- lambda.pois
    param.z[[i]] <- p
    param.ba[[i]] <- eta.ba
    hyper.xi[[i]] <- xi
    hyper.kappa[[i]] <- kappa
  }
  
  
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:length(idx.pred.pois)){
    pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
    pred.z[[i]] <- as.numeric(pred.z.mat[i,])
    pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
  }
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z
              
              # 'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.xi'=hyper.xi,
              # 'hyper.kappa'=hyper.kappa
              ))
}



post.pred.gamma.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
  pred.cnt.mat <- c()
  pred.z.mat <- c()
  pred.ba.mat <- c()
  param.cnt <- list()
  param.z <- list()
  param.ba <- list()
  hyper.ba <- list()
  
  for (i in 1:n.samples){
    eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
    lambda.pois <- exp(eta.pois)
    
    
    eta.z <- samples[[i]]$latent[idx.pred.z,1]
    p <- exp(eta.z)/(1 + exp(eta.z))
    
    eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
    mu.ba <- exp(eta.ba)
    prec.par <- samples[[i]]$hyperpar[1]
    
    a = prec.par
    b = mu.ba / a
    
    pred.z <- rbinom(length(idx.pred.pois), size=1, prob=p)
    pred.cnt <- rpois(length(idx.pred.z), lambda.pois )
    pred.ba <- rgamma(length(idx.pred.ba), shape = a, scale = b)
    
    zero.idx <- which(pred.z==0)
    pred.cnt[zero.idx] <- 0
    pred.ba[zero.idx] <- 0
    
    
    pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
    pred.z.mat <- cbind(pred.z.mat, pred.z)
    pred.ba.mat <- cbind(pred.ba.mat, pred.ba)
    
    param.cnt[[i]] <- lambda.pois
    param.z[[i]] <- p
    param.ba[[i]] <- b
    hyper.ba[[i]] <- a
  }
  
  
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:length(idx.pred.pois)){
    pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
    pred.z[[i]] <- as.numeric(pred.z.mat[i,])
    pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
  }
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z
              # 'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 
              # 'hyper.ba'=hyper.ba
              ))
}





# 
idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)
 

####################################################################
# Sample and save results
####################################################################

t1 <- Sys.time()
# pred.sp <- post.pred.egp.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
# pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
pred.sp <- post.pred.gamma.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)

t2 <- Sys.time()
print(t2-t1)


save(pred.sp, file=file.path(dir.out,'Model_gamma_bym2_pred_h1.RData'))
# save(samples, pred.sp, file=file.path(dir.out,'Model_egp_bym2_pred_h1.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_bym2_pred_h1.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_egp_bym2_pred_h1_noxgb.RData'))