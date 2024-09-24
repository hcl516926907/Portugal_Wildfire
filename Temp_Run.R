dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)

inla.setOption(scale.model.default = TRUE)

load(file=file.path(dir.out,'Final_Model_3.3.RData'))
result <- res3.3
n1 <- 192*108

post.pred.gpd.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, alpha,n.samples=200 ){
  rgp = function(n, sigma, eta, alpha, xi = 0.001)
  {
    if (missing(sigma)) {
      stopifnot(!missing(eta) && !missing(alpha))
      sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
    }
    return (sigma / xi * (runif(n)^(-xi) -1.0))
  }
  t1 <- Sys.time()
  res <- foreach(j = 1:n1) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      xi <- samples[[i]]$hyperpar[1]
      
      lambda <- exp(eta.pois)
      
      pred.cnt[i] <- rpois(1, lambda[j] )
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        pred.ba[i] <- rgp(1,eta=eta.ba[j], alpha=alpha, xi=xi)
      }else{
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:n1){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}


post.pred.gpd.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, alpha,n.samples=200 ){
  rgp = function(n, sigma, eta, alpha, xi = 0.001)
  {
    if (missing(sigma)) {
      stopifnot(!missing(eta) && !missing(alpha))
      sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
    }
    return (sigma / xi * (runif(n)^(-xi) -1.0))
  }
  t1 <- Sys.time()
  res <- foreach(j = 1:n1) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      xi <- samples[[i]]$hyperpar[1]
      
      lambda <- exp(eta.pois)
      

      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        pred.ba[i] <- rgp(1,eta=eta.ba[j], alpha=alpha, xi=xi)
        pred.cnt[i] <- rpois(1, lambda[j] )
      }else{
        pred.ba[i] <- 0
        pred.cnt[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:n1){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}





post.pred.gamma.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  
  t1 <- Sys.time()
  res <- foreach(j = 1:n1) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      mu.ba <- exp(eta.ba)
      prec.par <- samples[[i]]$hyperpar[1]
      
      a = prec.par
      b = mu.ba / a
      
      
      lambda <- exp(eta.pois)
      
      pred.cnt[i] <- rpois(1, lambda[j] )
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        pred.ba[i] <- rgamma(1, shape = a, scale = b[j])
      }else{
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:n1){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}



post.pred.gamma.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  
  t1 <- Sys.time()
  res <- foreach(j = 1:n1) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      mu.ba <- exp(eta.ba)
      prec.par <- samples[[i]]$hyperpar[1]
      
      a = prec.par
      b = mu.ba / a
      
      
      lambda <- exp(eta.pois)
      

      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        pred.ba[i] <- rgamma(1, shape = a, scale = b[j])
        pred.cnt[i] <- rpois(1, lambda[j] )
      }else{
        pred.ba[i] <- 0
        pred.cnt[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:n1){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}


post.pred.weibull.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  res <- foreach(j = 1:n1) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      
      alpha.ba <- samples[[i]]$hyperpar[1]
      
      lambda.ba = exp(eta.ba)
      
      
      lambda <- exp(eta.pois)
      
      pred.cnt[i] <- rpois(1, lambda[j] )
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        # pred.ba[i] <- rgamma(1, shape = a, scale = b[j])
        pred.ba[i] <- rweibull(1, shape = alpha.ba, scale = lambda.ba[j]^(-1/alpha.ba))
      }else{
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:n1){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}


post.pred.weibull.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  res <- foreach(j = 1:n1) %dopar%{
    set.seed(j)
    pred.cnt <- rep(NA, n.samples)
    pred.z <- rep(NA, n.samples)
    pred.ba <- rep(NA, n.samples)
    for (i in 1:n.samples){
      eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
      
      eta.z <- samples[[i]]$latent[idx.pred.z,1]
      p <- exp(eta.z)/(1 + exp(eta.z))
      
      
      eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
      
      alpha.ba <- samples[[i]]$hyperpar[1]
      
      lambda.ba = exp(eta.ba)
      
      
      lambda.pois <- exp(eta.pois)
      
      
      
      z <- rbinom(1, size=1, prob=p[j])
      pred.z[i] <- z
      if (z==1){
        # pred.ba[i] <- rgamma(1, shape = a, scale = b[j])
        pred.cnt[i] <- rpois(1, lambda.pois[j] )
        pred.ba[i] <- rweibull(1, shape = alpha.ba, scale = lambda.ba[j]^(-1/alpha.ba))
      }else{
        pred.cnt[i] <- 0
        pred.ba[i] <- 0
      }
    }
    res.list <- list('pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba)
    
  }
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:n1){
    pred.cnt[[i]] <- res[[i]][['pred.cnt']]
    pred.z[[i]] <- res[[i]][['pred.z']]
    pred.ba[[i]] <- res[[i]][['pred.ba']]
  }
  
  
  return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
}

library(doParallel)
cl <- detectCores()
print(cl)
registerDoParallel(cl)

n.samples = 200
# n.samples = 1000


idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)




samples = inla.posterior.sample(n.samples, result = result, seed=1234)


t1 <- Sys.time()
# pred.sp <- post.pred.gamma.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
# pred.sp <- post.pred.gamma.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
# pred.sp <- post.pred.gpd.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, alpha=0.5, n.samples=n.samples)
# pred.sp <- post.pred.gpd.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, alpha=0.5, n.samples=n.samples)
# pred.sp <- post.pred.weibull.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
t2 <- Sys.time()
print(t2-t1)

save(pred.sp, file=file.path(dir.out,'Final_Model_3.3_pred.sp_200.RData'))
