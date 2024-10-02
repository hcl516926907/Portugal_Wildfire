dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)

inla.setOption(scale.model.default = TRUE)

# load(file=file.path(dir.out,'Model_weibull_0.0625.RData'))
load(file=file.path(dir.out,'Model_weibull_bym2_0.125.RData'))
# load(file=file.pathC(dir.out,'Final_Model_3.1.RData'))
result <- res

# n.grid <- 2554
# n.grid <- 192
n.grid <- 681
n1 <- n.grid*108
pred.time <- c(97, 108)
subset.idx <- (1+(pred.time[1]-1)*n.grid) : (n.grid + (pred.time[2]-1)*n.grid )

post.pred.weibull.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
  t1 <- Sys.time()
  # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
  batch.size <- c(seq(1,length(idx.pred.pois),5000), length(idx.pred.pois))
  res <- c()
  for (k in 1:(length(batch.size)-1)){
    itr <- batch.size[k]:batch.size[k+1]
    res.batch <- foreach(j = itr ) %dopar%{
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
          pred.cnt[i] <- rpois(1, lambda.pois[j] )
          pred.ba[i] <- rweibull(1, shape = alpha.ba, scale = lambda.ba[j]^(-1/alpha.ba))
        }else{
          pred.cnt[i] <- 0
          pred.ba[i] <- 0
        }
      }
      res.list <- list('grid.idx'=j,'pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba,
                       'param.z' = p, 'param.cnt' = lambda.pois, 'param.ba'=lambda.ba^(-1/alpha.ba),
                       'hyper.ba'=alpha.ba)
      
    }
    res <- c(res,res.batch)
    rm(res.batch)
  }
  
  t2 <- Sys.time()
  print(t2-t1)
  pred.cnt <- list()
  pred.z <- list()
  pred.ba <- list()
  for (i in 1:length(idx.pred.pois)){
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


idx.pred.pois <- idx.pred.pois[subset.idx]
idx.pred.z <- idx.pred.z[subset.idx]
idx.pred.ba <- idx.pred.ba[subset.idx]


samples = inla.posterior.sample(n.samples, result = result, seed=1234)


t1 <- Sys.time()
pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
t2 <- Sys.time()
print(t2-t1)

save(pred.sp, file=file.path(dir.out,'Model_weibull_0.125_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.0625_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.25_pred_sp_200.RData'))