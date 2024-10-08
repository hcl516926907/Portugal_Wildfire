dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)
library(doParallel)

print(detectCores())


inla.setOption(scale.model.default = TRUE)

# load(file=file.path(dir.out,'Model_weibull_0.0625.RData'))
# load(file=file.path(dir.out,'Model_weibull_bym2_0.125.RData'))
load(file=file.path(dir.out, 'Model_weibull_spde_0.125.RData'))
# load(file=file.path(dir.out, 'Model_weibull_spde_0.25.RData'))
result <- res
rm(res)

load(file=file.path(dir.out, 'data.fit2.pred_0.125.RData'))

# n.grid <- 2554
# n.grid <- 192
# n.grid <- 681
# n1 <- n.grid*108
# pred.time <- c(97, 108)
# subset.idx <- (1+(pred.time[1]-1)*n.grid) : (n.grid + (pred.time[2]-1)*n.grid )
# 
# post.pred.weibull.hurdle.parallel <- function(samples, idx.pred.pois, idx.pred.z, idx.pred.ba,n.samples=200 ){
#   t1 <- Sys.time()
#   # res <- foreach(j = 1:length(idx.pred.pois) ) %dopar%{
#   batch.size <- c(seq(1,length(idx.pred.pois),5000), length(idx.pred.pois))
#   res <- c()
#   for (k in 1:(length(batch.size)-1)){
#     itr <- batch.size[k]:batch.size[k+1]
#     res.batch <- foreach(j = itr ) %dopar%{
#       set.seed(j)
#       pred.cnt <- rep(NA, n.samples)
#       pred.z <- rep(NA, n.samples)
#       pred.ba <- rep(NA, n.samples)
#       
#       for (i in 1:n.samples){
#         eta.pois <- samples[[i]]$latent[idx.pred.pois,1]
#         
#         eta.z <- samples[[i]]$latent[idx.pred.z,1]
#         p <- exp(eta.z)/(1 + exp(eta.z))
#         
#         
#         eta.ba <- samples[[i]]$latent[idx.pred.ba,1]
#         
#         alpha.ba <- samples[[i]]$hyperpar[1]
#         
#         lambda.ba = exp(eta.ba)
#         
#         lambda.pois <- exp(eta.pois)
#         
#         z <- rbinom(1, size=1, prob=p[j])
#         pred.z[i] <- z
#         if (z==1){
#           pred.cnt[i] <- rpois(1, lambda.pois[j] )
#           pred.ba[i] <- rweibull(1, shape = alpha.ba, scale = lambda.ba[j]^(-1/alpha.ba))
#         }else{
#           pred.cnt[i] <- 0
#           pred.ba[i] <- 0
#         }
#       }
#       res.list <- list('grid.idx'=j,'pred.cnt'=pred.cnt, 'pred.z'=pred.z,'pred.ba'=pred.ba,
#                        'param.z' = p, 'param.cnt' = lambda.pois, 'param.ba'=lambda.ba^(-1/alpha.ba),
#                        'hyper.ba'=alpha.ba)
#       
#     }
#     res <- c(res,res.batch)
#     rm(res.batch)
#   }
#   
#   t2 <- Sys.time()
#   print(t2-t1)
#   pred.cnt <- list()
#   pred.z <- list()
#   pred.ba <- list()
#   for (i in 1:length(idx.pred.pois)){
#     pred.cnt[[i]] <- res[[i]][['pred.cnt']]
#     pred.z[[i]] <- res[[i]][['pred.z']]
#     pred.ba[[i]] <- res[[i]][['pred.ba']]
#   }
#   
#   
#   return(list('pred.ba'=pred.ba,'pred.cnt'=pred.cnt, 'pred.z'=pred.z))
# }


post.pred.weibull.hurdle.parallel <- function(samples,n.samples ){
  t1 <- Sys.time()

  latent.names <- sub(":.*","",rownames(samples[[1]]$latent))
  idx.time.idx1 <- which(latent.names=='time.idx1')
  idx.time.idx2 <- which(latent.names=='time.idx2')
  idx.time.idx3 <- which(latent.names=='time.idx3')
  idx.score_1 <- which(latent.names=='score_1')
  idx.score_2 <- which(latent.names=='score_2')
  idx.score_3 <- which(latent.names=='score_3')
  idx.spat1 <- which(latent.names=='spat1')
  idx.spat2 <- which(latent.names=='spat2')
  idx.spat3 <- which(latent.names=='spat3')
  idx.Intercept1 <- which(latent.names=='Intercept1')
  idx.Intercept2 <- which(latent.names=='Intercept2')
  idx.Intercept3 <- which(latent.names=='Intercept3')
  
  pred.cnt.mat <- c()
  pred.z.mat <- c()
  pred.ba.mat <- c()
  param.cnt <- list()
  param.z <- list()
  param.ba <- list()
  hyper.ba <- list()
  for (i in 1:n.samples){
    set.seed(i)
    
    month.idx <- data.fit2.pred$month.idx
    
    eta.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2.pred$month.idx] + 
      A1.pred%*%samples[[i]]$latent[idx.score_1] + A.spat.1.pred%*%samples[[i]]$latent[idx.spat1]
    lambda.pois <- as.numeric(exp(eta.pois))

    eta.z <- samples[[i]]$latent[idx.Intercept2] + samples[[i]]$latent[idx.time.idx2][data.fit2.pred$month.idx] + 
      A2.pred%*%samples[[i]]$latent[idx.score_2] + A.spat.2.pred%*%samples[[i]]$latent[idx.spat2]
    p <- as.numeric(exp(eta.z)/(1 + exp(eta.z)))
    
    eta.ba <- samples[[i]]$latent[idx.Intercept3] + samples[[i]]$latent[idx.time.idx3][data.fit2.pred$month.idx] + 
      A3.pred%*%samples[[i]]$latent[idx.score_3] + A.spat.3.pred%*%samples[[i]]$latent[idx.spat3]
    
    alpha.ba <- as.numeric(samples[[i]]$hyperpar[1])
    lambda.ba = as.numeric(exp(eta.ba))
    scale.ba <- lambda.ba^(-1/alpha.ba)
    

    
    pred.z <- rbinom(nrow(data.fit2.pred), size=1, prob=p)


    pred.cnt <- rpois(nrow(data.fit2.pred), lambda.pois )
    pred.ba <- rweibull(nrow(data.fit2.pred), shape = alpha.ba, scale = scale.ba)
    
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

  for (i in 1:nrow(data.fit2.pred)){
    pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
    pred.z[[i]] <- as.numeric(pred.z.mat[i,])
    pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
  }
  
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
              'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.ba'=hyper.ba))
}


n.samples = 200
# n.samples = 1000

# 
# idx.pred.pois <- 1:n1
# idx.pred.z <- (n1+1):(2*n1)
# idx.pred.ba <- (2*n1+1):(3*n1)
# 
# 
# idx.pred.pois <- idx.pred.pois[subset.idx]
# idx.pred.z <- idx.pred.z[subset.idx]
# idx.pred.ba <- idx.pred.ba[subset.idx]


samples = inla.posterior.sample(n.samples, result = result, seed=1234)

# latent.names <- sub(":.*","",rownames(samples[[1]]$latent))
# idx.time.idx1 <- which(latent.names=='time.idx1')
# idx.time.idx2 <- which(latent.names=='time.idx2')
# idx.time.idx3 <- which(latent.names=='time.idx3')
# idx.score_1 <- which(latent.names=='score_1')
# idx.score_2 <- which(latent.names=='score_2')
# idx.score_3 <- which(latent.names=='score_3')
# idx.spat1 <- which(latent.names=='spat1')
# idx.spat2 <- which(latent.names=='spat2')
# idx.spat3 <- which(latent.names=='spat3')
# idx.Intercept1 <- which(latent.names=='Intercept1')
# idx.Intercept2 <- which(latent.names=='Intercept2')
# idx.Intercept3 <- which(latent.names=='Intercept3')
# 

# 
# test.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2$month.idx] +
#   A1%*%samples[[i]]$latent[idx.score_1] + A.spat.1%*%samples[[i]]$latent[idx.spat1]
# test.pois-samples[[i]]$latent[idx.pred.pois]
# 
# test.z <- samples[[i]]$latent[idx.Intercept2] + samples[[i]]$latent[idx.time.idx2][data.fit2$month.idx] +
#   A2%*%samples[[i]]$latent[idx.score_2] + A.spat.2%*%samples[[i]]$latent[idx.spat2]
# 
# test.z-samples[[i]]$latent[idx.pred.z]
# 
# test.ba <- samples[[i]]$latent[idx.Intercept3] + samples[[i]]$latent[idx.time.idx3][data.fit2$month.idx] +
#   A3%*%samples[[i]]$latent[idx.score_3] + A.spat.3%*%samples[[i]]$latent[idx.spat3]
# test.ba-samples[[i]]$latent[idx.pred.ba]

t1 <- Sys.time()
# pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
pred.sp <- post.pred.weibull.hurdle.parallel(samples, n.samples=n.samples)
t2 <- Sys.time()
print(t2-t1)



# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.125_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.0625_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.25_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.125_pred_0.0625_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_0.25_pred_0.0625_pred_sp_200.RData'))
save(pred.sp, file=file.path(dir.out,'Model_weibull_0.125_pred_0.125_pred_sp_200.RData'))


