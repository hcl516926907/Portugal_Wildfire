dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Urban_Fires"
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Urban_Fire_Model'


library(INLA)
library(doParallel)

print(detectCores())


inla.setOption(scale.model.default = TRUE)


# load(file=file.path(dir.out, 'Model_weibull_spde_0.125.RData'))
# load(file=file.path(dir.out, 'Model_egp_spde_0.125.RData'))
load(file=file.path(dir.out, 'Model_gamma_spde_0.125.RData'))

result <- res
rm(res)

# load(file=file.path(dir.out, 'data.fit2.pred_0.125.RData'))

# n.grid <- 2554
# n.grid <- 192
n.grid <- 681
n1 <- n.grid*108
pred.time <- c(97, 108)
subset.idx <- (1+(pred.time[1]-1)*n.grid) : (n.grid + (pred.time[2]-1)*n.grid )


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

  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
              'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.ba'=hyper.ba))
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
    # xi <- samples[[i]]$hyperpar[1]
    xi <- -0.488
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
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
              'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.xi'=hyper.xi,
              'hyper.kappa'=hyper.kappa ))
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
  
  return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
              'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.ba'=hyper.ba))
}


######SPDE prediction##############
# post.pred.weibull.hurdle.parallel <- function(samples,n.samples ){
#   t1 <- Sys.time()
# 
#   latent.names <- sub(":.*","",rownames(samples[[1]]$latent))
#   idx.time.idx1 <- which(latent.names=='time.idx1')
#   idx.time.idx2 <- which(latent.names=='time.idx2')
#   idx.time.idx3 <- which(latent.names=='time.idx3')
#   # idx.score_1 <- which(latent.names=='score_1')
#   # idx.score_2 <- which(latent.names=='score_2')
#   # idx.score_3 <- which(latent.names=='score_3')
#   idx.score_1 <- which(latent.names=='score_cnt_bin')
#   idx.score_2 <- which(latent.names=='score_z_bin')
#   idx.score_3 <- which(latent.names=='score_ba_bin')
#   idx.spat1 <- which(latent.names=='spat.cnt')
#   idx.spat2 <- which(latent.names=='spat.z')
#   idx.spat3 <- which(latent.names=='spat.ba')
#   idx.Intercept1 <- which(latent.names=='Intercept1')
#   idx.Intercept2 <- which(latent.names=='Intercept2')
#   idx.Intercept3 <- which(latent.names=='Intercept3')
# 
#   pred.cnt.mat <- c()
#   pred.z.mat <- c()
#   pred.ba.mat <- c()
#   param.cnt <- list()
#   param.z <- list()
#   param.ba <- list()
#   hyper.ba <- list()
#   for (i in 1:n.samples){
#     set.seed(i)
# 
#     month.idx <- data.fit2.pred$month.idx
# 
#     # eta.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2.pred$month.idx] +
#     #   A1.pred%*%samples[[i]]$latent[idx.score_1] + A.spat.1.pred%*%samples[[i]]$latent[idx.spat1]
#     eta.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2.pred$month.idx] +
#       samples[[i]]$latent[idx.score_1][data.fit2.pred$score_cnt_bin] + A.spat.1.pred%*%samples[[i]]$latent[idx.spat1]
#     lambda.pois <- as.numeric(exp(eta.pois))
# 
#     # eta.z <- samples[[i]]$latent[idx.Intercept2] + samples[[i]]$latent[idx.time.idx2][data.fit2.pred$month.idx] +
#     #   A2.pred%*%samples[[i]]$latent[idx.score_2] + A.spat.2.pred%*%samples[[i]]$latent[idx.spat2]
#     eta.z <- samples[[i]]$latent[idx.Intercept2] + samples[[i]]$latent[idx.time.idx2][data.fit2.pred$month.idx] +
#       samples[[i]]$latent[idx.score_2][data.fit2.pred$score_z_bin] + A.spat.2.pred%*%samples[[i]]$latent[idx.spat2]
# 
#     p <- as.numeric(exp(eta.z)/(1 + exp(eta.z)))
# 
#     # eta.ba <- samples[[i]]$latent[idx.Intercept3] + samples[[i]]$latent[idx.time.idx3][data.fit2.pred$month.idx] +
#     #   A3.pred%*%samples[[i]]$latent[idx.score_3] + A.spat.3.pred%*%samples[[i]]$latent[idx.spat3]
#     eta.ba <- samples[[i]]$latent[idx.Intercept3] + samples[[i]]$latent[idx.time.idx3][data.fit2.pred$month.idx] +
#       samples[[i]]$latent[idx.score_3][data.fit2.pred$score_ba_bin] + A.spat.3.pred%*%samples[[i]]$latent[idx.spat3]
#     
#     alpha.ba <- as.numeric(samples[[i]]$hyperpar[1])
#     lambda.ba = as.numeric(exp(eta.ba))
#     scale.ba <- lambda.ba^(-1/alpha.ba)
# 
# 
# 
#     pred.z <- rbinom(nrow(data.fit2.pred), size=1, prob=p)
# 
# 
#     pred.cnt <- rpois(nrow(data.fit2.pred), lambda.pois )
#     pred.ba <- rweibull(nrow(data.fit2.pred), shape = alpha.ba, scale = scale.ba)
# 
#     zero.idx <- which(pred.z==0)
#     pred.cnt[zero.idx] <- 0
#     pred.ba[zero.idx] <- 0
# 
# 
#     pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
#     pred.z.mat <- cbind(pred.z.mat, pred.z)
#     pred.ba.mat <- cbind(pred.ba.mat, pred.ba)
# 
#     param.cnt[[i]] <- lambda.pois
#     param.z[[i]] <- p
#     param.ba[[i]] <- scale.ba
#     hyper.ba[[i]] <- alpha.ba
#   }
# 
# 
#   t2 <- Sys.time()
#   print(t2-t1)
# 
# 
#   pred.cnt <- list()
#   pred.z <- list()
#   pred.ba <- list()
# 
#   for (i in 1:nrow(data.fit2.pred)){
#     pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
#     pred.z[[i]] <- as.numeric(pred.z.mat[i,])
#     pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
#   }
# 
# 
#   return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
#               'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.ba'=hyper.ba))
# }



#####BYM2 prediction##############
# post.pred.weibull.hurdle.parallel <- function(samples,n.samples ){
#   t1 <- Sys.time()
# 
#   latent.names <- sub(":.*","",rownames(samples[[1]]$latent))
#   idx.time.idx1 <- which(latent.names=='time.idx1')
#   idx.time.idx2 <- which(latent.names=='time.idx2')
#   idx.time.idx3 <- which(latent.names=='time.idx3')
#   idx.score_1 <- which(latent.names=='score_1')
#   idx.score_2 <- which(latent.names=='score_2')
#   idx.score_3 <- which(latent.names=='score_3')
#   idx.spat1 <- which(latent.names=='idarea1')[1:n.grid]
#   idx.spat2 <- which(latent.names=='idarea2')[1:n.grid]
#   idx.spat3 <- which(latent.names=='idarea3')[1:n.grid]
#   idx.Intercept1 <- which(latent.names=='Intercept1')
#   idx.Intercept2 <- which(latent.names=='Intercept2')
#   idx.Intercept3 <- which(latent.names=='Intercept3')
# 
#   pred.cnt.mat <- c()
#   pred.z.mat <- c()
#   pred.ba.mat <- c()
#   param.cnt <- list()
#   param.z <- list()
#   param.ba <- list()
#   hyper.ba <- list()
#   for (i in 1:n.samples){
#     set.seed(i)
# 
#     month.idx <- data.fit2.pred$month.idx
# 
#     eta.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2.pred$month.idx] +
#        samples[[i]]$latent[idx.spat1] + A1.pred%*%samples[[i]]$latent[idx.score_1] 
#     lambda.pois <- as.numeric(exp(eta.pois))
# 
#     eta.z <- samples[[i]]$latent[idx.Intercept2] + samples[[i]]$latent[idx.time.idx2][data.fit2.pred$month.idx] +
#       samples[[i]]$latent[idx.spat2] + A2.pred%*%samples[[i]]$latent[idx.score_2] 
#     p <- as.numeric(exp(eta.z)/(1 + exp(eta.z)))
# 
#     eta.ba <- samples[[i]]$latent[idx.Intercept3] + samples[[i]]$latent[idx.time.idx3][data.fit2.pred$month.idx] +
#        samples[[i]]$latent[idx.spat3] + A3.pred%*%samples[[i]]$latent[idx.score_3] 
# 
#     alpha.ba <- as.numeric(samples[[i]]$hyperpar[1])
#     lambda.ba = as.numeric(exp(eta.ba))
#     scale.ba <- lambda.ba^(-1/alpha.ba)
# 
# 
# 
#     pred.z <- rbinom(nrow(data.fit2.pred), size=1, prob=p)
# 
# 
#     pred.cnt <- rpois(nrow(data.fit2.pred), lambda.pois )
#     pred.ba <- rweibull(nrow(data.fit2.pred), shape = alpha.ba, scale = scale.ba)
# 
#     zero.idx <- which(pred.z==0)
#     pred.cnt[zero.idx] <- 0
#     pred.ba[zero.idx] <- 0
# 
# 
#     pred.cnt.mat <- cbind(pred.cnt.mat, pred.cnt)
#     pred.z.mat <- cbind(pred.z.mat, pred.z)
#     pred.ba.mat <- cbind(pred.ba.mat, pred.ba)
# 
#     param.cnt[[i]] <- lambda.pois
#     param.z[[i]] <- p
#     param.ba[[i]] <- scale.ba
#     hyper.ba[[i]] <- alpha.ba
#   }
# 
# 
#   t2 <- Sys.time()
#   print(t2-t1)
# 
# 
#   pred.cnt <- list()
#   pred.z <- list()
#   pred.ba <- list()
# 
#   for (i in 1:nrow(data.fit2.pred)){
#     pred.cnt[[i]] <- as.numeric(pred.cnt.mat[i,])
#     pred.z[[i]] <- as.numeric(pred.z.mat[i,])
#     pred.ba[[i]] <- as.numeric(pred.ba.mat[i,])
#   }
# 
# 
#   return(list('pred.ba'=pred.ba, 'pred.cnt'=pred.cnt, 'pred.z'=pred.z,
#               'param.cnt'=param.cnt, 'param.z'=param.z, 'param.ba'=param.ba, 'hyper.ba'=hyper.ba))
# }


n.samples = 200
# n.samples = 1000

# 
idx.pred.pois <- 1:n1
idx.pred.z <- (n1+1):(2*n1)
idx.pred.ba <- (2*n1+1):(3*n1)
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


# test.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2$month.idx] +
#   samples[[i]]$latent[idx.spat1] + samples[[i]]$latent[idx.score_1][data.fit2$score_cnt_bin]
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

# i <- 1
# latent.names <- sub(":.*","",rownames(samples[[1]]$latent))
# idx.time.idx1 <- which(latent.names=='time.idx1')
# idx.time.idx2 <- which(latent.names=='time.idx2')
# idx.time.idx3 <- which(latent.names=='time.idx3')
# # idx.score_1 <- which(latent.names=='score_1')
# # idx.score_2 <- which(latent.names=='score_2')
# # idx.score_3 <- which(latent.names=='score_3')
# idx.score_1 <- which(latent.names=='score_cnt_bin')
# idx.score_2 <- which(latent.names=='score_z_bin')
# idx.score_3 <- which(latent.names=='score_ba_bin')
# idx.spat1 <- which(latent.names=='spat.cnt')
# idx.spat2 <- which(latent.names=='spat.z')
# idx.spat3 <- which(latent.names=='spat.ba')
# idx.Intercept1 <- which(latent.names=='Intercept1')
# idx.Intercept2 <- which(latent.names=='Intercept2')
# idx.Intercept3 <- which(latent.names=='Intercept3')
# test.pois <- samples[[i]]$latent[idx.Intercept1] + samples[[i]]$latent[idx.time.idx1][data.fit2.pred$month.idx] +
#   samples[[i]]$latent[idx.score_1][data.fit2.pred$score_cnt_bin] + A.spat.1.pred%*%samples[[i]]$latent[idx.spat1]
# print(summary( as.numeric(test.pois-samples[[i]]$latent[idx.pred.pois][subset.idx])))

t1 <- Sys.time()
# pred.sp <- post.pred.weibull.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
# pred.sp <- post.pred.weibull.hurdle.parallel(samples, n.samples=n.samples)
# pred.sp <- post.pred.egp.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)
pred.sp <- post.pred.gamma.hurdle.parallel(samples, idx.pred.pois, idx.pred.z, idx.pred.ba, n.samples=n.samples)

t2 <- Sys.time()
print(t2-t1)



# save(pred.sp, file=file.path(dir.out,'Model_weibull_bym2_0.125_pred_0.125_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_weibull_spde_0.125_pred_0.125_pred_sp_200.RData'))
# save(pred.sp, file=file.path(dir.out,'Model_egp_spde_0.125_pred_0.125_pred_sp_200.RData'))
save(pred.sp, file=file.path(dir.out,'Model_gamma_spde_0.125_pred_0.125_pred_sp_200.RData'))
