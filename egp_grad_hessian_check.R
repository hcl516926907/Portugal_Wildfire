library(numDeriv)
##########xi !=0
func_1 <- function(x){
  sigma <- x[1]
  y <- x[2]
  kappa <- x[3]
  xi <- x[4]

  y1 <- 1+xi*y/sigma
  ll <- log(kappa) + (kappa-1)*log(1-y1^(-1/xi))-log(sigma)-(1+1/xi)*log(y1)
  return(-ll)
}

point <- c(2,1,2,0.2)

sigma <- point[1]
y <- point[2]
kappa <-point[3]
xi <- point[4]

y1 <- 1+xi*y/sigma
dl_dsigma <- (1 / sigma) - ((xi + 1) * y) / (sigma^2 * y1) + 
  ((kappa - 1) * y * y1^(-1 / xi - 1)) / (sigma^2 * (1 - y1^(-1 / xi)))
print(dl_dsigma)
grad(func_1, point)[1]




term1 <- -1 / sigma^2

term2 <- 2 * (xi + 1) * y / (sigma^3 * y1)
term3 <- (xi + 1) * xi * y^2 / (sigma^4 * y1^2)

term4_num <- 2 * y * y1^(-1 / xi - 1)
term4_den <- sigma^3 * (1 - y1^(-1 / xi))
term4 <- -(kappa - 1) * term4_num / term4_den

term5_num <- (1 +  xi) * y^2 * y1^(-1 / xi - 2)
term5_den <- sigma^4 * (1 - y1^(-1 / xi))
term5 <- (kappa - 1) * term5_num / term5_den

term6_num <-  y^2 * y1^(-2 / xi - 2)
term6_den <- sigma^4 * (1 - y1^(-1 / xi))^2
term6 <- (kappa - 1) * term6_num / term6_den

d2l_dsigma2 <- term1 + term2 - term3 + term4 + term5 + term6
print(d2l_dsigma2)
diag(hessian(func_1, point))[1]


##########xi ==0

func_2 <- function(x){
  sigma <- x[1]
  y <- x[2]
  kappa <- x[3]
  
  num <- log(kappa) + (kappa-1)*log(1-exp(-y/sigma))
  den <- log(sigma) + y/sigma
  
  ll <- num - den
  return(-ll)
}


point <- c(3,2,2)
sigma <- point[1]
y <- point[2]
kappa <-point[3]
y2 <- exp(-y/sigma)

dl_dsigma <- (1 / sigma) - y / sigma^2 + 
  ((kappa - 1) * y * y2) / (sigma^2 * (1 - y2))
print(dl_dsigma)
grad(func_2, point)[1]

term1 <- -1 / sigma^2
term2 <- 2 * y / sigma^3

term3 <- (kappa - 1) * y2*y^2 / (sigma^4 * (1 - y2))

term4_num <- 2*sigma*y2*(1-y2) - y2^2*y
term4_den <- sigma^4*(1-y2)^2

term4 <- (kappa-1)*y*term4_num/term4_den
d2l_dsigma2 <- term1 + term2 + term3 - term4
print(d2l_dsigma2)

diag(hessian(func_2, point))[1]