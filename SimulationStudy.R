library(gtools)
library(spc)
library(compositions)
SMW <- function(alpha, alpha1, n, lambda = 0.05, B = 5000){
  p <- length(alpha)
  A <- diag(trigamma(alpha))
  a <- trigamma(sum(alpha))
  one <- as.vector(rep(1,p))
  Sigma <- ( A-a*one%*%t(one) )/n
  SigmaE <- (lambda/(2-lambda))*Sigma
  ARL <- rep(NA, B)
  H <- mewma.crit(l = lambda, L0 = 200, p = p)
  for (b in 1:B){
    j <- 0
    Qj <- 0
    while (Qj<H) {
      x <- gtools::rdirichlet(n=n, alpha=alpha1 )
      xstar <- colMeans(log(x))
      Sj <- digamma(sum(alpha))-digamma(alpha)+xstar
      if (j==0){ 
        Ej<-lambda*Sj
      }else{
        Ej<-lambda*Sj+(1-lambda)*Ej  
      }
      Qj <- t(Ej)%*%solve(SigmaE)%*%Ej
      j <- j+1
    }
    ARL[b]<-j
  }
  mean(ARL)
}
MW <- function(alpha, alpha1, n, lambda = 0.05, H = 7.3565, B = 5000){
  ARL <- rep(NA, B)
  p <- length(alpha)
  mu0 <- exp(digamma(alpha))/sum(exp(digamma(alpha)))
  mu0star <- as.numeric(-1*ilr(mu0))
  ##### calculate covariance matrix
  var0 <- diag(((p-1)/p)^2*trigamma(alpha)+(1/p)^2*(sum(trigamma(alpha))-trigamma(alpha)))
  A <- matrix(0,p,p)
  for (i in 1:p){
    for(j in 1:p){
      A[i,j] = -((p-1)/p^2)*(trigamma(alpha)[i]+trigamma(alpha)[j])+(1/p)^2*sum(trigamma(alpha)[-c(i,j)])
    }
    A[i,i] = 0
  }
  var_clr <- A + var0
  U <- ilrBase(alpha)
  Sigmastar <- t(U)%*%var_clr%*%U
  SigmaYi <- lambda/(n*(2-lambda))*Sigmastar
  #####
  for (b in 1:B){
    Qi <- 0
    i <- 1
    while( Qi < H){
      x <- gtools::rdirichlet(n=n, alpha=alpha1 )
      z <- matrix( -1*as.numeric(ilr(x)), ncol = (p-1))
      xbarstar <- colMeans(z)
      if (i==1){ 
        Yi<-lambda*(xbarstar-mu0star)
      }else{
        Yi<-lambda*(xbarstar-mu0star)+(1-lambda)*Yi
      }
      Qi <- t(Yi)%*%solve(SigmaYi)%*%Yi
      i <- i+1
    }
    ARL[b] <- i  
  }
  return(mean(ARL))
}
MC <- function(alpha, alpha1, n, k = 0.5, H = 5.5, B = 5000){
  ARL <- rep(NA, B)
  p <- length(alpha)
  mu0 <- exp(digamma(alpha))/sum(exp(digamma(alpha)))
  mu0star <- as.numeric(-1*ilr(mu0))
  #####
  var0 <- diag(((p-1)/p)^2*trigamma(alpha)+(1/p)^2*(sum(trigamma(alpha))-trigamma(alpha)))
  A <- matrix(0,p,p)
  for (i in 1:p){
    for(j in 1:p){
      A[i,j] = -((p-1)/p^2)*(trigamma(alpha)[i]+trigamma(alpha)[j])+(1/p)^2*sum(trigamma(alpha)[-c(i,j)])
    }
    A[i,i] = 0
  }
  var_clr <- A + var0
  U <- ilrBase(alpha)
  Sigmastar <- t(U)%*%var_clr%*%U
  #####
  p <- p-1
  Sinv <-  solve(Sigmastar/n)
  #####
  for (b in 1:B){
    #Qt <- 0
    Yt <- 0
    i <- 1
    while( Yt < H){
      x <- gtools::rdirichlet(n=n, alpha=alpha1 )
      z <- matrix( -1*as.numeric(ilr(x)), ncol = p)
      if (n==1) xbarstar <- as.vector(z) else xbarstar <- colMeans(z)
      dif <- xbarstar-mu0star
      ###
      if (i==1) {
        Ct <- as.numeric(sqrt( t(dif)%*%Sinv%*%(dif) ))
        if (Ct  <= k)  St <- rep(0, p) else St <- dif*(1-k/Ct)
        Yt <- as.numeric(sqrt(( t(St)%*%Sinv%*%St)))
      } else {
        Ct <- as.numeric( sqrt(t(dif+St)%*%Sinv%*%(dif+St) ) )
        if (Ct  <= k)  St <- rep(0, p) else St <- (dif+St)*(1-k/Ct)
        Yt <- as.numeric(sqrt(( t(St)%*%Sinv%*%St )))
      }
      i <- i+1
    }
    ARL[b] <- i  
  }
  return(mean(ARL))
}
# in-control
B <- 5000
phi <- 150
alpha0 <- c(1/3, 1/3, 1/3)*phi
alpha0 <- c(1/3, 1/3, 1/3)*phi
n <- 5
SMW(alpha = alpha0, alpha1 = alpha0, n = n, lambda = 0.05, B = B)
MW(alpha = alpha0, alpha1 = alpha0, n = n, lambda = 0.05, H = 7.3565, B = B)
MC(alpha = alpha0, alpha1 =  alpha0, n = n, H = 5.5, k = 0.5, B = B)

# out-of-control
B <- 5000
phi0 <- 120
alpha0 <- c(1/3, 1/3, 1/3)*phi0
phi1 <- 150
alpha1 <- c(0.26, 0.294, 0.446)*phi1
n <- 5
SMW(alpha = alpha0, alpha1 = alpha1, n = n, lambda = 0.05, B = B)
MW(alpha = alpha0, alpha1 = alpha1, n = n, lambda = 0.05, H = 7.3565, B = B)
MC(alpha = alpha0, alpha1 = alpha1, n = n, H = 5.5, k = 0.5, B = B)
