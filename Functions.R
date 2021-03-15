################################################################################
#
#                                                                                     
#   Filename    :	      Functions.R    												  
#                                                                                     
#   Project     :       EJS article "Optimal designs for testing the efficacy 
#                       of heterogeneous experimental groups"                                                             
#   Authors     :       A. Baldi Antognini, R.Frieri, M. Novelli and M. Zagoraiou                                                                
#   Date        :       15.03.2021
#   Purpose     :       Functions used for implementing optimal designs and 
#                       their smoothed versions
#																				  
#   R Version   :       R-4.0.4 (2021-02-15)                                                               
#   Platform    :       x86_64-w64-mingw32/x64 (64-bit)
#
#   Input file  :       ---                                                       
#   Output data files :   ---
#
#   Required R packages :  ---
#
#
################################################################################


#IMPORTANT NOTE: Before using the functions, save the content of the code folder 
#                in the working directory so that the DLL can be properly loaded


## ----------------------------------------------------------------------------
## Description: Unconstrained optimal design \tilde{\rho} in Theorem 2.1
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## Usage: uncon_target(theta, v)
##
##        theta:     treatment effects
##        v:         corresponding variances
## 
## ----------------------------------------------------------------------------
## Value: Unconstrained optimal design maximizing NCP of Wald test
## ----------------------------------------------------------------------------
uncon_target=function(theta, v){
  K=length(theta)
  if(var(theta) == 0) return(rep(1/K, K))
  grid=t(combn(1:K,2))
  AUX=cbind(grid, NA)
  for(i in 1:nrow(grid)) AUX[i, 3]=(((theta[grid[i,1]]-theta[grid[i,2]])/(sqrt(v[grid[i,1]])+sqrt(v[grid[i,2]]))))^2
  idx=AUX[order(AUX[,3], decreasing = T),][1,1:2]
  neyman=sqrt(v[idx[1]])/(sqrt(v[idx[1]])+sqrt(v[idx[2]]))
  rho=rep(0, K)
  rho[idx[1]]=neyman
  rho[idx[2]]=1-neyman
  return(rho)
}

## ----------------------------------------------------------------------------
## Description: constrained optimal design \rho^* in Theorem 2.2
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## Usage: con_target(theta, v)
##
##        theta:     treatment effects
##        v:         corresponding variances
## 
## ----------------------------------------------------------------------------
## Value: constrained optimal design maximizing NCP of Wald test
## ----------------------------------------------------------------------------
con_target <- function(theta, v){
  K <- Kt <- length(theta)
  if(var(theta)==0){return(rep(1/K, K))}
  ord=order(theta,decreasing=TRUE)
  theta=theta[ord]
  v=v[ord]
  x <- v^-.5
  y <- (theta[1] - theta)*x
  T1 <- (mean(x^2) * mean(y^2)) / mean(x*y)^2
  A <- P1 <- numeric(K-1)
  P1b <- matrix(F, nrow=K-1, ncol = K-1)
  sj <- bj <- gj <- numeric()
  for(i in 1:(K-1)){
    sj[i] <- mean(x[1:i]^2)/mean(x^2)
    bj[i] <- mean(x[1:i]*y[1:i])/mean(x*y)
    gj[i] <- mean(y[1:i]^2) / mean(y^2)
    det <- (1 - bj[i])^2 - T1*(1 - sj[i])* (1 - gj[i]) 
    if( det  < 0 | is.nan(det) ){
      A[i] = 0
    }else{
      A[i] <- T1*(1 - gj[i]) /  ((1 - bj[i]) + sqrt(det) )
      P1[i] <- A[i]>=bj[i]/sj[i] & A[i] <= 1 
    }
  }
  P1[is.na(P1)] <- 0
  sel <- 1:(K-1)
  for(j in sel){
    for(l in sel[-j])
      P1b[j,l] <- A[j]^2*(1 - sj[l]) - 2*A[j]*(1 - bj[l]) + T1*(1 - gj[l])
  }
  P1b[abs(P1b) < 10^-15] <- 0
  # check for balanced allocation
  bal <- numeric(K-1)
  for(i in sel){# condition P2
    bal[i] <- T1*(1 - gj[i]) - 1 - sj[i] + 2*bj[i] > 0
  }
  if(sum(bal) == K - 1){
    rho <- rep(1/K, K)# balanced allocation 
    return(rho)
  }else{ 
    
    if(sum(P1)>0){# check condition P1
      index <- which(P1>0)
      cc1 <- numeric()
      for(i in index){
        cc1[i] <- sum(P1b[i,] >= 0)} # check condition P1b
      ind <- which(cc1 == K-1)[1]  
      if(!is.na(ind)){
        tau <- (A[ind]*sj[ind] - bj[ind])/(K*(1-bj[ind]-A[ind]*(1 - sj[ind])))
        rho <- c(rep((1-(K-ind)*tau)/ind, ind),rep(tau, K-ind))
        rho[ord]=rho # constrained optimal design: equation 2.7
        return(rho)
      }
    }
    # check for boundary solutions
    f <- K
    alt <- 0
    while(alt == 0 & f > 2){
      theta <- theta[-f]
      v <- v[-f]
      K <- length(theta)
      x <- v^-.5
      y <- (theta[1] - theta)*x
      T1 <- (mean(x^2) * mean(y^2)) / mean(x*y)^2
      A <- P1 <- numeric(K-1)
      P1b <- matrix(F, nrow=K-1, ncol = K-1)
      sj <- bj <- gj <- numeric()
      for(i in 1:(K-1)){
        sj[i] <- mean(x[1:i]^2)/mean(x^2)
        bj[i] <- mean(x[1:i]*y[1:i])/mean(x*y)
        gj[i] <- mean(y[1:i]^2) / mean(y^2)
        det <- (1 - bj[i])^2 - T1*(1 - sj[i])* (1 - gj[i])
        if( det  < 0  ){
          A[i] = 0
        }else{
          A[i] <- T1*(1 - gj[i]) /  ((1 - bj[i]) + sqrt(det) )
          P1[i] <- A[i]>=bj[i]/sj[i] & A[i] <= 1
        }
      }
      sel <- 1:(K-1)
      for(j in sel){
        for(l in sel[-j])
          P1b[j,l] <- A[j]^2*(1 - sj[l]) - 2*A[j]*(1 - bj[l]) + T1*(1 - gj[l])
      }
      # check for balanced allocation of boundary solutions
      bal <- numeric(K-1)
      for(i in sel){
        bal[i] <- T1*(1 - gj[i]) - 1 - sj[i] + 2*bj[i] > 0
      }
      if(sum(bal) == K - 1){
        rho <- c(rep(1/K, K), rep(0, Kt-f+1))
        rho[ord]=rho
        return(rho)
        alt <- 1
      }else{
        if(sum(P1) == K-2){
          ind <- which(P1==1)[1]
          tau <- (A[ind]*sj[ind] - bj[ind])/(K*(1-bj[ind]-A[ind]*(1 - sj[ind])))
          rho <- c(rep((1-(K-ind)*tau)/ind, ind),rep(tau, K-ind), rep(0,Kt-f+1))
          rho[ord]=rho
          return(rho)
          alt <- 1
        }else if(sum(P1) > 0){
          index <- which(P1>0)
          cc1 <- numeric()
          for(i in index){
            cc1[i] <- sum(P1b[i,] >= 0)}
          ind <- which(cc1 == K-sum(P1))[1]  
          if(!is.na(ind)){
            tau <- (A[ind]*sj[ind] - bj[ind])/(K*(1-bj[ind]-A[ind]*(1 - sj[ind])))
            rho <- c(rep((1-(K-ind)*tau)/ind, ind),rep(tau, K-ind), rep(0,Kt-f+1))
            rho[ord]=rho
            alt <- 1
            return(rho)
          }
        }
      }
    }
  }
}


## ----------------------------------------------------------------------------
## Description: Discrete Gaussian smoothing kernel for K = 3 and K = 4
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## Usage: kernel.gauss(theta, v)
##
##        ksize:  number of smoothing points
##        sigma:  covariance matrix of the Gaussian kernel
## 
## ----------------------------------------------------------------------------
## Value: the Gaussian smoothin kernel
## ----------------------------------------------------------------------------
# smoothing kernel for K = 4
kernel.gauss<-function(voxdim = c(1, 1, 1, 1), ksize = 7, sigma = diag(3, 4))  {
  a <- array(0, dim = c(ksize, ksize, ksize, ksize))
  centre <- (ksize + 1) / 2
  sig.inv <- solve(sigma)
  sig.det <- abs(det(sigma))
  for(i in 1:ksize) {
    for(j in 1:ksize) {
      for(k in 1:ksize) {
        for(z in 1:ksize) {
          x <- (c(i, j, k, z) - centre) * voxdim
          a[i, j, k, z] <- ((2 * pi)^(-4 / 2)) * exp(-.5 * (t(x) %*% sig.inv %*% x)) / sqrt(sig.det)
        }
      }
    }
  }  
  a <- a / sum(a)
  return(a)
}
# smoothing kernel for K = 3
kernel.gauss_3<-function(voxdim = c(1, 1, 1), ksize = 7, sigma = diag(3, 3))  {
  a <- array(0, dim = c(ksize, ksize, ksize))
  centre <- (ksize + 1) / 2
  sig.inv <- solve(sigma)
  sig.det <- abs(det(sigma))
  for(i in 1:ksize) {
    for(j in 1:ksize) {
      for(k in 1:ksize) {
        x <- (c(i, j, k) - centre) * voxdim
        a[i, j, k] <- ((2 * pi)^(-3 / 2)) * exp(-.5 * (t(x) %*% sig.inv %*% x)) / sqrt(sig.det)
      }
    }
  }  
  a <- a / sum(a)
  return(a)
}




## ----------------------------------------------------------------------------
## Description: Smoothed optimal design for RAR implementation
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## Usage: targets_smoother(theta.hat, target, ksize, mod, v, m1)
##
##    theta.hat:  estimated treatment effects
##
##    target:     target allocation to be smoothed, either uncon_target or con_target.
##                The user may provide a different target to be smoothed in the form
##                of a function with two arguments, theta, and v
##                
##    ksize:      odd number of smoothing points, higher values provide a more 
##                accurate smoothing at the expense of higher computational time; 
##                a value between 11 and 25 can provide a good compromise
##
##    mod:        observations distribution, one of "bin", "norm", "exp" or "pois"
##
##    v:          variances, to be used only with mod = "norm"
##
##    m1:         smoothing range, to be used only for mod = ""norm"
## 
## ----------------------------------------------------------------------------
## Value: smoothed optimal designs
## ----------------------------------------------------------------------------

targets_smoother = function(theta.hat,
                            target,
                            ksize,
                            mod,
                            v = NULL,
                            m1 = NULL, 
                            ...) {
  models <- c("bin", "norm", "exp", "pois")
  K <- length(theta.hat)
  if(!(mod %in% models)) stop("mod mus be on of bin, norm, exp or pois")
  if(K != 3 & K != 4) stop("K must be either 3 or 4")
  if(!is.loaded("conv3")){
    if(Sys.info()['sysname'] == "Windows"){dyn.load(paste0(getwd(), "/conv.dll"))
  }else if (Sys.info()['sysname'] == "Linux"){dyn.load(paste0(getwd(), "/conv.so"))
  }else{dyn.load(paste0(getwd(), "/conv1.so"))}}
  if(!is.loaded("conv3")) stop("Operating system not supported")
  if(length(formals(target))<2) stop("the target to be smothed must be a function of at least two arguments")
  if(ksize %% 2 == 0)ksize = ksize + 1
  cent <- (ksize + 1) / 2
  xx <- yy <- list()
  init <- c(0.02, 0.01, 0.015, 0.017)
  if (mod == "bin") {
    cc1 <- (ksize - 1) / 2
    for(i in 1:K){
      xx[[i]] <- c(seq(-theta.hat[i] + init[i],-theta.hat[i]/cc1,length.out = cc1), 0, seq(theta.hat[i]/cc1 ,theta.hat[i],length.out = cc1))
      yy[[i]] <- ifelse(theta.hat[i]+xx[[i]] <= 0, init[i], ifelse(theta.hat[i]+xx[[i]]>=1, 1-init[i], theta.hat[i]+xx[[i]]))
    }
    grid <- expand.grid(yy)
    res <- t(apply(grid, 1, function(x) target(x, x * (1-x), ...)))
  } else if (mod == "norm") {
    if(is.null(v)) v <- rep(1,K)
    if(is.null(m1)) m1 <- max(theta.hat) - min(theta.hat)
    if(length(v) != K) stop("the dimension of v must match that of theta.hat")
    x <- seq(-m1, m1, length.out = ksize)
    for(i in 1:K){
      xx[[i]] <- theta.hat[i] + x + init[i]
    }
    grid <- as.matrix(expand.grid(xx))
    res <- t(apply(grid, 1, function(x) target(x, v, ...) ))
  } else{
    cc1 <- (ksize - 1) / 2
    for(i in 1:K){
      xx[[i]] <- theta.hat[i] + c(seq(-theta.hat[i] + init[i],-theta.hat[i]/cc1,length.out = cc1), 0, seq(theta.hat[i]/cc1 , theta.hat[i],length.out = cc1))
    }
    grid <- as.matrix(expand.grid(xx))
    if (mod == "exp") {
      res <- t(apply(grid, 1, function(x) target(x, x^2, ...))) 
    } else if (mod == "pois") {
      res <- t(apply(grid, 1, function(x) target(x, x, ...))) 
    }
  }
  if(K == 3){
    filtmat <- kernel.gauss_3(ksize = ksize)
    out <-
      array(
        .Fortran(
          "conv3",
          as.double(res),
          as.integer(ksize),
          as.integer(ksize),
          as.integer(ksize),
          as.integer(K),
          as.double(filtmat),
          as.integer(ksize),
          as.double(res)
        )[[8]],
        dim = c(rep(ksize, K), K)
      )[cent, cent, cent,] 
  }else{
    filtmat <- kernel.gauss(ksize = ksize)
    out <- array(
      .Fortran(
        "conv4",
        as.double(res),
        as.integer(ksize),
        as.integer(ksize),
        as.integer(ksize),
        as.integer(ksize),
        as.integer(K),
        as.double(filtmat),
        as.integer(ksize),
        as.double(res)
      )[[9]],
      dim = c(rep(ksize, K), K)
    )[cent, cent, cent, cent, ]
  }
  return(out)
}
