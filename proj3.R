#### Wenqi Ni s1792412

#### address of github repo:
# 

### Input the mcycle data from the MASS library
library(MASS)
mcycle_data <- force(mcycle)
#data <- c(mcycle_data["times"], mcycle_data["accel"])
x <- mcycle_data$times
y <- mcycle_data$accel


### 

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  n <- length(x) ## the size of the data or the number of observations
  ## Construct the D matrix
  D <- diff(diag(k),differences=pord)
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  ## Form the QR decomposition
  qrx <- qr(X)
  R <- qr.R(qrx); Q <- qr.Q(qrx)
  ## Eigen-decomposition of the covariance matrix
  S <- solve(t(R))%*%crossprod(D)%*%solve(R)
  ec <- eigen(S)
  U <- ec$vectors
  lambda <- ec$values
  lambda <- diag(lambda)
  ## Create an identity matrix
  I <- diag(k)
  ## Give a list of the values of lambda within the range
  lamb <- rep(0, ngrid)  
  effdegree <- rep(0, ngrid)
  sig2_list <- rep(0, ngrid)
  gcv <- rep(0, ngrid)
  coefval <- array(rep(0,k*ngrid), dim=c(k,ngrid))
  fittedval <- array(rep(0,n*ngrid), dim=c(n,ngrid))
  for (i in 1:ngrid) {
    lamb[i] <- logsp[1] + (i-1)*(logsp[2]-logsp[1])/(ngrid-1)
    lamb[i] <- exp(lamb[i]) ## transform from log lambda to lambda
    ## Find the effective degrees of freedom for each lambda
    effdegree[i] <- sum(diag(solve(I+lamb[i]*lambda)))
    coefval[,i] <- solve(R)%*%U%*%solve(I+lamb[i]*lambda)%*%t(U)%*%t(Q)%*%y
    fittedval[,i] <- X%*%coefval[,i]
    sig2_list[i] <- sum((y-fittedval[,i])^2) / (n-effdegree[i])
    gcv[i] <- sig2_list[i]/(n-effdegree[i]) 
  }  
  ## To find the optimal smoothing parameter
  ind <- order(gcv)[1]
  opt_lamb <- lamb[ind] 
  sig2 <- sig2_list[ind]
  coef <- coefval[ind]
  fitted <- fittedval[ind]
  ## Store all the elements into a vector
  result <- c(coef, fitted, sig2)
  return(result)
}

###
print.pspline <- function(m) {
  
  
}

