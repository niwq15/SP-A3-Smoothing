#### Wenqi Ni s1792412

#### address of github repo:
# 

#### Overview

## Input the mcycle data from the MASS library
library(MASS)
mcycle_data <- force(mcycle)
#data <- c(mcycle_data["times"], mcycle_data["accel"])
x <- mcycle_data$times
y <- mcycle_data$accel


## The 'pspline' function
## 

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  n <- length(x) ## the size of the data x,y, or the number of observations
  ## Construct the D matrix
  D <- diff(diag(k),differences=pord)
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  ## Apply the splineDesign function from R’s built in splines package
  ## to set up the basis, X matrix
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  ## Form the QR decomposition of the X matrix
  qrx <- qr(X)
  R <- qr.R(qrx); Q <- qr.Q(qrx) ## exract components Q, R
  ## Then define the eigen-decomposition
  S <- solve(t(R))%*%crossprod(D)%*%solve(R) ## covariance matrix
  ec <- eigen(S)## eigen decompose the S matrix
  U <- ec$vectors; lambda <- ec$values ## exract eigenvectors and values
  lambda <- diag(lambda) ## construct a diagonal matrix whose diagonal elements are eigenvalues
  I <- diag(k) ## create an identity matrix
  ## Give empty lists for the smoothing parameter, lamb, and its corresponding effective degrees of freedom (k), 
  lamb <- c(); 
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
  ## To find the index for the optimal smoothing parameter which minimises the 'gcv' values
  ind <- which(gcv %in% sort(gcv)[1])
  opt_lamb <- lamb[ind] 
  sig2 <- sig2_list[ind]
  coef <- coefval[,ind]
  fitted <- fittedval[,ind]
  ## Store all the elements into a vector
  #result <- c(coef, fitted, sig2, knots, effdegree[ind] )
  m <- list(coef, fitted, sig2, knots, bord, lamb[ind],effdegree[ind], gcv[ind] )
  names(m) <- c("coef", "fitted", "sig2", "knots", "bord","lamb", "effdegree", "gcv")
  return(m)
}

m <- pspline(x,y)

r_squared <- 1-(n-1)*sig2/sum((y-sum(y)/length(y))^2)

###
print.pspline <- function(m) {
  edf <- m[6] ## effective degrees of freedom
  sig2 <- unlist(m[3])
  rsd <- sqrt(sig2) ## residual std dev
  r2 <- 1-(n-1)*sig2/sum((y-sum(y)/length(y))^2) ## r-squared
  gcv <- m[7]
  coeff <- length(unlist(m[1])) ## ? k if edf?
  values <- list(edf, rsd,r2, gcv)
  output <- cat("\n", "Order 3 p-spline with order 2 penalty")
  #cat("\n", "Effective degrees of freedom: ", edf, "     Coefficients: ")
  #cat("\n", "residual std dev: ", rsd, "     r-squared: ", r2, "     GCV: ", gcv)
  output <- list("Effective degrees of freedom: "= edf, "     Coefficients: "= coeff,"residual std dev: "=rsd, "     r-squared: "= r2, "  GCV: "= gcv)
  return(output)
}

###
## new x values which are within the range of the original data
predict.pspline <- function(m,x,se=TRUE) {
  D <- diff(diag(k),differences=pord)
  #bord <- unlist(m[5])
  #coef <- unlist(m[1])
  ## Create a new model matrix for the new x data
  Xp <- splines::splineDesign(knots,xp,ord=bord+1,outer.ok=TRUE)
  ## Covariance matrix
  ## Form the QR decomposition
  qrxp <- qr(Xp)
  Rp <- qr.R(qrxp); Qp <- qr.Q(qrxp)
  ## Eigen-decomposition of the covariance matrix
  V <- solve(t(Rp))%*%crossprod(D)%*%solve(Rp)
  ## Standard errors 
  #se <- Xp %*% V %*% t(Xp)
  stanerr <- rowSums(Xp*(Xp%*%V))^.5
  ## Predicted values
  Y <- Xp %*% coef
  if (se==TRUE) {
    ## Return a 2-item named list
    itemlist <- list(fit= Y, se=stanerr) 
    return(itemlist)
  } else { # if se= FALSE
    ## Show the vector of predictions for the new x
    return(Y)
  }
}

### check for the new x which is a subset of original x
## Create new x values within the range of original data
xp <- x[1:25]
xp <- sample(min(x):max(x), length(x), replace = TRUE)
mp <- predict.pspline(m,xp,se=TRUE)

t <- rowSums(Xp*(Xp%*%V))^.5

sqrt( sum( diag(se)) )


###

plot.pspline <- function(m) {
  ## The first plot: 
  plot1 <- plot(x,y, xlab="x: times", ylab="y: accel",xlim=c(0,60),ylim=c(-150,80), main="Original 'mcycle' data")
  ## 
  EY <- X%*%(m$coef)
  ## Add the estimated smooth function
  abline(lm(EY~x), col="blue") 
  
  mu = unlist(m[2]); sig2 <- unlist(m[3])
  th <- cbind(x,y)
  ## The 95% credible intervals
  ll <- mu - 1.95*EY
  ul <- mu + 1.95*EY
  plot1 <-plot(x,ll, xlim=c(0,60)) 
  plot1 <- ggplot2(data=th, aes(x=times, y=accel)) + geom_line(aes(y=ll), color="blue") +geom_line(aes(y=ul), color="blue")
  #ci <- apply(th,1,quantile,prob=c(.025,.975))
  
  ## The second plot: the model residuals against fitted values
  plot(EY, sig2)
  
  ## The third should be a qqplot of the residuals
  
  return(plot1, plot2, plot3)
}
