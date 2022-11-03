#### Wenqi Ni s1792412

#### address of github repo:
# 

#### Overview




## 'pspline'
## The 'pspline' function has inputs, x and y (data), k (the number of basis functions), logsp (the ends of the interval of smoothing parameters),
## bord (the B-spline order, 3 for cubic), pord(the order of difference to use in the penalty), and ngrid(the number of smoothing parameters).
## The function adopts QR decomposition (components: R, Q) and eigen decomposition (components: U, lambda) to fit B-splines to x, y data 
## and then use for-loop to find the optimal smoothing parameter, lamb
## The 'pspline' function will return a list of elements, including x and y(the data), bord, pord, D, X(the basis model matrix), knots(the vector to set up the B-spline basis),
## lamb (the best fit spline smoother), and the corresponding coef (the coefficients), fitted (the fitted values), sig2 (residual variance), edf(), and gcv(generalized cross validation criterion)

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
  R <- qr.R(qrx); Q <- qr.Q(qrx) ## extract components Q, R
  I <- diag(k) ## create an identity matrix
  ## Then define the eigen-decomposition
  S <- solve(t(R))%*%crossprod(D)%*%solve(R) ## covariance matrix
  ec <- eigen(S)## eigen decompose the S matrix
  U <- ec$vectors; lambda <- ec$values ## exract eigenvectors and values
  lambda <- diag(lambda) ## construct a diagonal matrix whose diagonal elements are eigenvalues
  ## Give empty lists for the smoothing parameter, lamb, and its corresponding effective degrees of freedom (effdegree), the residual covariance (sig2_list), 
  ## the generalized cross validation criterion (gcv), the coefficient (coefval), and the fitted values (fittedval).
  lamb <- c(); edf <- c(); sig2_list <- c(); gcv <- c()
  coefval <- matrix(0, nrow=k, ncol=ngrid); fittedval <- matrix(0, nrow=n, ncol=ngrid)
  for (i in 1:ngrid) { ## the number of lamb that needs to try is ngrid, 100
    lamb[i] <- logsp[1] + (i-1)*(logsp[2]-logsp[1])/(ngrid-1) ## log lamb values are within the interval whose ends are specified by 'logsp'
    lamb[i] <- exp(lamb[i]) ## transform from log lamb to lamb
    ## Find the effective degrees of freedom for each lamb and store in the 'effdegree' vector
    edf[i] <- sum(diag(solve(I+lamb[i]*lambda))) 
    coefval[,i] <- solve(R)%*%U%*%solve(I+lamb[i]*lambda)%*%t(U)%*%t(Q)%*%y
    fittedval[,i] <- X%*%coefval[,i]
    sig2_list[i] <- sum((y-fittedval[,i])^2) / (n-edf[i])
    gcv[i] <- sig2_list[i]/(n-edf[i]) 
  }  
  ## To find the index for the optimal smoothing parameter which minimises the 'gcv' values
  ind <- which(gcv %in% sort(gcv)[1])
  ## find the values of other variables corresponding to the best smoothing parameter
  opt_lamb <- lamb[ind]; sig2 <- sig2_list[ind]; coef <- coefval[,ind]; fitted <- fittedval[,ind]
  ## Store the best values into a named list
  m <- list(x,y, lamb[ind],bord, pord,D,X,U,R,I,knots,lambda,coef,fitted, sig2, edf[ind], gcv[ind] )
  names(m) <- c("x", "y","lamb", "bord","pord","D","X","U","R","I","knots","lambda","coef", "fitted", "sig2", "edf", "gcv")
  return(m)
}

m <- pspline(x,y)

## 'print.pspline'
## The 'print.pspline' function has inputs, m, which is a list of the results of the model fit
## The function reports some results of the model fit which can be obtained by using the 'pspline' function. 
## The 'print.pspline' function will return the order of the p-spline (bord), the order of difference in penalty (pord), the effective degrees of freedom (edf),
## the size of the coefficients (coeff), the residual standard deviation (rsd), the r-squared (r2), the GCV (gcv).

print.pspline <- function(m) {
  ## Extract the components from m
  edf <- m$edf; y <- m$y
  sig2 <- unlist(m$sig2)
  gcv <- unlist(m$gcv)
  rsd <- sqrt(sig2) ## calculate residual std dev
  n <- length(y)
  r2 <- 1-(n-1)*sig2/sum((y-sum(y)/length(y))^2) ## calculate the r-squared
  coeff <- length(unlist(m$coef)) ## the length of the coefficients obtained in m 
  out <- list(edf, r2,gcv)
  output <- cat("\n Order ", unlist(m$bord), " p-spline with order ", unlist(m$pord), " penalty \n")
  output <- cat(" Effective degrees of freedom: ",unlist(out[1]), " Coefficients: ",coeff, " \n residual std dev: ",rsd, " r-squared: ", unlist(out[2]), "  GCV: ", unlist(out[3]))
  #return(out) # invisible
  invisible(out) 
}

 
## 'predict.pspline'
## The function has inputs, m (the model fit), x (new x values within the range of the original x), and 'se=TRUE' 
## The function uses the results of the smooth fit (obtained using the 'pspline' function), including D, knots, and coef, to make predictions for new x values. 
## After extracting useful components from m, we firstly use 'splineDesign' to create a new model matrix based on the new x, 
## If 'se=FALSE' then the 'predict.pspline' function will return the vector of the predicted values.
## If 'se=TRUE' then a 2-item named list will be returned, the 'fit' item shows the predicted values and the 'se' item gives the corresponding standard errors

## new x values which are within the range of the original data


predict.pspline <- function(m,x,se=TRUE) {
  ## Extract elements from m
  D <- m$D; bord <- m$bord; knots <- m$knots; coef <- m$coef; R <- m$R; U <- m$U
  #D <- diff(diag(k),differences=pord)
  #bord <- unlist(m[5])
  #coef <- unlist(m[1])
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)## create a new model matrix for the new x data
  ## Covariance matrix
  ## Form the QR decomposition
  #qrx <- qr(X)
  #R <- qr.R(qrx); Q <- qr.Q(qrx) ## extract the components, Rp and Qp
  ## Eigen-decomposition of the covariance matrix
  ## the covariance matrix
  V <- solve(R)%*%(U)%*%solve(m$I+ (m$lamb)*(m$lambda))%*%t(U)%*%solve(t(R))*(m$sig2)
  #V <- solve(t(R))%*%crossprod(D)%*%solve(R)
  ## To calculate the standard errors: stanerr <- Xp %*% V %*% t(Xp)
  stanerr <- rowSums(X*(X%*%V))^.5 ## a more computable efficiently method
  Y <- X %*% coef   ## calculate the predicted values
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

yp <- predict.pspline(m,xp,se=TRUE)
t <- rowSums(Xp*(Xp%*%V))^.5

sqrt( sum( diag(se)) )

EY <- predict.pspline(m,x,se=TRUE)


## 'plot.pspline'
## The function has input, m (the model fit of the original data)
## what is done
## The 'plot.pspline' function will return three plots, 

plot.pspline <- function(m) {
  ## Extract useful components from m
  x <- m$x; y <- m$y
  ## new x values
  #x <- 
  ### alternative
  EY <- predict.pspline(m,x,se=TRUE) ## use the above function to
  #fitted EY
  ll <- EY$fit - 1.96*EY$se; ul <- EY$fit + 1.96*EY$se
  ## Store all the EY, ll, ul, x, and y values together
  values <- cbind(x,y,EY$fit,ll,ul)
  values <- values[order(x), ] ## reorder the data based on the x values.
  plot(x, values[,2], xlab="x: times",ylab="accel",pch=16,col='blue',main=" data")
  ## Add the estimated smooth function
  lines(x,values[,3], col="red", type="l")
  lines(x,values[,4], col="red", type="l")
  lines(x,values[,5],col="red")
  legend("topright", legend = c("fitted","ll","ul", "y"), col=c("red","red","red","blue"), lty=c(1,1,1,NA), pch=c(NA,NA,NA,16))
  
  ## The second plot: the model residuals against fitted values
  res <- y - EY$fit ## calculate the model residuals
  plot(EY$fit, res, xlab=" fitted values ",ylab="the model residuals ", pch=16)
  ## The third should be a qqplot of the residuals
  qqnorm(res) ## the normal Q-Q plot for the residuals
  #return(plot1, plot2, plot3)
}



plot(a,b,a,d, xlim=c(0,2), ylim=c(0,15))

plot(x, values[,2], xlab="x: times",xlim=c(0,60), ylim=c(-150,60),col='blue',main="Original x, y data")
plot(x, values[,3], xlab="x: times",xlim=c(0,60), ylim=c(-150,60),col='blue',main="Original x, y data", type="l", lty=2)
plot(x, values[,4], xlab="x: times",xlim=c(0,60), ylim=c(-150,60),col='blue',main="Original x, y data",type="l", lty=2)

## plot 1
plot(x, values[,2], xlab="x: times", col='blue',main="Original x, y data") 
abline(lm(EY~x)) 
abline(lm(values[,3]~x)) 
abline(lm(values[,4]~x)) 


## Input the mcycle data from the MASS library
library(MASS)
mcycle_data <- force(mcycle)
#data <- c(mcycle_data["times"], mcycle_data["accel"])
x <- mcycle_data$times
y <- mcycle_data$accel
