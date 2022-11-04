#### Wenqi Ni s1792412

#### address of github repo:
# https://github.com/niwq15/SP-A3-Smoothing

#### Overview
# This file holds the self-contained code to write R functions for smoothing x, y data.
# The main logic is to combine basis expansions and penalized regression, as follows:
# Consider a model, y_i = f(x_i) + error_i, with the observed x, y data, f as an unknown smooth function, and error_i as a zero mean 
# error term with some variance. To estimate f, we will use evenly spaced B-spline basis functions which consist of 2 parts that need
# to be estimated, basis functions and the coefficients. The basis functions are bell shaped functions, and they are well constructed that 
# the evaluations at any x value within the range of the data sum to 1. To avoid over fitting, we use a smoothing penalty to encourage the 
# coefficients corresponding to the neighbouring basis functions to vary smoothly from one to the next. The model is estimated by penalized 
# least squares and we choose a smoothing parameter, lamb, which minimized the generalized cross validation criterion (GCV). 
# To obtain the best lamb, we mainly need to use the knots (the vector to to set up the B-spline basis), X (the basis), some matrices such as 
# Q and R (can be obtained from the QR decomposition X = QR), U and lambda (obtained from eigen-decomposition), D (will be defined later), 
# I (an identity matrix), S (covariance matrix), edf (the effective degrees of freedom), sig2 (the residual variance), coef (the coefficients), 
# and fitted (the fitted values). 
# After obtaining the optimal lamb, we will use the best fit spline smoother to make predictions on new x values and draw some plots. 

# The following are some useful formulas to calculate those parameters/variables. 1. X = QR; 2. fitted = Xcoef
# 3. S = sig2*((X^T)X+lamb*(D^T)D)^-1 = sig2*(R^-1)U((I+lamb*lambda)^-1)(U^T)(R^-T));
# 4. coef = (((X^T)X+lamb*(D^T)D)^-1)(X^T)y = (R^-1)U((I+lamb*lambda)^-1)(U^T)(Q^T)y = (((X^T)X+lamb*(D^T)D)^-1)(R^T)(Q^T)y
# 5. UΛ(U^T) = (R^−T)(D^T)D(R^−1); 6. sig2 = (||y-fitted||^2)/(n-edf); 7. edf = tr{(I+lamb*lambda)^-1} = tr{(((X^T)X+lamb*(D^T)D)^-1)(X^T)X}.
# More details can be seen in Professor Simon Wood's notes, "Basis Penalty Smoothers"

# Firstly, we will show a 'pspline' function which returns an object (a list) defining the best fit spline smoother, then create a 
# 'print.pspline' function to report some details such as edf and r-squared of the model fit which are obtained in the 'pspline' function,
# then use a 'predict.pspline' to make predictions for new x values (within the range of the original data) based on the smooth fit. 
# Finally, a 'plot.pspline' function is created to give 3 plots: Plot1 shows the original x,y data and the overlaid estimated smooth function,
# together with a 95% credible intervals for the smooth; Plot2 gives the plot of the model residuals against the fitted values; Plot3 is 
# a qqplot of the model residuals. 


## 'pspline'
## The 'pspline' function has inputs, x and y (data), k (the number of basis functions), logsp (the ends of the interval of smoothing parameters),
## bord (the B-spline order, 3 for cubic), pord(the order of difference to use in the penalty), and ngrid(the number of smoothing parameters).
## The function uses a splineDesign function to set up the basis (X matrix), then adopts QR decomposition (components: R, Q) and eigen decomposition 
## (components: U, lambda), then use a for-loop to calculate the values of those variables corresponding to different smoothing paramters, 
## finally find the optimal smoothing parameter, lamb, which minimizes the generalized cross validation criterion. 
## The 'pspline' function will return a list of elements, including x and y, bord, pord, D, X(the basis model matrix), U, R, I(an identity matrix), 
## knots(the vector to set up the B-spline basis), lamb (the best smoothing paramter), and the corresponding coef (the coefficients), 
## fitted (the fitted values), sig2 (residual variance), edf, gcv(generalized cross validation criterion), and S (the covariance matrix).

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  n <- length(x) ## the size of the data x,y, or the number of observations
  D <- diff(diag(k),differences=pord) ## Construct the D matrix
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  ## Apply the splineDesign function from R’s built in splines package to set up the basis, X matrix
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  ## Form the QR decomposition of the X matrix
  qrx <- qr(X)
  R <- qr.R(qrx); Q <- qr.Q(qrx) ## extract components Q, R
  I <- diag(k) ## create an identity matrix
  ## Then define the eigen-decomposition, using the formula UΛ(U^T) = (R^−T)(D^T)D(R^−1)
  Sed <- solve(t(R))%*%crossprod(D)%*%solve(R)
  ec <- eigen(Sed)## eigen decompose the S matrix
  U <- ec$vectors; lambda <- ec$values ## extract eigenvectors and values
  lambda <- diag(lambda) ## construct a diagonal matrix whose diagonal elements are eigenvalues
  ## Create empty lists to store the smoothing parameter, lamb, and its corresponding edf, sig2, gcv, coef, fitted
  lamb <- c(); edf <- c(); sig2_list <- c(); gcv <- c()
  coefval <- matrix(0, nrow=k, ncol=ngrid); fittedval <- matrix(0, nrow=n, ncol=ngrid)
  for (i in 1:ngrid) { ## the number of lamb that needs to try is ngrid, 100
    ## The log lamb values are within the interval whose ends are specified by 'logsp'
    lamb[i] <- logsp[1] + (i-1)*(logsp[2]-logsp[1])/(ngrid-1) ## even spacing of values on the log scale.
    lamb[i] <- exp(lamb[i]) ## transform from log lamb to lamb
    ## The following calculations are based on the formulas mentioned in the overview
    edf[i] <- sum(diag(solve(I+lamb[i]*lambda))) ## find the effective degrees of freedom for each lamb and store in the 'effdegree' vector
    coefval[,i] <- solve(R)%*%U%*%solve(I+lamb[i]*lambda)%*%t(U)%*%t(Q)%*%y  ## calculate the coefficients
    fittedval[,i] <- X%*%coefval[,i] ## the fitted values, fitted = Xcoef
    sig2_list[i] <- sum((y-fittedval[,i])^2) / (n-edf[i]) ## calculate the residual variance 
    gcv[i] <- sig2_list[i]/(n-edf[i]) ## calculate the generalized cross validation criterion
  }  
  ## To find the index for the optimal smoothing parameter which minimizes the 'gcv' values
  ## The 'sort' function helps to reorder the element from smallest to largest and hence the first element after reordering is the smallest
  ind <- which(gcv %in% sort(gcv)[1])
  ## find the values of other variables corresponding to the best smoothing parameter
  opt_lamb <- lamb[ind]; sig2 <- sig2_list[ind]; coef <- coefval[,ind]; fitted <- fittedval[,ind] 
  S <- solve(R)%*%U%*%solve(I+lamb*lambda)%*%t(U)%*%solve(t(R))*sig2   ## calculate the covariance matrix in a computable efficiently way
  #S <- solve(crossprod(m$X)+m$lamb*crossprod(m$D))*m$sig2 (alternative formula)
  ## Store the values into a named list, pspline
  pspline <- list(x,y, opt_lamb,bord, pord,D,X,U,R,I,knots,lambda,coef,fitted, sig2, edf[ind], gcv[ind],S)
  names(pspline) <- c("x", "y","lamb", "bord","pord","D","X","U","R","I","knots","lambda","coef", "fitted", "sig2", "edf", "gcv","S")
  return(pspline)
}

## 'print.pspline'
## The 'print.pspline' function has inputs, m, which is a list of the results of the model fit
## The function reports some results of the model fit which can be obtained by using the 'pspline' function. 
## The 'print.pspline' function will return the order of the p-spline (bord), the order of difference in penalty (pord), 
## the effective degrees of freedom (edf), the size of the coefficients (coeff), the residual standard deviation (rsd), 
## the r-squared (r2), the GCV (gcv).

print.pspline <- function(m) {
  ## Extract the components from m
  y <- m$y
  sig2 <- m$sig2
  rsd <- sqrt(sig2) ## calculate residual std dev
  n <- length(y) ## the number of x or y observations
  r2 <- 1-(n-1)*sig2/sum((y-sum(y)/length(y))^2) ## calculate the r-squared
  coeff <- length(m$coef) ## the length of the coefficients obtained in m 
  ## Add a named list containing edf, r2, and gcv
  out <- list(m$edf, r2,m$gcv)
  names(out) <- c("edf", "r2", "gcv")
  output <- cat("\n Order ", m$bord, " p-spline with order ", m$pord, " penalty \n")
  output <- cat(" Effective degrees of freedom: ",unlist(out[1]), " Coefficients: ",coeff, " \n residual std dev: ",rsd, " r-squared: ", unlist(out[2]), "  GCV: ", unlist(out[3]))
  #return(out) 
  invisible(out) # silently return the 'out' list
}

 
## 'predict.pspline'
## The function has inputs, m (the model fit), x (new x values within the range of the original x), and 'se=TRUE' 
## The function uses the results of the smooth fit (obtained using the 'pspline' function), such as D, knots, and coef, to make predictions 
## for new x values. After extracting useful components from m, we firstly adopt 'splineDesign' to create a new model matrix based on the 
## new x, then use the results of the smooth fit to obtain the covariance matrix for the coefficients, then use a computable efficiently method
## to calculate the standard errors.
## If 'se=FALSE' then the 'predict.pspline' function will return the vector of the predicted values. If 'se=TRUE' then a 2-item named list 
## will be returned, the 'fit' item shows the predicted values and the 'se' item gives the corresponding standard errors

predict.pspline <- function(m,x,se=TRUE) {
  ## Create a new model matrix for the new x data
  Xp <- splines::splineDesign(m$knots,x,ord=m$bord+1,outer.ok=TRUE)
  V <- solve(m$R)%*%(m$U)%*%solve(m$I+ (m$lamb)*(m$lambda))%*%t(m$U)%*%solve(t(m$R))*(m$sig2) ## the covariance matrix; can also directly from the pspline function
  ## A general method to calculate the standard errors: stanerr <- Xp %*% V %*% t(Xp)
  stanerr <- rowSums(Xp*(Xp%*%V))^.5 ## a more computable efficiently method
  Y <- Xp %*% (m$coef)   ## calculate the predicted values
  if (se==TRUE) {
    ## Return a 2-item named list
    itemlist <- list(fit= Y, se=stanerr) # the 'fit' item stores the predicted values; the 'se' item stores the corresponding standard errors
    return(itemlist)
  } else { # if se= FALSE
    ## Show the vector of predictions for the new x
    return(Y)
  }
}

### check for the new x which is a subset of original x
## Create new x values within the range of original data
## xp <- runif(length(x), min=min(x), max=max(x) ) 


## 'plot.pspline'
## The function has input, m (the model fit of the original data)
## The function utilizes the components of m and the 'predict.pspline' funtion to make predictions on the fitted values, EY (a 2-item list 
## of 'fit' and 'se'), and uses the EY's 'fit' and 'se' values to calculate the lower and upper limits for the 95% credible intervals, 
## ll and ul. Before plotting, we also reorder the data of x, y, fit, ll and ul according to the x values. Then the model residuals will be 
## calculated by the difference between the original y values and the predicted values. The residuals will used for both the second and the third plots. 
## The 'plot.pspline' function will return 3 plots: The first is a plot of the original x,y data, with the estimated smooth function overlaid 
## as a line, together with a 95% credible intervals for the smooth. The second is for the model residuals and the fitted values. 
## The third is a qqplot of the model residuals. 

## Note: for the first plot, we can also use some new x values (within the range of the original x) to plot the estimated smooth function. 

plot.pspline <- function(m) {
  x <- m$x; y <- m$y ## extract useful components from m
  #xnew <- runif(length(x), min=min(x), max=max(x) ) ## new x values within the range of the original x
  ## The first plot
  ## Use the 'predict.pspline' funtion to make predictions
  EY <- predict.pspline(m,x,se=TRUE) ## 2 parts: EY$fit (predicted values), EY$se (the corresponding standard errors)
  ll <- EY$fit - 1.96*EY$se ## the lower confidence limit
  ul <- EY$fit + 1.96*EY$se ## the upper confidence limit
  values <- cbind(x,y,EY$fit,ll,ul) ## store all the x, y, EY, ll, ul values together
  values <- values[order(x), ] ## reorder the data based on the x values.
  ## Create a named list only containing elements, ll, ul, and x
  out <- list(values[,1], values[,4], values[,5])
  names(out) <- c("x","ll","ul")
  plot(values[,1], values[,2], xlab="x: times",ylab="y: accel",pch=16,col='blue',main="Original x, y data with the estimated \n smooth function and 95% credible intervals")
  ## Add the estimated smooth function and the credible intervals
  lines(values[,1],values[,3], col="black", type="l")
  lines(values[,1],values[,4], col="green", type="l",lty=2)
  lines(values[,1],values[,5],col="red",type="l",lty=2)
  legend("bottomright", legend = c("Original y values","Fitted values","Lower confidence limit","Upper confidence limit"), col=c("blue","black","green","red"), lty=c(NA,1,2,2), cex=0.8, pch=c(16,NA,NA,NA))
  ## The plot illustrates that the model is a good fit and it often has relatively small 95% credible intervals,
  ## but the intervals become larger as x tends to its boundary values
  
  ## The second plot: the model residuals against fitted values
  res <- values[,2] - values[,3] ## calculate the model residuals using the data after the reordering process
  plot(values[,3], res, xlab="Fitted values ",ylab="Residuals ", pch=16, main="Residuals vs Fitted values")
  ## It suggests no change to the current model.
  
  ## The third plot: a qqplot of the residuals
  qqnorm(res) ## the normal Q-Q plot for the residuals
  ## It indicates that the model does not violate the normality assumption. 
  invisible(out) ## silently return a list containing elements, ll, ul, and x
}




######## The following code can input the mcycle data from the MASS library
##library(MASS)
##mcycle_data <- force(mcycle)
##data <- c(mcycle_data["times"], mcycle_data["accel"])
##x <- mcycle_data$times
##y <- mcycle_data$accel
##m <- pspline(x,y)
