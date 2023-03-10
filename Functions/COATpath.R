#----------------------------------------------------------------------------------------
#  COAT estimate
#  Input:
#           x ------ n x p composition data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding
#        time ------ execution time
#----------------------------------------------------------------------------------------

coatPath <- function(X, Xval, soft = 1){
  startTime <- proc.time()
  p <- ncol(X)
  clrX <- log(X) - rowSums(log(X)) %*%matrix(1,1,p) / p
  clrXval <- log(Xval) - rowSums(log(Xval)) %*%matrix(1,1,p) / p
  coatPred <- adaptThresoldCov(clrX, clrXval, soft = soft)
  sigma <- coatPred$sigma
  corr <- coatPred$corr
  exeTimeClass <- proc.time() - startTime
  exeTime <- as.numeric(exeTimeClass[3])
  return(list(sigma = sigma, corr = corr, time = exeTime))
}

#----------------------------------------------------------------------------------------
#  Adaptive thresholding estimation of cov(x)
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding 
#----------------------------------------------------------------------------------------

adaptThresoldCov <- function(X, Xval, soft = 1){
  n <- nrow(X)
  p <- ncol(X)
  # Set the grid for the choice of tuning parameter
  nGrid <- 100
  gridInfo <- adaptThresholdRange(X)
  grid <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nGrid)/nGrid
  # Multi-folder cross validation
  error <- rep(0, nGrid)
  xTest <- Xval
  xTrain <- X
  gridInfoTrain <- adaptThresholdRange(xTrain)
  covTest <- cov(xTest)*(n-1)/n
  for (j in 1:nGrid){
    sigmaTrain <- adaptThreshold(gridInfoTrain$cov,gridInfoTrain$theta,grid[j],soft)
    error[j] <- (norm(sigmaTrain-covTest, "F"))
  }
  lambda <- grid[which(error == min(error))][1]
  sigma <- adaptThreshold(gridInfo$cov,gridInfo$theta,lambda,soft)
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))
  return(list(sigma = sigma, corr = corr))
}

#----------------------------------------------------------------------------------------
#  Range of the tuning parameter
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#  Output:
#      A list structure contains:
#       upper ------ upper bound of tuning parameter
#       lower ------ lower bound of tuning parameter
#         cov ------ sample covariance of x
#       theta ------ sample variance of covariance
#----------------------------------------------------------------------------------------

adaptThresholdRange <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  cov <- cov(x)*(n-1)/n
  centered.x <- scale(x, scale = FALSE)
  theta <- (t(centered.x)^2)%*%(centered.x^2)/n - cov^2
  delta <- cov/(theta^0.5)
  delta <- abs(delta - diag(diag(delta)))
  upper <- max(delta)
  lower <- min(delta[which(delta != 0)])
  return(list(upper = upper, lower = lower, theta = theta, cov = cov))
}

#----------------------------------------------------------------------------------------
#  Apply adaptive thresholding to the sample covariance
#  Input:
#           cov ------ p x p covariance matrix
#         theta ------ p x p variance of covariance matrix
#        lambda ------ tuning parameter
#          soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#         sigma ------ p x p matrix, adaptive thresholding result
#----------------------------------------------------------------------------------------

adaptThreshold <- function(cov,theta,lambda,soft){
  covOffDiag <- cov - diag(diag(cov))
  thetaOffDiag <- theta - diag(diag(theta))
  sigmaTmp <- abs(covOffDiag) - lambda*thetaOffDiag^0.5
  sigmaTmp[which(sigmaTmp < 0)] <- 0
  if (soft == 1){
    sigma <- diag(diag(cov)) + sigmaTmp*sign(covOffDiag)
  }else{
    sigma <- cov
    sigma[which(sigmaTmp < 1e-10)] <- 0
    sigma <- sigma + diag(diag(cov))
  }
  return(sigma)
}