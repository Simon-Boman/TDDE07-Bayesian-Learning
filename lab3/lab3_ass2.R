library("mvtnorm")
library("MASS")


data = read.table("eBayNumberOfBidderData.dat", header=TRUE)
n = nrow(data)

#########################a

glmModel <- glm(nBids ~ . - Const, data = data, family = poisson)
glmModel
#minBidShare, sealed, verifyID



#########b

y = data[,1]
X = as.matrix(data[,2:ncol(data)])

mu0 = as.matrix(c(rep(0,9)))
sigma0 = 100*solve(t(X)%*%X)


#exp(x*betas) är vårt theta i slidsen L2 poisson likelihood.
LogPostPoisson <- function(betas,y,X,mu,sigma){
  linPred <- X%*%betas;
  lambda = as.matrix(exp(X%*%betas)) #theta samma som lambda
 # logLikelihood <- sum( linPred*y - log(1 + exp(linPred)) ); #eftersom log, skiljer sig från slides? 
  #logLikelihood = sum(y) * linPred - linPred*n
  
  logLikelihood = sum(y*log(lambda) - lambda - log(factorial(y)))
  
  
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, steer the optimizer away from here!
  logPrior <- dmvnorm(betas, mu0, sigma0, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates
  return(logLikelihood + logPrior)
}

#logPrior <- dmvnorm(t(initVal), mu0, sigma0, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates

initVal <- matrix(0,9,1) #default just 7 zeros, if 7 covariates/variables



OptimRes <- optim(initVal,LogPostPoisson,gr=NULL,y,X,mu0,sigma0,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# posterior mode, essentially what is optimized by the optimization. 
# posterior mode is same as posterior mean and posteriro variance for multivariate normal distribution
betaTilde = OptimRes$par
betaTilde

# J_y_inv(betaTilde)
J_betaTilde = -solve(OptimRes$hessian)
J_betaTilde

betaTilde
glmModel$coefficients

#######c







