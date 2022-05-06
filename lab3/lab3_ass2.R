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
  print(dim(X))
  print(dim(betas))
  #linPred <- X%*%betas;
  lambda = as.matrix(exp(X%*%betas)) #theta samma som lambda
 # logLikelihood <- sum( linPred*y - log(1 + exp(linPred)) ); #eftersom log, skiljer sig från slides? 
  #logLikelihood = sum(y) * linPred - linPred*n
  
  logLikelihood = sum(y*log(lambda) - lambda - log(factorial(y)))
  
  
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, steer the optimizer away from here!
  logPrior <- dmvnorm(betas, mu0, sigma0, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates
  return(logLikelihood + logPrior)
}

#logPrior <- dmvnorm(t(initVal), mu0, sigma0, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates

initVal <- matrix(0,9,1) #default just 9 zeros, if 9 covariates/variables
OptimRes <- optim(initVal,LogPostPoisson,gr=NULL,y,X,mu0,sigma0,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
# posterior mode, essentially what is optimized by the optimization. 
# posterior mode is same as posterior mean and posteriro variance for multivariate normal distribution
betaTilde = OptimRes$par
betaTilde
# J_y_inv(betaTilde)
sigma = -solve(OptimRes$hessian)
sigma
betaTilde
glmModel$coefficients

#######c
RWMSampler = function(c, sigma, prevBeta, logPostFunc, ...) {
  #sample proposal
  proposalBeta = rmvnorm(n=1, mean = prevBeta, sigma = c*sigma) #rmvnorm - Draw from MultiVariate Normal Distribution 
  proposalBeta = t(proposalBeta)
  #compute acceptance probability, using user selected function 
  alpha = min(1, exp(logPostFunc(proposalBeta, ...)-logPostFunc(prevBeta, ...)) )
  #with probability alpha, set theta(i)=sample proposal, otherwise theta(i) = theta(i-1) (prev beta)
  prob = runif(1)
  if (alpha >= prob) {
    return (proposalBeta)
  }
  return (prevBeta)
}


# initialize theta0 
# initial beta0 
initVal <- matrix(0,9,1) #default just 9 zeros, if 9 covariates/variables
c = 3
# sigma from b)
# logPostFunc = LogPostPoisson

test = RWMSampler(c, sigma, initVal, LogPostPoisson, y, X, mu, sigma)



