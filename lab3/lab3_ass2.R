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
 # print("x, betas")
  #print(dim(X))
  #print(dim(betas))
  #linPred <- X%*%betas;
  lambda = as.matrix(exp(X%*%betas)) #theta samma som lambda
 # logLikelihood <- sum( linPred*y - log(1 + exp(linPred)) ); #eftersom log, skiljer sig från slides? 
  #logLikelihood = sum(y) * linPred - linPred*n
  
  logLikelihood = sum(y*log(lambda) - lambda - log(factorial(y)))
  
  
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, steer the optimizer away from here!
#  print("betas, mu")
 # print(dim(betas))
  #print(dim(mu0))
  logPrior <- dmvnorm(t(betas), t(mu), sigma, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates
  #print("klart")
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
J_y_inv_betaTilde = -solve(OptimRes$hessian)
sigma = -solve(OptimRes$hessian)
sigma
betaTilde
glmModel$coefficients




#######c
RWMSampler = function(c, sigma, prevBeta, logPostFunc, ...) {
  #sample proposal
 # print(dim(prevBeta))
  #print(dim(c*sigma))
  proposalBeta = rmvnorm(n=1, mean = prevBeta, sigma = c*sigma) #rmvnorm - Draw from MultiVariate Normal Distribution 
  proposalBeta = t(proposalBeta)
  #print("DIM PROPOSAL:")
  #print(dim(proposalBeta))
  #compute acceptance probability, using user selected function 
  alpha = min(1, exp(logPostFunc(proposalBeta, ...)-logPostFunc(prevBeta, ...)) )
  #with probability alpha, set theta(i)=sample proposal, otherwise theta(i) = theta(i-1) (prev beta)
  prob = runif(1)
  #print("proposal, prev")
  #print(dim(proposalBeta))
  #print(dim(prevBeta))
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

#test = RWMSampler(c, sigma, initVal, LogPostPoisson, y, X, mu0, sigma0)
#test#OptimRes <- optim(initVal,LogPostPoisson,gr=NULL,y,X,mu0,sigma0,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

nDraws = 10000
posteriorDraws = matrix(nrow=9, ncol=nDraws)
posteriorDraws[,1] = RWMSampler(c, sigma, initVal, LogPostPoisson, y, X, mu0, sigma0)
posteriorDraws[,1]

#for (i in 2:nDraws-1) {
  #posteriorDraws[,i] = RWMSampler(c, sigma, (posteriorDraws[,i-1]), LogPostPoisson, y, X, mu0, sigma0)
 # posteriorDraws[,i] = RWMSampler(c, sigma, matrix(posteriorDraws[,i-1],9,1), LogPostPoisson, y, X, mu0, sigma0)
  #
#}



for(i in 1:(nDraws-1)){
  #print(dim(posteriorDraws))
  posteriorDraws[,i+1]= RWMSampler(c, sigma, posteriorDraws[,i], LogPostPoisson, y, X, mu0, sigma0)
  
}


glmEstimation = glmModel$coefficients
par(mfrow=c(3,3))
for (betaIndex in 1:9) {
  plot(c(1:nDraws), abs(posteriorDraws[betaIndex,]-glmEstimation[betaIndex]), type="l", col="red", xlab = "draw", ylab="deviance true and approximated", main = paste("Absolute deviance between actual and approximated Beta", betaIndex) )
  
}



##################################### 2d
#intercept and x
x= c(1, 1, 0, 1, 0, 1, 0, 1.2, 0.8)

nrBids = c(1:nDraws)
for (i in 1:nDraws) {
  nrBids[i] = rpois(1, lambda=exp(x%*%posteriorDraws[,i]))
}
nrBids
#nrBids = rpois(nDraws, lambda=t(x%*%posteriorDraws))
par(mfrow=c(1,1))
hist(nrBids)



probNoBids = sum(nrBids==0)/nDraws
probNoBids


