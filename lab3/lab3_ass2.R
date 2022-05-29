############################### ASSIGNMENT 2 #####################################
# Metropolis Random Walk for Poisson regression
# Consider the following Poisson regression model: yi|Beta iid ~ Poisson[exp(t(xi) * Beta)], i = 1,...,n,
# where yi is the count for the ith observation in the sample and xi is the p-dimensional vector with 
# covariate observations for the ith observation. Use the data set eBayNumberOfBidderData.dat. 
# This dataset contains observations from 1000 eBay auctions of coins. The response variable is nBids and records the number of bids
# in each auction. The remaining variables are features/covariates (x): Const (for the intercept), some binary 1/0 vars, 
# LogBook (logarithm of the book value of the auctioned coin according to expert sellers. Standardized)
# MinBidShare (ratio of the minimum selling price (starting price) to the book value. Standardized).


############################### 2a #####################################
# Obtain the maximum likelihood estimator of Beta in the Poisson regression model for the eBay data 
# [Hint: glm.R, don't forget that glm() adds its own intercept so don't input the covariate Const]. 
# Which covariates are significant?

library("mvtnorm")
library("MASS")

data = read.table("eBayNumberOfBidderData.dat", header=TRUE)
n = nrow(data)

# fitting a generalized linear model using our data, nBids depends on all parameter except Const (the intercept). 
# family = poisson to specify Poisson regression model
set.seed(12345)
glmModel <- glm(nBids ~ . - Const, data = data, family = poisson)
glmModel



############################### 2b #####################################
# Let's do a Bayesian analysis of the Poisson regression. Let the prior be Beta ~ N[0, 100 * (t(X)X)^-1],
# where X is the n × p covariate matrix. This is a commonly used prior, which is called Zellner's g-prior. 
# Assume first that the posterior density is approximately multivariate normal: 
# Beta|y ~ N(Beta-tilde, J_y^-1(Beta-tilde)) , where Beta-tilde is the posterior mode and 
# Jy^-1(Beta-tilde) is the negative Hessian at the posterior mode. 
# These can be obtained by numerical optimization (optim.R) exactly like you already did for the logistic regression in Lab 2 
# (but with the log posterior function replaced by the corresponding one for the Poisson model, which you have to code up.).



#Lambda = exp(xB) in the Po-distr
y = data[,1]
X = as.matrix(data[,2:ncol(data)])

mu0 = as.matrix(c(rep(0,9)))
sigma0_sqr = 100*solve(t(X)%*%X)

# function for calculating the log posterior for the Poisson mode. 
LogPostPoisson <- function(betas,y,X,mu,sigma){
  #linPred <- X%*%betas;
  # lambda = exp(x*Beta) is our theta, i.e. the parameter we are interested in estimating (see L2 Poisson LH)
  # could also call this variable theta to be consistent with course notations
  lambda = as.matrix(exp(X%*%betas))
  logLikelihood = sum(y*log(lambda) - lambda - log(factorial(y)))
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, steer the optimizer away from here!
  logPrior <- dmvnorm(t(betas), t(mu), sigma, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates
  return(logLikelihood + logPrior)
}

# initialize Beta to be 9 zeros, since 9 covariates/variables
initVal <- matrix(0,9,1) 

# numeric optimization for Beta using above function
set.seed(12345)
OptimRes <- optim(initVal,LogPostPoisson,gr=NULL,y,X,mu0,sigma0_sqr,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# posterior mode, essentially what is optimized by the optimization. 
# posterior mode is same as posterior mean and posterior variance for multivariate normal distribution
betaTilde = OptimRes$par

# J_y_inv(betaTilde)
Jy_inv_betaTilde = -solve(OptimRes$hessian)
sigma = -solve(OptimRes$hessian)

betaTilde
glmModel$coefficients
Jy_inv_betaTilde



############################### 2c #####################################
# Let's simulate from the actual posterior of Beta using the Metropolis algorithm and compare the results with 
# the approximate results in b). Program a general function that uses the Metropolis algorithm to generate random draws from an
# arbitrary posterior density. In order to show that it is a general function for any model, 
# we denote the vector of model parameters by theta. 
# Let the proposal density be the multivariate normal density mentioned in Lecture 8 (random walk Metropolis):
# theta_p|theta^(i-1) ~ N(theta^(i-1), c * sigma) 
# i.e. above is: theta_proposal | theta previous ~ N()... 
#, where sigma = Jy^-1(beta-Tilde) wasobtained in b). 
# The value c is a tuning parameter and should be an input to your Metropolis function. 
# The user of your Metropolis function should be able to supply her own posterior density function, 
# not necessarily for the Poisson regression, and still be able to use your Metropolis function. 
# This is not so straightforward, unless you have come across function objects in R. 
# The note HowToCodeRWM.pdf in Lisam describes how yo can do this in R.
# Now, use your new Metropolis function to sample from the posterior of Beta in the Poisson regression for the eBay dataset. 
# Assess MCMC convergence by graphical methods.


# RWM sampler function
# prevBeta could also be called prevTheta to be more consistent with lecture notation
RWMSampler = function(c, sigma, prevBeta, logPostFunc, ...) {

  # sample proposal beta
  proposalBeta = rmvnorm(n=1, mean = prevBeta, sigma = c*sigma) 
  proposalBeta = t(proposalBeta)

  # compute acceptance probability, using the user selected function 
  # alpha = min (1, p(theta-p|y)/(p(theta-prev|y))) , = 1, p(y|theta-p)*p(theta-p) / p(y|theta-prev)*p(theta-prev)
  alpha = min(1, exp(logPostFunc(proposalBeta, ...) - logPostFunc(prevBeta, ...)) )
  
  # with probability alpha, set theta(i) = sample proposal, otherwise set theta(i) = theta(i-1) (i.e. prev beta)
  prob = runif(1)
  if (alpha >= prob) {
    return (proposalBeta)
  }
  return (prevBeta)
}

# initial beta0 (i.e. theta0 in the general case general)
initVal <- matrix(0,9,1) 
c = 3
# sigma was generated in b) 

nDraws = 10000
posteriorDraws = matrix(nrow=9, ncol=nDraws)
# For the first draw, use initVal
set.seed(12345)
posteriorDraws[,1] = RWMSampler(c, sigma, initVal, LogPostPoisson, y, X, mu0, sigma0_sqr)
# for draw 2 - nDraws
set.seed(12345)
for (i in 2:(nDraws)) {
  posteriorDraws[,i] = RWMSampler(c, sigma, (posteriorDraws[,i-1]), LogPostPoisson, y, X, mu0, sigma0_sqr)
}

# compare with approximation from b) and plot deviance
approximation_b = betaTilde
par(mfrow=c(3,3))
for (betaIndex in 1:9) {
  plot(c(1:nDraws), abs(posteriorDraws[betaIndex,]-glmEstimation[betaIndex]), type="l", col="red", 
  xlab = "draw", ylab="deviance true and approximated", 
  main = paste("Absolute deviance between actual and approximated Beta", betaIndex))
}



############################### 2d #####################################
# Use the MCMC draws from c) to simulate from the predictive distribution of
# the number of bidders in a new auction with the characteristics below. 
# Plot the predictive distribution. What is the probability of no bidders in this new auction?

# dont forget to add 1 for the intercept
x= c(1, 1, 0, 1, 0, 1, 0, 1.2, 0.8)

# for each of our posterior draws, perform a draw from the poisson distribution with lambda = exp(x*Beta)
# this draw represents the estimated number of bids
set.seed(12345)
nrBids = rpois(nDraws, lambda=exp(x%*%posteriorDraws))

#nrBids = c(1:nDraws)
#for (i in 1:nDraws) {
#  nrBids[i] = rpois(1, lambda=exp(x%*%posteriorDraws[,i]))
#}

par(mfrow=c(1,1))
hist(nrBids)

#count the probability the the nr of bids is 0
probNoBids = sum(nrBids==0)/nDraws
probNoBids


