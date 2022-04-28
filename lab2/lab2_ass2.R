############################### ASSIGNMENT 2 #####################################
# Posterior approximation for classification with logistic regression
# The dataset WomenAtWork.dat contains n = 168 observations on the following eight variables related to women

############################### 2a #####################################
#Consider the logistic regression model: Pr(y = 1|x, beta) = exp(x*beta) /( 1 + exp(x*beta)),
# where y equals 1 if the woman works and 0 if she does not. 
# x is a 7-dimensional vector containing the seven features (including a 1 to model the intercept).
# The goal is to approximate the posterior distribution of the parameter vector
# beta with a multivariate normal distribution beta|y,x ~ N(betaTilde, J_y^-1(betaTilde)), 
# where betaTilde is the posterior mode and J(betaTilde) = - 2nd derivative lol, is the negative of the observed Hessian 
# evaluated at the posterior mode.
# Note that lol isa 7 × 7 matrix with second derivatives on the diagonal and cross-derivative lol2 on the off-diagonal. 
# You can compute this derivative by hand, but we will let the computer do it numerically for you. 
# Calculate both betaTilde and J(betaTilde) using the optim function in R. 
# Use the prior beta ??? N (0, tao^2 * I), where tao = 5.

# Present the numerical values of betaTilde and J(betaTilde) for the WomenAtWork data. 
# Compute an approximate 95% equal tail posterior probability interval for the regression coefficient to the variable
# NSmallChild. Would you say that this feature is of importance for the probability that a woman works?
# [Hint: You can verify that your estimation results are reasonable by comparing  the posterior means to the 
# maximum likelihood estimates, given by: glmModel <- glm(Work ~ 0 + ., data = WomenAtWork, family = binomial).]

library("mvtnorm")

data = read.table("WomenAtWork.dat", header=TRUE)
n = nrow(data)

y = data$Work
X = as.matrix(data[,2:8])
nr_parameters <- ncol(X)

# prior: Beta ~ N(mu0, tao^2 * I)
mu0 <- as.matrix(rep(0, nr_parameters))
tao = 5
sigma0 = tao^2*diag(nr_parameters)

# if we want to standardize 
standardize <- FALSE 
if (standardize){
  Index <- 2:7 #dont normalize intercept
  X[,Index] <- scale(X[,Index]) #thus only scale columns 2-7
}

LogPostLogistic <- function(betas,y,X,mu,sigma){
  linPred <- X%*%betas;
  logLikelihood <- sum( linPred*y - log(1 + exp(linPred)) ); #eftersom log, skiljer sig från slides? 
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, steer the optimizer away from here!
  logPrior <- dmvnorm(betas, mu0, sigma0, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates
  return(logLikelihood + logPrior)
}


# Select the initial values for beta
initVal <- matrix(0,nr_parameters,1) #default just 7 zeros, if 7 covariates/variables


# first parameter is initVal, which is what we are optimizing, using function LogPostLogistic,
# which also takes in the parameters y, X, mu0, sigma0. 
# y is response variable (working or not) 
# X are our variables/covariates we worked with before.
# mu0, sigma0 were given at the start. 
# inimize using BFGS. 
# The argument control is a list of options to the optimizer optim, where fnscale=-1 means that we minimize 
# the negative log posterior. Hence, we maximize the log posterior. 
# hessian=TRUE cause we want the hessian to be computed as well
OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu0,sigma0,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# posterior mode, essentially what is optimized by the optimization. 
# posterior mode is same as posterior mean and posteriro variance for multivariate normal distribution
betaTilde = OptimRes$par
betaTilde

# J_y_inv(betaTilde)
J_betaTilde = -solve(OptimRes$hessian)
J_betaTilde


# simulate 10000 draws from posterior
nDraws = 10000
betaDraws = rmvnorm(n=nDraws, mean=betaTilde, sigma=J_betaTilde)

# Computing approximate 95% equal tail posterior probability interval for the regression coefficient 
# to the variable NSmallChild, which is betaDraws[,6]
interval = quantile(as.matrix(betaDraws[,6]), probs=c(0.025,0.975))
lowerBound = interval[1]
upperBound = interval[2]
lowerBound
upperBound
# alternatively, we can compute the lower and upper bounds like this: 
# coeff_smallkids = betaTilde[6]
# lowerBound = coeff_smallkids - 1.96*sd(betaDraws[,6])
# upperBound = coeff_smallkids + 1.96*sd(betaDraws[,6])









###############3


plot(density(betaDraws[,6]))
abline(v=lowerBound, col="red")
abline(v=upperBound, col="red")

plot(density(betaDraws))
abline(v=lowerBound, col="red")
abline(v=upperBound, col="red")


glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)
glmModel


############################2b
x_new = as.matrix(c(1, 20, 12, 8, 43, 0, 2))
log_reg = exp(t(x_new)%*%t(betaDraws)) / (1 + exp(t(x_new)%*%t(betaDraws)))
plot(density(log_reg))



############################2c
prob_working = mean(log_reg)

#nr_working = c(1:nDraws)
#for (i in 1:nDraws) {
 # nr_working[i] = sum(rbinom(11, 1, prob_working))
#}
#hist(nr_working, freq = 11)

nr_working2 = rbinom(nDraws, 11, prob_working)
hist(nr_working2)