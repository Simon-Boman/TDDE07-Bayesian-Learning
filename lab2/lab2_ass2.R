library("mvtnorm")

data = read.table("WomenAtWork.dat", header=TRUE)
n = nrow(data)

### Prior and data inputs ###
Covs <- c(2:7) # Select which covariates/features to include
standardize <- TRUE # If TRUE, covariates/features are standardized to mean 0 and variance 1
lambda <- 1# scaling factor for the prior of beta 

y = data$Work
X = as.matrix(data[,2:8])
n_par <- ncol(X)


mu0 = 0
mu0 <- as.matrix(rep(0,n_par))
tao = 5
sigma0 = tao^2*diag(7)
#Beta ~ N(0, tao^2* I)


LogPostLogistic <- function(betas,y,X,mu,sigma){
  linPred <- X%*%betas;
  logLikelihood <- sum( linPred*y - log(1 + exp(linPred)) ); #eftersom log, skiljer sig från slides? 
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, steer the optimizer away from here!
  logPrior <- dmvnorm(betas, mu0, sigma0, log=TRUE); #dmvnorm - Multivariate Normal Density and Random Deviates
  return(logLikelihood + logPrior)
}


# Select the initial values for beta
initVal <- matrix(0,n_par,1) #default just 11 zeros, if 11 covs.

# The argument control is a list of options to the optimizer optim, where fnscale=-1 means that we minimize 
# the negative log posterior. Hence, we maximize the log posterior. 
# y = 1 if working. X = our variables we worked with before. 
OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu0,sigma0,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)




#posterior mode
betaTilde = OptimRes$par
betaTilde

#J_y_inv(betaTilde)
J_betaTilde = -solve(OptimRes$hessian)
J_betaTilde

coeff_smallkids = betaTilde[6]
coeff_smallkids

nDraws = 10000
betaDraws = rmvnorm(n=nDraws, mean=betaTilde, sigma=J_betaTilde)

interval = apply(X=as.matrix(betaDraws[,6]), MARGIN=2, FUN=function(x) quantile(x, probs=c(0.025,0.975)))

#plot(density(betaDraws[,6]))
plot(density(betaDraws))

abline(v=interval[1,], col="red")
abline(v=interval[2,] ,col="red")


glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)





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
hist(nr_working, freq = 11)

nr_working2 = rbinom(nDraws, 11, prob_working)
hist(nr_working2)