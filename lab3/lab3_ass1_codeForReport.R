library("mvtnorm")
library(mvtnorm)
# library(MASS) # To access the mvrnorm() function

###################### 1 A ###########################
# Gibbs sampler for a normal model
# The dataset Precipitation.rds consists of daily records of weather with rain or snow (in units of mm) from the beginning of
# 1962 to the end of 2008 in a certain area. 
# Assume the natural log of the daily precipitation {y1, ..., yn} to be independent normally distributed,
# ln y1, ...,ln yn|mu, sigma^2 ~ iid N(mu, sigma^2) where both mu and sigma^2 are unknown. 
# Let mu ~ N(mu0, tau0squared) independently of sigma^2 ~ Inv-chi^2(v0, sigma0^2)

# Implement (code!) a Gibbs sampler that simulates from the joint posterior p(mu, sigma^2 | ln y1,...,ln yn). 
# The full conditional posteriors are given on the slides from Lecture 7. 
# Evaluate the convergence of the Gibbs sampler by calculating the Ineffciency Factors (IFs) 
# and by plotting the trajectories of the sampled Markov chains.

data=readRDS("Precipitation.rds")
n = length(data)
log_data = log(data)

# Initial parameters
mu0 = mean(log_data)
sigma0 = 10
v0 = 250
tau0 = 100
nDraws = 10000

# Priors independent --> use Gibbs sampling to sample from the posterior 
# posterior distribution of mu and sigma

# Want to sample from our two full conditional posteriors (mu and sigma respectively)

# Gibbs sampling
gibbsDraws <- matrix(0,nDraws,2)
colnames(gibbsDraws) = c("mu", "sigma")
sigma <- 1 # Initial value for sigma
v_n = v0+n # Degrees of freedom

set.seed(12345)
for (i in 1:nDraws){
  
  # Calculate posterior parameters
  w = (n/sigma) / ((n/sigma)+1/tau0)
  mu_n = w*mean(log_data) + (1-w)*mu0
  tau_n = 1/((n/sigma)+(1/tau0))
  
  # Update mu given sigma
  mu <- rnorm(1, mean = mu_n, sd = tau_n)
  gibbsDraws[i,1] <- mu
  
  # Update sigma given mu
  parameter_squared = (v0*sigma0+sum(log_data-mu)^2)/(n+v0)
  sigma = v_n*parameter_squared/rchisq(1,v_n)
  gibbsDraws[i,2] <- sigma
}

# Calculate Inefficiency Factors (IFs)
a_Gibbs = acf(gibbsDraws[,1])
IF_Gibbs = 1+2*sum(a_Gibbs$acf[-1])

#Plot Gibbs sampling
plot(1:nDraws, gibbsDraws[,1], type="l", col="orange")
hist(gibbsDraws[,1], col="orange")

a_Gibbs = acf(gibbsDraws[,2])
IF_Gibbs = 1+2*sum(a_Gibbs$acf[-1])

plot(1:nDraws, gibbsDraws[,2], type="l", col="orange")
hist(gibbsDraws[,2], col="orange")


###################### 1 B ###########################
# Plot the following in one figure: 
# 1) a histogram or kernel density estimate of the daily precipitation {y1, ..., yn}. 
# 2) The resulting posterior predictive density p(ytilde|y1,...,yn) using the simulated posterior draws from (a)
# How well does the posterior predictive density agree with this data?

# Want to use our sigma and mu from the Gibbs sampling to make draws of y

postY = c(1:nDraws)

set.seed(12345)
postY = rnorm(n=nDraws, mean=exp(gibbsDraws[,1]), sd=exp(gibbsDraws[,2]))

plot(density(data), lwd=2, axes=FALSE)
axis(1)
axis(2)
lines(density(postY), col="orange", lwd=2)

