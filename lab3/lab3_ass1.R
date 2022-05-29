############################### ASSIGNMENT 1 #####################################
# Gibbs sampler for a normal model
# The dataset Precipitation.rds consists of daily records of weather with rain or snow (in units of mm) from the beginning of
# 1962 to the end of 2008 in a certain area. 
# Assume the natural log of the daily precipitation {y1, ..., yn} to be independent normally distributed,
# ln y1, ...,ln yn|mu, sigma^2 ~ iid N(mu, sigma^2) where both mu and sigma^2 are unknown. 
# Let mu ~ N(mu0, tau0squared) independently of sigma^2 ~ Inv-chi^2(v0, sigma0_sq^2)


############################### 1a #####################################
# Implement (code!) a Gibbs sampler that simulates from the joint posterior p(mu, sigma^2 | ln y1,...,ln yn). 
# The full conditional posteriors are given on the slides from Lecture 7. 
# Evaluate the convergence of the Gibbs sampler by calculating the Inefficiency Factors (IFs) 
# and by plotting the trajectories of the sampled Markov chains.

library("mvtnorm")
library(mvtnorm)

data_original = readRDS("Precipitation.rds")
n = length(data)
data = log(data_original)

mu0 = 15 #expected rainfall
tao0_sq = 180 #variance of expected rainfaill
v0 = 140 #expected value of variance
sigma0_sq = 17 #expected 

nDraws <- 10000
gibbsDraws <- matrix(0,nDraws,2)
colnames(gibbsDraws) = c("mu", "sigma")

vn = v0 + n 

sigma_sq <- 1 # Initial value for sigma for the first sampling 
set.seed(12345)
for (i in 1:nDraws){  
  
  #update tao_n^2, w, mu_n
  #alt dela med tao0^2 om vi anger den som ej squared istället.
  tao_n_sq = 1/(n/sigma_sq + 1/tao0_sq)  
  w = (n/sigma_sq) / (n/sigma_sq + 1/tao0_sq) 
  mu_n = w*mean(data) + (1-w)*mu0
  
  # Update mu given sigma
  mu = rnorm(n = 1, mean = mu_n, sd = tao_n_sq)
  gibbsDraws[i,1] <- mu
  
  # Update sigma given mu
  # tao_param2_sq = variance? 
  tao_param2_sq = (v0*sigma0_sq + sum(data-mu)^2) / (n + v0)
  sigma_sq = (vn*tao_param2_sq/rchisq(1, vn))
  gibbsDraws[i,2] <- sigma_sq
  
#vi har X^2 (vn, param2)
      # vn*param2 + 
  #     X^2 (v0, sigma0_sq^2 )
  #    ((v0)*sigma0_sq)/rchisq(nDraws,v0) 
}

# acf - autocorrelation
a_Gibbs_mu <- acf(gibbsDraws[,1])
IF_Gibbs_mu <- 1+2*sum(a_Gibbs_mu$acf[-1])
IF_Gibbs_mu

a_Gibbs_sigma <- acf(gibbsDraws[,2])
IF_Gibbs_sigma <- 1+2*sum(a_Gibbs_sigma$acf[-1])
IF_Gibbs_sigma


plot(1:nDraws, gibbsDraws[,1], type = "l",col="red", main = "Markov chain trajectory, mu") # traceplot of Gibbs draws
hist(gibbsDraws[,1],col="red") # histogram of Gibbs draws
#cusumData =  cumsum(gibbsDraws[,1])/seq(1,nDraws) # Cumulative mean value of mu, Gibbs draws
#plot(1:nDraws, cusumData, type = "l", col="red")
#barplot(height = a_Gibbs_mu$acf[-1],col="red") # acf for Gibbs draws

plot(1:nDraws, gibbsDraws[,2], type = "l",col="red", main = "Markov chain trajectory, sigma") # traceplot of Gibbs draws
hist(gibbsDraws[,2],col="red") # histogram of Gibbs draws
#cusumData =  cumsum(gibbsDraws[,2])/seq(1,nDraws) # Cumulative mean value of mu, Gibbs draws
#plot(1:nDraws, cusumData, type = "l", col="red")
#barplot(height = a_Gibbs_sigma$acf[-1],col="red") # acf for Gibbs draws

# Plotting the cumulative path of estimates of Pr(theta1>0, theta2>0)
# par(mfrow=c(2,1))
# plot(cumsum(directDraws[,1]>0 & directDraws[,2]>0)/seq(1,nDraws),type="l", main='Direct draws', xlab='Iteration number', ylab='', ylim = c(0,1))
# plot(cumsum(gibbsDraws[,1]>0 & gibbsDraws[,2]>0)/seq(1,nDraws),type="l", main='Gibbs draws', xlab='Iteration number', ylab='', ylim = c(0,1))



############################### 1b #####################################
# Plot the following in one figure: 
# 1) a histogram or kernel density estimate of the daily precipitation {y1, ..., yn}. 
# 2) The resulting posterior predictive density p(ytilde|y1,...,yn) using the simulated posteriro draws from (a)
# How well does the posterior predictive density agree with this data?
  
  
hist(data_original, main = "original data")
posteriorDraws = rnorm(n = nDraws, mean = gibbsDraws[,1], sd = gibbsDraws[,2])
#exp to convert back to normal scale
posteriorDraws = exp(posteriorDraws)
hist(posteriorDraws, main = "posterior draws")

plot(density(posteriorDraws), col = "red", xlim=c(0,50), main = "green is original data, red is draws from posterior")
lines(density(data_original), col = "green")
#lines(density(data), col = "blue")
