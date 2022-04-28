############################### ASSIGNMENT 1 #####################################
# The dataset TempLambohov.txt contains daily average temperatures (in degree Celcius) at Lambohov, Linköping over the course of the year 2019. 
# The response variable is temp and the covariate is time = the number of days since the beginning of the year/365
# Quadratic regression model: temp = Beta0 + Beta1 · time + Beta2 · time2 + epsilon, epsilon iid ~ N(0, sigma^2)

############################### 1a #####################################
# Use the conjugate prior for the linear regression model. The prior hyperparameters mu0, omega0, ??0 and sigma0^2 shall be set to sensible values. 
# Check if this prior agrees with your prior opinions by simulating draws from the joint prior of all parameters and for every draw compute 
# the regression curve. This gives a collection of regression curves; one for each draw from the prior. Does the collection of curves look reasonable?
# If not, change the prior hyperparameters until the collection of prior regression curves agrees with your prior beliefs about the regression curve.
# [Hint: R package mvtnorm can be used and your Inv-chi^2 simulator of random draws from Lab 1.]

library("mvtnorm")
library(mvtnorm)

data = read.table("TempLambohov.txt", header=TRUE)
n = nrow(data)


# parameters for the Multivariate Normal Distribution (N_3)
mu0 = c(-10,100,-100) # the means for the 3 values we get from the multivariate normal distribution.
omega0 = 0.05*diag(3) # omega0 * drawn sigma (sigma^2) kind of gives us variance (and covariance) for the multivariate normal distr? 

# parameters for sigma distribution
v0 = 10 # parameter selected based on prior beliefs, related to degrees of freedom.
sigma0 = 0.1 # variance for the distribution of sigma.  



# simulate draws from joint prior of all parameters, and for every draw compute the regression curve. 
# i.e. we get one regression curve for EACH draw from the prior. 
nDraws = 1000
set.seed(12345)
# drawing 1000 sigma^2 from the inv-chi^2 distribution (prior)
sigma2Draws = c(1:nDraws)
sigma2Draws = ((v0)*sigma0)/rchisq(nDraws,v0) #in last lab, we used n (instead of v0), cause vn = v0 + n, and if v0 = 0 we get n in the posterior. 

# drawing Betas from the mutivariate normal distribution (prior), using above sigma^2 draws
betaDraws = matrix(ncol=3, nrow=nDraws)
for (i in 1:nDraws) {
  betaDraws[i, ] = rmvnorm(n=1, mean=t(mu0), sigma=sigma2Draws[i]*solve(omega0)) #rmvnorm - Draw from MultiVariate Normal Distribution 
}

# for each draw of Betas, compute and plot the regression curve
plot(NULL, xlim=c(0,1), ylim=c(-20,30), ylab="temp", xlab="time", main="Plot of regression curve for each draw of Betas")
legend(x = 0.6, y = 31, legend = c("Regression curve for each draw of Betas", "Mean of the Beta draws", "True temperature"), col = c("black","red", "green"), cex = 0.6, lwd = 3)

for (i in 1:nDraws) {
  time = data$time
  beta0 = betaDraws[i, 1]
  beta1 = betaDraws[i, 2]
  beta2 = betaDraws[i, 3]
  regressionCurve = beta0 + beta1*time + beta2*time^2
  lines(time, regressionCurve, col=rgb(0,0,0,0.15))
  
  # for visualizing how well the regression lines seem to fit, plot mean curve 
  if (i == nDraws) {
    lines(time, (mean(betaDraws[,1]) + mean(betaDraws[,2])*time + mean(betaDraws[,3])*time^2 + rnorm(1)), col="red" )
  }
}

# plot the actual true temperatures
lines(data$time, data$temp, col="green")





############################### 1b #####################################
# Write a function that simulate draws from the joint posterior distribution of beta0, beta1, beta2, sigma^2.
# i. Plot a histogram for each marginal posterior of the parameters.
# ii. Make a scatter plot of the temperature data and overlay a curve for the posterior median of the regression function 
# f(time) = E [temp|time] = beta0 + beta1 * time + beta2 * time^2,
# i.e. the median of f (time) is computed for every value of time. 
# In addition, overlay curves for the 95% equal tail posterior probability intervals of f(time),
# i.e. the 2.5 and 97.5 posterior percentiles of f (time) is computed for every value of time. 
# Does the posterior probability intervals contain most of the data points? Should they?


# X is the covariates matrix, with columns representing the values for 
# the intercept (1), time, time^2, for all 365 times, i.e. 3 col x 365 row matrix. 
# y are the 365 temperatures
X = matrix(ncol = 3, nrow = nrow(data))
X[,1] = 1
X[,2] = data$time
X[,3] = (data$time)^2
y = data$temp

#calculating posterior distribution parameters
v_n = v0 + n
omega_n = t(X) %*% X + omega0

beta_hat = solve(t(X)%*%X) %*% t(X)%*%y
mju_n = solve(t(X)%*%X + omega0) %*% (t(X)%*%X%*%beta_hat + omega0%*%mu0)
sigma_n = ( v0*sigma0 + ( t(y %*% y) + mu0%*%omega0%*%mu0 - t(mju_n)%*%omega_n%*%mju_n ) ) / v_n


# simulating draws from the posterior
nDraws = 1000
set.seed(12345)
# drawing sigma_n^2 from the posterior inv-chi^2 distribution
postSigmaDraws = c(1:nDraws)
postSigmaDraws = ((v_n)*sigma_n)/rchisq(nDraws,v_n) 

# drawing Betas from the posterior multivariate normal distribution, using above sigma_n^2
postBetaDraws = matrix(ncol=3, nrow=nDraws)
for (i in 1:nDraws) {
  postBetaDraws[i, ] = rmvnorm(n=1, mean=mju_n, sigma=postSigmaDraws[i]*solve(omega_n))
}

# i. plot histograms of each marginal posterior of the parameters 
hist(postBetaDraws[,1])
hist(postBetaDraws[,2])
hist(postBetaDraws[,3])
hist(postSigmaDraws)

# ii. scatter plot of the input data, and a curve for the median posterior regression curve using our posterior draws
plot(data$time, data$temp, cex = 0.8, ylab="temp", xlab="time", main="Regression curve for each posterior draw of Betas, median posterior curve, and posterior probability interval")
temps = median(postBetaDraws[,1]) +  median(postBetaDraws[,2]) * data$time + median(postBetaDraws[,3]) * data$time^2
points(data$time, temps, type="l", col = "red")
legend(x = 0.3, y = -9, legend = c("Regression curve for each posterior draw of Betas", "Posterior median of the regression curve for Beta draws", "Posterior probability interval"), col = c("black","red", "green"), cex = 0.6, lwd = 3)

# post_temps is a 1000 (rows) x 365 (columns) matrix. The columns correspond to the 365 different times,
# and each row corresponds to the temperature values for the times, for the different posterior beta draws. 
post_temps = matrix(ncol=nrow(data), nrow=nDraws)
for (i in 1:nDraws) {
  time = data$time
  beta0 = postBetaDraws[i, 1]
  beta1 = postBetaDraws[i, 2]
  beta2 = postBetaDraws[i, 3]
  regressionCurve = beta0 + beta1*time + beta2*time^2
  post_temps[i, ] = regressionCurve
  #for visualizing if the confidence interval seems to be correct
  lines(time, regressionCurve, col=rgb(0,0,0,0.2))
}

# for the confidence interval, calculate the lower and upper bound for each of the 365 times
# and then plot the lower and upper interval curves

# interval matrix to store the lower and upper bound for each of the 365 times
# for each of the 365 times, calculate lower and upper bound using the quantile function and our post_temps matrix
# which for each of the 365 times contains different values from our beta draws
# e.g. for col 1, where time = 0, calculate the posterior probability interval using our 1000 beta draws and their values 
# for time= 0, and repeat this for each of the 365 times. 
interval = matrix(nrow = 2, ncol = 365)
for (i in 1:nrow(data)) {
  interval[,i] = quantile(post_temps[,i], probs=c(0.025, 0.975)) 
}
# apply "loops" over the vector, thus we get a 2 (row) x 365 (column) matrix, where in this instance each column
# corresponds to the lower and upper bound for a time. So, row 1 corresponds to all lower bounds, and row 2 to all upper bounds. 
# interval = apply(X=post_temps, MARGIN=2, FUN=function(x) quantile(x, probs=c(0.025,0.975)))

# plot the equal tailposterior probability intervals
lines(data$time, interval[1,], type="l", col = "green")
lines(data$time, interval[2,], type="l", col = "green")
# plot the original data for visualization
lines(data$time, temps, type="l", col = "red")






############################### 1c #####################################
# It is of interest to locate the time with the highest expected temperature (i.e. the time where f(time) is maximal). 
# Let's call this value xtilde. Use the simulated draws in (b) to simulate from the posterior distribution of xtilde. 
# [Hint: the regression curve is a quadratic polynomial. Given each posterior draw of beta10, beta1, beta2,
# you can find a simple formula for xtilde.]


# we have the following quadratic polynomial regression curve:
# temp(time) = b0 + b1*time + b2*time^2
# take the derivative of the response variable (temp) w.r.t the covariate (time):
# d(temp) /d(time) = b1 + 2*b2*time
# setting this to 0 and bryt ut time gives us the following expression for time: 
# b1 + 2*b2*time = 0 <=> time = - b1 / 2*b2

# calculate the times which correspond to highest measured temperature
times_highest_temp = -postBetaDraws[,2] / (2*postBetaDraws[,3]) 

plot(density(times_highest_temp))
print(density(times_highest_temp)) #
hist(times_highest_temp, breaks = 20) 
# The time with most highest temperature recording has mean and median 0.5448,
# i.e. day 198 (18 july), which seems reasonable.

# the highest temperature values
highest_temps = postBetaDraws[,1] + postBetaDraws[,2]*times_highest_temp + postBetaDraws[,3]*times_highest_temp^2
hist(highest_temps, breaks = 20)

# the times which had highest temperatures, and the temperatures
plot(data$time, data$temp, cex = 0.5, ylab="highest temperatures", xlab="time for highest temperatures")
points(times_highest_temp, highest_temps, col = "blue")
# the points are very clustered, which makes sense since generally the highest temperature and the time
# for the highest temperature should be close to the same for our posterior draws.






############################### 1d #####################################
# Say now that you want to estimate a polynomial regression of order 8,
# but you suspect that higher order terms may not be needed, and you worry about overfitting the data. 
# Suggest a suitable prior that mitigates this potential problem. You do not need to compute the posterior. 
# Just write down your prior. [Hint: the task is to specify mu0 and omega0 in a suitable way.]



# We can use spline regression to solve this problem, and to avoid overfitting we can use a 
# smoothness/shrinkage/regularization prior. It can be selected using the same prior as previously, 
# but setting mu0 = 0, and sigma0 = lambda*I, where lambda determines the amount of shrinkage.
# a larger lambda means more shrinkage/regularization. 
# mu = 0 results in small values of the betas, and a large lambda results in small variance (since in the prior
# we have sigma^2/lambda), which results in in many values being around 0.
# thus, using these values gives us shrinkage. 

#beta0 (intercept), beta1, ... , beta8
mu0 = c(rep(0,9))
lambda = 10
sigma0 = lambda*diag(9)
#prior: B_j|sigma^2 ~ N(0, sigma^2/lambda)

