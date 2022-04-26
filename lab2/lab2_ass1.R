library("mvtnorm")
library(mvtnorm)

data = read.table("TempLambohov.txt", header=TRUE)
n = nrow(data)

#parameters for multivarate distribution
mju0 = t(c(-10,100,-100)) #the means for the values we get from the multivariate normal distribution.  första -8 
omega0 = 0.05*diag(3) #omega * drawn sigma kind of gives us variance? 
omega0_inv = solve(0.05*diag(3)) 

#parameters for sigma distribution
v0 = 10 #parameter selected based on prior beliefs, related to degrees of freedom. #300 ist?
sigma0 = 0.1 #variance for the distribution of sigma.  



# simulate draws from joint prior, for every draw compute the regression curve. 
# i.e. we get one regression curve for each draw from the prior. 
nDraws = 1000
sigma2Draws = c(1:nDraws)
set.seed(12345)
#drawing 1000 sigma^2 from the inv-chi^2 distribution
sigma2Draws = ((v0)*sigma0)/rchisq(nDraws,v0) #in last lab, we used n, cause vn = v0 + n, and if v0 = 0 we get n in the posterior. 

#drawing Betas from the mutivariate normal distribution, using above sigma^2 draws
betaDraws = matrix(ncol=3, nrow=nDraws)
for (i in 1:nDraws) {
  betaDraws[i, ] = rmvnorm(n=1, mean=mju0, sigma=sigma2Draws[i]*omega0_inv) #mvnorm - MultiVariate Normal 
}

plot(NULL, xlim=c(0,1), ylim=c(-20,30), ylab="temp", xlab="time")
for (i in 1:nDraws) {
  time = data$time
  beta0 = betaDraws[i, 1]
  beta1 = betaDraws[i, 2]
  beta2 = betaDraws[i, 3]
  regressionLine = beta0 + beta1*time + beta2*time^2
  lines(time, regressionLine, col=rgb(0,0,0,0.15))
  
  #for visualizing how well the regression lines seem to fit
  if (i == nDraws) {
    lines(time, (mean(betaDraws[,1]) + mean(betaDraws[,2])*time + mean(betaDraws[,3])*time^2 + rnorm(1)), col="red" )
  }
  
}

lines(data$time, data$temp, col="green")


################1b
#X is the covariates matrix matrix, with columns representing the values for 1(for intercept), time, time^2, for all 365 times, i.e. 3 col x 365 row
X = matrix(ncol = 3, nrow = nrow(data))
X[,1] = 1
X[,2] = data$time
X[,3] = (data$time)^2
y = data$temp

#calculating posterior distribution variable values
v_n = v0 + n
omega_n = t(X) %*% X + omega0

beta_hat = solve(t(X)%*%X) %*% t(X)%*%y
mju_n = solve(t(X)%*%X + omega0) %*% (t(X)%*%X%*%beta_hat + omega0%*%t(mju0))
sigma_n = ( v0*sigma0 + ( t(y %*% y) + mju0%*%omega0%*%t(mju0) - t(mju_n)%*%omega_n%*%mju_n ) ) / v_n


nDraws = 1000
postSigmaDraws = c(1:nDraws)
set.seed(12345)
#drawing 1000 sigma_n^2 from the posterior inv-chi^2 distribution
postSigmaDraws = ((v_n)*sigma_n)/rchisq(nDraws,v_n) 

#drawing Betas from the posterior mutivariate normal distribution, using above sigma_n^2
postBetaDraws = matrix(ncol=3, nrow=nDraws)
for (i in 1:nDraws) {
  postBetaDraws[i, ] = rmvnorm(n=1, mean=mju_n, sigma=postSigmaDraws[i]*solve(omega_n))
}

hist(postBetaDraws[,1])
hist(postBetaDraws[,2])
hist(postBetaDraws[,3])
hist(postSigmaDraws)

#scatter plot of the data, with median posterior 
plot(data$time, data$temp, cex = 0.8)
temps = median(postBetaDraws[,1]) +  median(postBetaDraws[,2]) * data$time + median(postBetaDraws[,3]) * data$time^2
points(data$time, temps, type="l", col = "red")

#1000 (rows) x 365 (columns) matrix. The columns correspond to the 365 different times,
#and each row corresponds to the temperature values for the times, for the specific posterior beta draw. 
post_temps = matrix(ncol=nrow(data), nrow=nDraws)
for (i in 1:nDraws) {
  time = data$time
  beta0 = postBetaDraws[i, 1]
  beta1 = postBetaDraws[i, 2]
  beta2 = postBetaDraws[i, 3]
  regression = beta0 + beta1*time + beta2*time^2
  post_temps[i, ] = regression
  #for visualization, the confidence interval seems to be correct
  lines(time, regression, col=rgb(0,0,0,0.2))
}

#for the confidence interval, gives us lower and upper bound for each of the 365 times
#apply "loops" over the vector, thus we get a 2 (row) x 365 (column) matrix, where in this instance each column
#corresponds to the lower and uppwer bound for a time. So, row 1 corresponds to all lower bounds, and row 2 to all upper bounds. 
interval = apply(X=post_temps, MARGIN=2, FUN=function(x) quantile(x, probs=c(0.025,0.975)))

lines(data$time, interval[1,], type="l", col = "green")
lines(data$time, interval[2,], type="l", col = "green")
lines(data$time, temps, type="l", col = "red")





##############################1c

# temp(time) = b0 + b1*time + b2*time^2
# d(temp) /d(time) = b1 + 2*b2*time = 0 <=> time = - b1 / 2*b2

#times_highest_temp = c(1:nDraws)
times_highest_temp = -postBetaDraws[,2] / (2*postBetaDraws[,3]) 

plot(times_highest_temp)
hist(times_highest_temp, breaks = 20) #0.540 - 0.542 -> 0.541, i.e. 197 days (17 july)

#highest_temps = c(1:nDraws)
highest_temps = postBetaDraws[,1] + postBetaDraws[,2]*times_highest_temp + postBetaDraws[,3]*times_highest_temp^2
hist(highest_temps, breaks = 20)

plot(data$time, data$temp, cex = 0.5)
points(times_highest_temp, highest_temps, col = "blue")


##############################1d
#we can use spline regression, where we to avoid overfittign use a smoothness/shrinkage/regulariation
#prior. It can be selected using the same prior as previously, but setting mu0 = 0, sigma0 = lambda*I,
#where lambda determines the amount of shrinkage.


mu0 = 0
lambda = 1
sigma0 = lambda*diag(9)
#prior: B_j|sigma^2 ~N(0, sigma^2/lambda)