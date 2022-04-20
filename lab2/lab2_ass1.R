library("mvtnorm")
library(mvtnorm)


data = read.table("TempLambohov.txt", header=TRUE)
n = nrow(data)

mju0 = t(c(-10,100,-100)) #-8 1st
omega0 = 0.05*diag(3) #0.05
omega0_inv = solve(0.05*diag(3)) #0.05
v0 = 10 #300
sigma0 = 0.1


nDraws = 1000
sigmaDraws = c(1:nDraws)
set.seed(12345)
#drawing 10000 sigma^2 from the inv-chi^2 distribution
sigmaDraws = ((v0)*sigma0)/rchisq(nDraws,v0) 

#drawing Betas from the mutivariate normal distribution
betaDraws = matrix(ncol=3, nrow=nDraws)
for (i in 1:nDraws) {
  betaDraws[i, ] = rmvnorm(n=1, mean=mju0, sigma=sigmaDraws[i]*omega0_inv)
}

plot(NULL, xlim=c(0,1), ylim=c(-20,30), ylab="temp", xlab="time")
for (i in 1:nDraws) {
  time = data$time
  beta0 = betaDraws[i, 1]
  beta1 = betaDraws[i, 2]
  beta2 = betaDraws[i, 3]
  regression = beta0 + beta1*time + beta2*time^2
  lines(time, regression, col=rgb(0,0,0,0.2))
  
}

lines(data$time, data$temp, col="red")


################1b
X = matrix(ncol = 3, nrow = nrow(data))
X[,1] = 1
X[,2] = data$time
X[,3] = (data$time)^2
y = data$temp

#X är: 1(för intercept), time, time^2, för alla 365 , dvs 3 col x 365 row
v_n = v0 + n
omega_n = t(X) %*% X + omega0

beta_hat = solve(t(X)%*%X) %*% t(X)%*%y
mju_n = solve(t(X)%*%X + omega0) %*% (t(X)%*%X%*%beta_hat + omega0%*%t(mju0))
sigma_n = ( v0*sigma0 + ( t(y %*% y) + mju0%*%omega0%*%t(mju0) - t(mju_n)%*%omega_n%*%mju_n ) ) / v_n


nDraws = 1000
postSigmaDraws = c(1:nDraws)
set.seed(12345)
#drawing 10000 sigma_n^2 from the posterior inv-chi^2 distribution
postSigmaDraws = ((v_n)*sigma_n)/rchisq(nDraws,v_n) 

#drawing Betas from the posterior mutivariate normal distribution
postBetaDraws = matrix(ncol=3, nrow=nDraws)
for (i in 1:nDraws) {
  postBetaDraws[i, ] = rmvnorm(n=1, mean=mju_n, sigma=postSigmaDraws[i]*solve(omega_n))
}

hist(postBetaDraws[,1])
hist(postBetaDraws[,2])
hist(postBetaDraws[,3])
hist(postSigmaDraws)

plot(data$time, data$temp, cex = 0.8)
temps = median(postBetaDraws[,1]) +  median(postBetaDraws[,2]) * data$time + median(postBetaDraws[,3]) * data$time^2
points(data$time, temps, type="l", col = "red")

#post_temps = postBetaDraws[,1] +  postBetaDraws[,2] * data$time + postBetaDraws[,3] * data$time^2

post_temps = matrix(ncol=nrow(data), nrow=nDraws)
for (i in 1:nDraws) {
  time = data$time
  beta0 = postBetaDraws[i, 1]
  beta1 = postBetaDraws[i, 2]
  beta2 = postBetaDraws[i, 3]
  regression = beta0 + beta1*time + beta2*time^2
  post_temps[i, ] = regression
  #lines(time, regression, col=rgb(0,0,0,0.2))
}

#for the confidence interval, gives us lower and uppwer bound for each of the 365 times
#apply "loops" over the vector
interval = apply(X=post_temps, MARGIN=2, FUN=function(x) quantile(x, probs=c(0.025,0.975)))

#lines(data$time, temps, type="l", col = "red")
lines(data$time, interval[1,], type="l")
lines(data$time, interval[2,], type="l")





##############################1c




