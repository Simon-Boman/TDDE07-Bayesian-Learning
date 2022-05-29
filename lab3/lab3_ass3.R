############################### ASSIGNMENT 3 #####################################
# Time series models in Stan

############################### 3a #####################################
# Write a function in R that simulates data from the AR(1)-process xt = mu + phi(xt-1 - mu) + epslion_t, epsilon_t ~N(0, sigma^2),
# for given values of mu, phi and sigma^2. Start the process at x1 = mu and then simulate values for xt for t = 2, 3,...,T 
# and return the vector x1:T containing all time points. Use my = 13, sigma^2 = 3 and T = 300 and look at some different realizations
# (simulations) of x1:T for values of phi between ???1 and 1 (this is the interval of phi where the AR(1)-process is stationary). 
# Include a plot of at least one realization in the report. What effect does the value of phi have on x1:T ?
  
library(rstan)

mu = 13
sigma2 = 3
T = 300

# function that simulate values for x2, x3, ..., xT given a phi
ARSim = function(phi) {
  x = c(1:T)
  x[1] = mu
  for (i in 2:T) {
    # the rnorm() term corresponds to epsilon.
    x[i] = mu + phi * (x[i-1] - mu) + rnorm(n=1, mean = 0, sd = sqrt(sigma2))
  }
  return(x)
}
  

# the different phi values to be used to get different simulations
phiVals = seq(-1, 1, 0.1)

# perform the simulation for the different phis using the ARSim function
xSimulations = matrix(ncol = T, nrow = length(phiVals))
set.seed(12345)
for (i in 1:length(phiVals)) {
  x = ARSim(phiVals[i])
  xSimulations[i,] = x
}

par(mfrow=c(3,3))
plot(xSimulations[1,], type="l", ylab = "x", main = "Traceplot, phi = -1")
plot(xSimulations[2,], type="l", ylab = "x", main = "Traceplot, phi = -0.9")
plot(xSimulations[5,], type="l", ylab = "x", main = "Traceplot, phi = -0.6")
plot(xSimulations[8,], type="l", ylab = "x", main = "Traceplot, phi = -0.3")
plot(xSimulations[11,], type="l", ylab = "x", main = "Traceplot, phi = 0")
plot(xSimulations[14,], type="l", ylab = "x", main = "Traceplot, phi = 0.3")
plot(xSimulations[17,], type="l", ylab = "x", main = "Traceplot, phi = 0.6")
plot(xSimulations[20,], type="l", ylab = "x", main = "Traceplot, phi = 0.9")
plot(xSimulations[21,], type="l", ylab = "x", main = "Traceplot, phi = 1")

# alternative for more efficient plotting
par(mfrow=c(3,3))
plot(xSimulations[1,], type="l", ylab = "x", main = "Traceplot, phi = -1")
for (i in 1:(length(phiVals)/3)) {
  print(3*i)
  plot(xSimulations[3*i,], type="l", ylab = "x", main = paste("Traceplot, phi = ", phiVals[3*i]))
  
}

par(mfrow = c(1,1))

############################### 3b #####################################
# Use your function from a) to simulate two AR(1)-processes, x1:T with phi = 0.2 and y1:T with phi = 0.95.
# Now, treat your simulated vectors as synthetic data, and treat the values of mu, phi and sigma^2 as unknown parameters.

# Implement Stan-code that samples from the posterior of the three parameters, using suitable non-informative priors of your choice.
# [Hint: Look at the time-series models examples in the Stan user's guide/reference manual, 
# and note the different parameterization used here.]

# i. Report the posterior mean, 95% credible intervals,
# and the number of effective posterior samples for the three inferred parameters for each of the simulated AR(1)-process.
# Are you able to estimate the true values? 

# ii. For each of the two data sets, evaluate the convergence of the samplers 
# and plot the joint posterior of mu and phi. Comments?


# Stan model which has uninformative priors, taken from Stan documentation for AR(1)-process
StanModel = '
data { 
  int<lower=0> N; 
  vector[N] y;    
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}
model {
  for (n in 2:N)
    y[n] ~ normal(mu + phi * (y[n-1] - mu) , sigma);
}
'

# x: phi = 0.2
x = ARSim(0.2)
data_x = list(N=T, y=x)
set.seed(12345)
#default warmup (1000), default iterations (2000), default chains = 4
fit_x = stan(model_code = StanModel, data=data_x)

# Print the fitted model
print(fit_x,digits_summary=3)

# Extract posterior samples
postDraws_x <- extract(fit_x)

# Do automatic traceplots of all chains
traceplot(fit_x)

# plot(postDraws_x$mu, type="l",ylab="mu",main="Traceplot") # Traceplots of the first chain
# pairs(fit_x) # Bivariate posterior plots


# y: phi = 0.95
y = ARSim(0.95)
data_y = list(N=T, y=y)
set.seed(12345)
fit_y =stan(model_code = StanModel, data=data_y)

# Print the fitted model
print(fit_y,digits_summary=3)

# Extract posterior samples
postDraws_y <- extract(fit_y)

# Do automatic traceplots of all chains
traceplot(fit_y)

# plot joint posterior of mu and phi
plot(postDraws_x$mu, postDraws_x$phi,ylab="phi", xlab="mu")
plot(postDraws_y$mu, postDraws_y$phi,ylab="phi", xlab="mu")

