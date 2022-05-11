library(rstan)


##3a 
mu = 13
sigma2 = 3
T = 300


#simulate values for x2, x3, ..., xT
ARSim = function(phi) {
  x = c(1:T)
  x[1] = mu
  for (i in 2:T) {
    # where the rnorm() term corresponds to epsilon.
    x[i] = mu + phi * (x[i-1] - mu) + rnorm(n=1, mean = 0, sd = sqrt(sigma2))
  }
  return(x)
}
  

#look at some different realizations (simulations) of x1:T for values of PHI between -1 and 1 
#(this is the interval where the AR(1)-process is stationary)
phiVals = seq(-1, 1, 0.1)

xSimulations = matrix(ncol = T, nrow = length(phiVals))
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


par(mfrow=c(3,3))
plot(xSimulations[1,], type="l", ylab = "x", main = "Traceplot, phi = -1")
for (i in 1:(length(phiVals)/3)) {
  print(3*i)
  plot(xSimulations[3*i,], type="l", ylab = "x", main = paste("Traceplot, phi = ", phiVals[3*i]))
  
}
#plot(..., main=paste("Jaccard vs. Tsim for depths",  min.depth, "to",max.depth,"m", sep=" "))


######3b
# Use your function from a) to simulate two AR(1)-processes, x1:T with phi = 0.2 and y1:T with phi = 0.95.
# Now, treat your simulated vectors as synthetic data, and treat the values of mu, phi and sigma^2 as unknown parameters.
# Implement Stan-code that samples from the posterior of the three parameters, using suitable non-informative priors of your choice.
# [Hint: Look at the time-series models examples in the Stan user's guide/reference manual, 
# and note the different parameterization used here.]



