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



######3b