###############################3a#####################################
degrees = c(285, 296, 314, 20, 299, 296, 40, 303, 326, 308)
radians = c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)

#?Bessel #mju - mean direcction, K concentration parameter
#Large k -> small variance around mju, and vice versa. 


mju = 2.51
lambda = 1
besselDegree = 0


#posterior proportioanl to likelihood * prior
#likelihood proportional to exp(k*cos(y-mju)) / I0(k), mju=2.51
#K~Exp(lambda=1) which gives the prior: 
#prior = lambda*exp(-lambda*y), lambda = 1


posterior = function(k, y, mju, degree, lambda) {
  #posterior proportional to this:
  #posterior = sum( exp(k*cos(y-mju))/besselI(k, degree) ) * lambda*exp(-lambda*k)
  #
  posterior = prod( exp(k*cos(y-mju))/besselI(k, degree) ) * dexp(k)
}

 
#gör om neg radians till positiva: 

k = seq(0.1, 8, 0.005)
g = posterior(0.1, radians, mju, besselDegree, lambda)

xD = c(1:length(k))
for(i in 1:length(k)) {
  xD[i] = posterior(k[i], radians, mju, besselDegree, lambda)
  
}


plot(k, xD, type="l")
o
xD = xD/sum(xD)

plot(k, xD, type="l")


###############################3b#####################################

approx_post_K = k[which.max(xD)]
approx_post_K
#ctrl shift c