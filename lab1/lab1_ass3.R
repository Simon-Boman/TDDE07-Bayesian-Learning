###############################1a#####################################

#north is 0
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
  posterior = exp(k*cos(y-mju))/besselI(k, degree) * exp(-lambda*y)
}

 
#gör om neg radians till positiva: 
g = posterior(1,2, mju, besselDegree, 1)
