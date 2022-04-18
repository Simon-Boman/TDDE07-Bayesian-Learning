###############################3a#####################################
degrees = c(285, 296, 314, 20, 299, 296, 40, 303, 326, 308)
radians = c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)

mju = 2.51
lambda = 1
besselDegree = 0

#posterior proportional to likelihood * prior
#likelihood proportional to exp(k*cos(y-mju)) / I0(k)
#K~Exp(lambda=1) which gives the prior: lambda*exp(-lambda*k)
posterior = function(k, y, mju, degree, lambda) {
  #posterior proportional to:
  posterior = prod( exp(k*cos(y-mju))/besselI(k, degree) ) * lambda*exp(-lambda*k)
  #posterior = prod( exp(k*cos(y-mju))/besselI(k, degree) ) * dexp(k, rate=lambda)
}

#grid of k-values to plot over
gridstep = 0.005
k = seq(0, 7, gridstep)

#calculate and store the value for each of the k-values
post = c(1:length(k))
for(i in 1:length(k)) {
  post[i] = posterior(k[i], radians, mju, besselDegree, lambda)
}

#normalizing the posterior distirbution of k so that it integrates to 1
post = post/(sum(post)*gridstep)
plot(k, post, type="l", main="Posterior distribution of k", col="blue")




###############################3b#####################################
#the approximate posterior mode of k is the largest value in the posterior
approx_post_K = k[which.max(post)]
approx_post_K
plot(k, post, type="l", main="Posterior distribution of k", col="blue")
abline(v=approx_post_K, col="red")




