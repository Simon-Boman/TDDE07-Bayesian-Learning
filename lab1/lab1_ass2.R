###############################2a#####################################
data = c(33, 24, 48, 32, 55, 74, 23, 76, 17)
n <- length(data)

mju = 3.5
tao2 = sum((log(data)-mju)^2)/n

NDraws = 10000
PostDraws = c(1:NDraws)
set.seed(123465)
#drawing 10000 sigma^2 from the inv-chi^2 distribution
PostDraws = ((n)*tao2)/rchisq(NDraws,n)  

plot(density(PostDraws), xlim=c(0,2), main = expression(paste('Posterior distribution of ',sigma^2,)))




###############################2b#####################################
#phi(z) is the cumulative distribution function (CDF) for the standard normal distribution with mean 0 and variance 1.
#pnorm gives the distribution function for the Normal distribution with default mean=1 and sd=1 (the CDF),
#i.e. the accumulated probability mass up to some value (area under curve less or equal to the value)
G = 2*pnorm(sqrt(PostDraws/2)) - 1
plot(density(G), main = "Posterior distribution of the Gini coefficient G for the current data set ")




###############################2c#####################################
#estimation of the the 95% equal tail credible interval using mean and sd
mean = mean(G) 
sd = sd(G)
lowerBound = mean - 1.96*sd
upperBound = mean + 1.96*sd 
lowerBound
upperBound

#more exact values using qnorm
lowerBound = qnorm(0.025, mean, sd)
upperBound = qnorm(0.975, mean, sd)
lowerBound
upperBound

plot(density(G), main = "Posterior distribution of the Gini coefficient G for the current data set, \n with 95% equal tail credible interval. ")
abline(v=lowerBound, col="red")
abline(v=upperBound, col="red")




###############################2d#####################################
#extract x- and y-values from the density estimate and sort these on y-value 
#density = density(G, n=10000)
density = density(G)
x = density$x
y = density$y
xy = cbind(x,y)
xy_sorted = xy[order(xy[,2],decreasing=TRUE),]

#calculate what value 95% of the accumulated y-values corresponds to, i.e. the accumulated value to stop at 
#when deciding which (x,y)-pairs to look at
stop_value = sum(xy_sorted[,2])*0.95

#accumulated sum of the sorted y-values
#find out at which index we reach 95% mass (i.e. the stop_value),
#then index 1 to stop_index corresponds to 95% mass
acumulatedMass <- function(stop_val) {
  for (i in 1:length(xy_sorted)) {
    if (sum(xy_sorted[1:i,2]) > stop_val) {
      return (i-1)
    }
  }
}
stop_index = acumulatedMass(stop_value)

#alternative to using above function:
#cumSum = cumsum(xy_sorted[,2])
#stop_index = min(which(cumSum > stop_value)) - 1

#only keep the 95% of the mass corresponding to highest densities
xy_sorted_95 = xy_sorted[1:stop_index,]

#find the indexes corresponding to lower and upper bound
lower_index = which.min(xy_sorted_95[,1])
upper_index = which.max(xy_sorted_95[,1])
lower_index
upper_index

#extract actual lower and upper bounds
HPDIlower = xy_sorted_95[lower_index,1]
HPDIupper = xy_sorted_95[upper_index,1]
HPDIlower
HPDIupper

plot(density(G), main = "Posterior distribution of the Gini coefficient G for the current data set, \n with 95% Equal Tail Credible Interval, \n and 95% Highest Posterior Density Interval HPDI ")
legend(x = 0.6, y = 5, legend = c("Equal Tail Credible Interval", "HPDI"), col = c("red","blue"), lwd = 3)
abline(v=lowerBound, col="red")
abline(v=upperBound, col="red")
abline(v=HPDIlower, col="blue")
abline(v=HPDIupper, col="blue")



