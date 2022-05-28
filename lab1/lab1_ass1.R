###############################1a#####################################
options(scipen=999) # for removing scientific notation on the plots

# number of observations and number of successes
n = 50
s = 13
# alpha and beta in prior Beta distribution
alpha0 = 5
beta0 = 5
# alpha and beta in posterior Beta distribution
alpha = alpha0 + s
beta = beta0 + (n - s)

#true mean and standard deviation of the posterior Beta distribution (formulas from wikipedia for Beta distribution)
true_mean = alpha / (alpha + beta)
true_sd = sqrt(alpha * beta / ((alpha + beta)^2 * (alpha+beta+1)))


draws = seq(5000, 500000, 5000) #different number of draws, i.e. sample sizes to use 
means = c(1:length(draws)) #vector to store the means of the different sample sizes
sds = c(1:length(draws)) #vector to store the sd's of the different sample sizes

#for each of the different number of draws, compute and store the estimated mean and sd.
#storing the absolute value for easier comparison.
for (i in 1:length(draws)) {
  set.seed(12345)
  data = rbeta(draws[i], alpha, beta)
  mean = mean(data)
  sd = sd(data)
  means[i] = abs(true_mean - mean)
  sds[i] = abs(true_sd - sd)
}

#As we see in the plot, the difference between the true mean/sd and the estimated mean/sd from the simulations
#go to 0 as we increase the number of draws, i.e. the posterior mean and sd converges to the true values of 0.3 and 0.05867387. 
max_y = max(means, sds) 
plot(draws, means, type="l", lwd = 2, col="blue", ylim=c(0,max_y), main = "Deviance from true mean and sd for different number
     of draws", xlab="Number of draws", ylab="Deviance between true mean/sd and simulated mean/sd" )
points(draws, sds, type="l", lwd = 2, col="green")
legend(x = max(draws)*0.8, y = max_y*0.95, legend = c("mean", "sd"), col = c("blue","green"), lwd = 3)




#################################1b###################################
#pbeta gives the distribution function for the Beta distribution,
#i.e. the accumulated probability mass up to some value (area under curve less or equal to the value)
true_prob = pbeta(0.3, alpha, beta)
true_prob

nDraws = 10000
set.seed(12345)
data = rbeta(nDraws, alpha, beta)
#using our 10000 draws from the posterior,
#we can calculate the probability by taking number of samples < 0.3 divided by total number of samples
sim_prob = length(data[data<0.3]) / nDraws 
sim_prob

#The probability from simulation with 10000 draws, 0.5156, is very close to the exact probability of 0.5150. 

#Plot of density estimates of our data, a value of around 0.515 looks reasonable. 
plot(density(data))
abline(v=0.3, col="red")




###############################1c#####################################
nDraws = 10000
set.seed(12345)
data = rbeta(nDraws, alpha, beta)

#posterior distribution
hist(data, breaks=25)
plot(density(data), main = "Density estimation of posterior distribution")

#simulating draws from the posterior distribution of the log-odds
#(logit maps the probability values from (0,1) to (-inf, inf))
log_odds = log(data / (1-data)) 

#posterior distribution of the log-odds
hist(log_odds, breaks=25)
plot(density(log_odds), main = "Density estimation of posterior distribution of log_odds")




