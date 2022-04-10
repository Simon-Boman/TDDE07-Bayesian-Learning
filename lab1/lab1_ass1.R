###############################1a#####################################
options(scipen=999) #for removing scientific notation on the plots

n = 50
s = 13
#alpha and beta in prior
alpha0 = 5
beta0 = 5
#alpha and beta in posterior
alpha = alpha0 + s
beta = beta0 + (n - s)

#true mean and sd
true_mean = alpha / (alpha + beta)
true_sd = sqrt(alpha * beta / ((alpha + beta)^2 * (alpha+beta+1)))


draws = seq(1000, 100000, 1000) #different values for number of sample size
means = c(1:length(draws)) #vector to store the means
sds = c(1:length(draws)) #vector to store the sds

#for each of the different number of draws, compute and store mean and sd
for (i in 1:length(draws)) {
  set.seed(12345)
  data = rbeta(draws[i], alpha, beta)
  mean = mean(data)
  sd = sd(data)
  means[i] = abs(true_mean - mean)
  sds[i] = abs(true_sd - sd)
}

#As we see in the plot, the difference of the mean/sd between the true value and the value from the simulations
#go to 0 as we increase the number of draws, i.e. the posterior mean and sd converges to the true values of 0.3 and 0.05867387. 
max_y = max(means, sds) 
plot(draws, means, type="l", lwd = 2, col="blue", ylim=c(0,max_y), main = "Deviance from true mean and sd for different number of draws",
     xlab="Number of draws", ylab="Deviance between true mean/sd and simulated mean/sd" )
points(draws, sds, type="l", lwd = 2, col="green")
legend(x = max(draws)*0.8, y = max_y*0.95, legend = c("mean", "sd"), col = c("blue","green"), lwd = 3)



#################################1b###################################

#pbeta gives the distribution function, i.e. the accumulated probability mass up to some value
true_prob = pbeta(0.3, alpha, beta)
true_prob

nDraws = 10000
set.seed(12345)
data = rbeta(nDraws, alpha, beta)
#using our draws, calculate the probability by taking the #samples < 0.3 divided by total #samples
sim_prob = length(data[data<0.3]) / nDraws 
sim_prob

#just for visualization 
plot(density(data))
abline(v=0.3, col="red")

#The probability from simulation with 10000 draws, 0.5156, is very close to the true probability of 0.5150226. 



###############################1c#####################################

nDraws = 10000
set.seed(12345)
data = rbeta(nDraws, alpha, beta)
#hist(data, breaks=25)
#plot(density(data))
log_odds = log(data / (1-data)) #logit maps the probability values from (0,1) to (-inf, inf)
hist(log_odds, breaks=25)
plot(density(log_odds), main = "Density estimation of log_odds")



