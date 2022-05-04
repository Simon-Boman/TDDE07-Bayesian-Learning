

data_original = readRDS("Precipitation.rds")

n = length(data)
data = log(data_original)

mu0 = 20
tao0 = 20

v0 = 300
sigma0 = 2.2


nDraws <- 10000
gibbsDraws <- matrix(0,nDraws,2)
colnames(gibbsDraws) = c("mu", "sigma")

vn = v0 + n 

sigma_squared <- 1 # Initial value for sigma
for (i in 1:nDraws){  
  
  # Update mu given sigma
  tao_n_squared = 1/(n/sigma_squared + 1/tao0^2)
  w = (n/sigma_squared) / (n/sigma_squared + 1/tao0^2) 
  mu_n = w*mean(data) + (1-w)*mu0
  
  mu = rnorm(n = 1, mean = mu_n, sd = tao_n_squared)
  gibbsDraws[i,1] <- mu
  
  # Update sigma given mu
  param2_squared = (v0*sigma0^2 + sum(data-mu)^2) / (n + v0)
  sigma_squared = (vn*param2_squared/rchisq(1, vn))
  gibbsDraws[i,2] <- sigma_squared
  
#vi har X^2 (vn, param2)
      # vn*param2 + 
  #     X^2 (v0, sigma0^2 )
  #    ((v0)*sigma0)/rchisq(nDraws,v0) 
}

gibbsDraws

a_Gibbs_mu <- acf(gibbsDraws[,1])
IF_Gibbs_mu <- 1+2*sum(a_Gibbs_mu$acf[-1])
IF_Gibbs_mu

a_Gibbs_sigma <- acf(gibbsDraws[,2])
IF_Gibbs_sigma <- 1+2*sum(a_Gibbs_sigma$acf[-1])
IF_Gibbs_sigma



plot(1:nDraws, gibbsDraws[,1], type = "l",col="red") # traceplot of Gibbs draws
hist(gibbsDraws[,1],col="red") # histogram of Gibbs draws
cusumData =  cumsum(gibbsDraws[,1])/seq(1,nDraws) # Cumulative mean value of mu, Gibbs draws
#plot(1:nDraws, cusumData, type = "l", col="red")
#barplot(height = a_Gibbs_mu$acf[-1],col="red") # acf for Gibbs draws


plot(1:nDraws, gibbsDraws[,2], type = "l",col="red") # traceplot of Gibbs draws
hist(gibbsDraws[,2],col="red") # histogram of Gibbs draws
cusumData =  cumsum(gibbsDraws[,2])/seq(1,nDraws) # Cumulative mean value of mu, Gibbs draws
#plot(1:nDraws, cusumData, type = "l", col="red")
#barplot(height = a_Gibbs_sigma$acf[-1],col="red") # acf for Gibbs draws



##################1b
hist(data_original, main = "original data")
plot(density(data_original))
posteriorDraws = rnorm(n = nDraws, mean = gibbsDraws[,1], sd = gibbsDraws[,2])
posteriorDraws = exp(posteriorDraws)
hist(posteriorDraws, main = "posterior draws")
plot(density(posteriorDraws))

plot(density(posteriorDraws), col = "red")
lines(density(data_original), col = "green")

