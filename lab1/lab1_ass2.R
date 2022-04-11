###############################1a#####################################

data = c(33, 24, 48, 32, 55, 74, 23, 76, 17)
n <- length(data)

mju = 3.5
tao = sum((log(data)-mju)^2)/n

NDraws = 10000
PostDraws <- c(1:NDraws)
set.seed(12345)
PostDraws <- ((n)*tao)/rchisq(NDraws,n) 
#PostDraws[,1] <- rlnorm(NDraws, mean=tao, sd=sqrt(PostDraws[,2]/n))

plot(density(PostDraws))

###############################1b#####################################

#phi(z)  cdf för vanliga normal distribution with mean 0 and variance 1
#pnorm - gives distribution function
G = 2*pnorm(sqrt(PostDraws/2)) - 1
plot(density(G))

###############################1c#####################################

mean = mean(G)
sd = sd(G)

lowerBound = mean - 1.96*sd
upperBound = mean + 1.96*sd 
lowerBound
upperBound

abline(v=lowerBound, col="red")
abline(v=upperBound, col="red")
#varför olika massa på de olika sidorna?

qnorm(0.025, mean, sd)
qnorm(0.975, mean, sd)

###############################1d#####################################

density = density(G)
density

sortedG = sort(G)
densitySorted = density(sortedG)
densitySorted
