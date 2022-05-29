par(mfrow=c(2,2))

### Plot the exact Beta posterior vs the approximate normal posterior
n <- 8
s <- 2
f <- n - s
p <- s/n
Grid <- seq(0,1,0.001)

plot(Grid,dbeta(Grid,1+s,1+f),col="blue")
lines(Grid,dnorm(Grid,mean=p,sd=sqrt(p*(1-p)/n)),col="red")

### Plot the exact Beta posterior vs the approximate normal posterior
n <- 2*8
s <- 2*2
f <- n - s
p <- s/n

plot(Grid,dbeta(Grid,1+s,1+f),col="blue")
lines(Grid,dnorm(Grid,mean=p,sd=sqrt(p*(1-p)/n)),col="red")

### Plot the exact Beta posterior vs the approximate normal posterior
n <- 4*8
s <- 4*2
f <- n - s
p <- s/n

plot(Grid,dbeta(Grid,1+s,1+f),col="blue")
lines(Grid,dnorm(Grid,mean=p,sd=sqrt(p*(1-p)/n)),col="red")

### Plot the exact Beta posterior vs the approximate normal posterior
n <- 100*8
s <- 100*2
f <- n - s
p <- s/n

plot(Grid,dbeta(Grid,1+s,1+f),col="blue")
lines(Grid,dnorm(Grid,mean=p,sd=sqrt(p*(1-p)/n)),col="red")