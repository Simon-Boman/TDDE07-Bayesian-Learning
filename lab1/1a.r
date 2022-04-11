###### 1 a)

# Parameters - Prior distribution and data
s <- 13
n <- 50
p <- s/n                      #0.26
alpha <- 5
beta <- 5

# Parameters - Posterior distribution
alpha_post <- alpha+n*p        #5+13=18
beta_post <- beta + n*(1-p)    #5+37=42


# True values for mean and standard deviation
mean <- alpha_post / (alpha_post + beta_post)        # 0.3
sd <- sqrt((alpha_post*beta_post)/((alpha_post+beta_post)*(alpha_post+beta_post)*(alpha_post+beta_post+1)))    # 0.06
           
# Draw random numbers from posterior, want to draw from multiple sample sizes
# and plot a graph to see that the mean and sd converges to the true mean and sd
set.seed(13579)  
N = 100

v_mean <- c(0:4)
v_sd <- c(0:4)
v_sample <- c(0:4)
for (k in 1:5){
  v_sample[k] <- N
  
  y_rbeta <- rbeta(N, shape1 = alpha_post, shape2 = beta_post) 
  temp_mean = mean(y_rbeta)
  temp_sd = sd(y_rbeta)
  v_mean[k] <- abs(mean - temp_mean)
  v_sd[k] <- abs(sd - temp_sd)
 # plot(density(y_rbeta),                             
     #  main = "beta Distribution in R")
  N = N*5
}
v_mean

plot(v_sample, v_mean, col= "blue", ylab="Y", xlab="Sample size", type="l", ylim= c(0,0.002))
points(v_sample, v_sd, col= "green", type="l")






###### 1 b)

N = 10000
y_rbeta <- rbeta(N, shape1 = alpha_post, shape2 = beta_post)
y_rbeta <- ifelse(y_rbeta < 0.3, 1, 0)
true_cases = sum(y_rbeta)
prob_true_cases = true_cases/N
prob_true_cases

x_pbeta <- seq(0, 1, by = 0.0001)   

y_pbeta <- pbeta(x_pbeta, shape1 = alpha_post, shape2 = beta_post)

plot(y_pbeta)
print(y_pbeta[3000])



###### 1 c)

N = 10000
y_rbeta <- rbeta(N, shape1 = alpha_post, shape2 = beta_post)
log_odds = log(y_rbeta/(1-y_rbeta))
log_odds
hist(log_odds)
plot(density(log_odds))



