# R Code comparing Delta Method, Parametric Bootstrap, and Nonparametric Bootstrap
#
# This script performs the following steps:
# 1. Compute SE and 95% CI for \theta using the Delta Method.
# 2. Compute SE and 95% CI using the Parametric Bootstrap.
# 3. Compute SE and 95% CI using the Nonparametric Bootstrap.
# 4. Compare the results from all methods.
# 5. Plot histograms of bootstrap replications.


# ________________________________________________set the coding environment
rm(list = objects()); gc()   # empty the global environment

options(
  java.parameters = c(
    "-XX:+UseG1GC",
    "-Xms3072m"
  )
) # Allocates more heap space for Java and uses a different garbage collector.


if(!require( "pacman")){
  install.packages("pacman")
}    # prepare package which helps package loading

pacman::p_load(
  boot # for nonparametric bootstrap
) # load necessary packages

# __________________________________________________set the seed for reproducibility
set.seed(1122)


# __________________________________________________generate dataset
mu <- 5
n <- 100
sigma <- 1
data <- rnorm(n,mu,sigma)

#_______________________________________________get SE and 95% CI using delta method
se.delta <- exp(mean(data))/sqrt(n)

z<- qnorm(0.975)

ci.delta <- c(exp(mean(data))-z*se.delta,exp(mean(data))+z*se.delta)

theta_hat <- exp(mean(data))
# _______________________________________________________nonparametric bootstrap
theta.statistic<-function(Pop,indices){
  Pop <- Pop[indices]
  theta<-exp(mean(Pop))
  return(theta)
}

boot.results<-boot(data, theta.statistic, R=500)
se.boot<-sqrt(var(boot.results$t))

ci.boot.nonparametric <- boot.ci(boot.results, conf=0.95, type="perc")$percent[4:5]

# __________________________________________________________parametric bootstrap
mu.estimate <- mean(data)

B <- 500
theta.estimate <-c()
for(b in 1:B){
  b.data <-rnorm(100, mu.estimate, 1)
  theta.estimate [b] <-exp(mean(b.data))
}

q.lo =quantile(theta.estimate, prob=0.05/2)
q.hi =quantile(theta.estimate, prob=1-0.05/2)

se.nonparam <- sd(theta.estimate)

ci.boot.parametric <- c(2*exp(mu.estimate)-q.hi,2*exp(mu.estimate)-q.lo)

#______________________________________________________ compare results
results <- data.frame(
  Method = c("Delta Method", "Nonparametric Bootstrap", "Parametric Bootstrap"),
  Theta_Hat = c(theta_hat, mean(boot.results$t), mean(theta.estimate)),
  SE = c(se.delta, se.boot, se.nonparam),
  CI_Lower = c(ci.delta[1], ci.boot.nonparametric[1], ci.boot.parametric[1]),
  CI_Upper = c(ci.delta[2], ci.boot.nonparametric[2], ci.boot.parametric[2])
)
print(results)

# _________________________________________________histogram of estimates of theta
par(mfrow=c(1,2), family="Times")
hist(boot.results$t, col="orange",
     main="Samples from nonparametric bootstrap",
     xlab=expression(hat(theta)))
abline(v=exp(mu), lwd=3, lty=2)

hist(theta.estimate, col="lightblue",
     main="Samples from parametric bootstrap",
     xlab=expression(hat(theta)))
abline(v=exp(mu), lwd=3, lty=2)

ci.true <- c(exp(mu-z*1/10),exp(mu+z*1/10))
