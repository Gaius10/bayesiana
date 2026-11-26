
#
set.seed(123)
u <- runif(10000)
head(u)
hist(u, breaks = 30, probability = TRUE,
     main = "Uniforme(0,1)", xlab = "u")
abline(h = 1, col = "red", lwd = 2)
mean(u); var(u)

# 
set.seed(123)
u <- runif(10000)
x <- ifelse(u <= 0.25, 0,
            ifelse(u <= 0.75, 1, 2))
prop.table(table(x))
barplot(prop.table(table(x)),
        main = "Probabilidades simuladas")

# 
set.seed(123)
n <- 10000; p <- 0.6
u <- runif(n)
x <- as.integer(u > 1 - p)
mean(x)                  # deve se aproximar de p
mean(rbinom(n, 1, p))    # funcao pronta do R

# 
set.seed(123)
N <- 1000; m <- 10; p <- 0.4
b <- matrix(rbinom(N*m, 1, p), nrow = N)
x <- rowSums(b)
mean(x); var(x)
cbind(k = 0:m,
  sim = tabulate(x + 1, nbins = m+1)/N,
  teor = dbinom(0:m, m, p))

# 
set.seed(123)
p <- 0.3; u <- runif(10000)
x <- floor(log1p(-u)/log1p(-p))
mean(x)                    # (1-p)/p
table(x)[1:8] / length(x)
plot(0:10, dgeom(0:10, p), type = "h",
     xlab = "falhas", ylab = "Probabilidade")

# 
rpois_inv <- function(n, lambda) {
  u <- runif(n)
  vapply(u, function(t) {
    k <- 0; pr <- exp(-lambda); acum <- pr
    while (t > acum) {
      k <- k + 1
      pr <- pr * lambda / k
      acum <- acum + pr
    }
    k
  }, numeric(1))
}
set.seed(123)
x <- rpois_inv(2000, 3)
mean(x); mean(rpois(2000, 3))

# 
set.seed(123)
lambda <- 2; n <- 10000
u <- runif(n)
x <- -log1p(-u)/lambda
mean(x)                   # 1/lambda = 0.5
hist(x, breaks = 40, probability = TRUE,
     main = "Exponencial por inversao")
curve(dexp(x, rate = lambda), add = TRUE,
      col = "red", lwd = 2)

#
set.seed(123)
n <- 10000
u1 <- runif(n); u2 <- runif(n)
z1 <- sqrt(-2*log(u1))*cos(2*pi*u2)
z2 <- sqrt(-2*log(u1))*sin(2*pi*u2)
c(mean(z1), sd(z1), cor(z1,z2))
hist(z1, breaks = 40, probability = TRUE)
curve(dnorm(x), add = TRUE, col = "red", lwd=2)

#
set.seed(123)
Sigma <- matrix(c(1, .5, .5, 1), 2, 2)
L <- t(chol(Sigma))
z <- matrix(rnorm(2*5000), nrow=2)
y <- t(L %*% z)
cor(y[,1], y[,2])
plot(y, pch=16, cex=.3, xlab="X", ylab="Y")
# Alternativa: MASS::mvrnorm(5000,c(0,0),Sigma)

# 
set.seed(123)
n <- 5000; saida <- numeric(0); tent <- 0
while (length(saida) < n) {
  y <- runif(1); u <- runif(1); tent <- tent+1
  if (u < 4*y*(1-y)) saida <- c(saida,y)
}
length(saida)/tent
hist(saida, probability=TRUE, breaks=30)
curve(dbeta(x,2,2),add=TRUE,col="red",lwd=2)

#
set.seed(123)
int_exp <- function(N, a=1, b=3) {
  u <- runif(N, a, b)
  (b-a)*mean(exp(-u))
}
exato <- exp(-1)-exp(-3)
c(exato=exato, MC100=int_exp(100),
  MC10000=int_exp(10000), MC1e5=int_exp(1e5))
integrate(function(x) exp(-x),1,3)

#
set.seed(123)
N <- 10000
u <- runif(N, 1, 3)
y <- 2*exp(-u)
est <- mean(y)
ep <- sd(y)/sqrt(N)
c(estimativa=est, EP=ep,
  LI=est-1.96*ep, LS=est+1.96*ep,
  valor_exato=exp(-1)-exp(-3))

# 
set.seed(123)
N <- c(10, 20, 50, 100, 200, 500, 1000)
reps <- 50
res <- sapply(N, function(nn) {
  replicate(reps, int_exp(nn))
})
boxplot(as.data.frame(res), names=N,
        xlab="Numero de simulacoes",
        ylab="Estimativas da integral")
abline(h=exp(-1)-exp(-3), col="red",lty=2)

# 
set.seed(123)
Sigma <- matrix(c(1,.5,.5,1),2,2)
y <- MASS::mvrnorm(100000,mu=c(0,0),Sigma)
p_hat <- mean(y[,1] < 1 & y[,2] < 1)
ep <- sqrt(p_hat*(1-p_hat)/nrow(y))
c(probabilidade=p_hat, EP=ep)
# Com o pacote pbivnorm instalado:
# pbivnorm::pbivnorm(1, 1, rho=0.5)

# 
set.seed(123)
theta <- rbeta(100000, shape1=14, shape2=4)
c(media_MC=mean(theta),
  media_exata=14/18,
  prob_MC=mean(theta>0.6),
  prob_exata=1-pbeta(0.6,14,4))
quantile(theta,c(.025,.5,.975))
qbeta(c(.025,.5,.975),14,4)

# 
set.seed(123)
x <- rnorm(50, mean=2, sd=1)
n <- length(x); xb <- mean(x)
amostras <- numeric(1000); total <- 0
for (i in seq_along(amostras)) {
  repeat {
    total <- total+1
    y <- rcauchy(1)
    if (runif(1) < exp(-n*(y-xb)^2/2)) {
      amostras[i] <- y
      break
    }
  }
}
mean(amostras); length(amostras)/total
