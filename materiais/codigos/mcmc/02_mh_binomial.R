# Metropolis-Hastings independente para Binomial-Beta
set.seed(20261015)
y <- 7; n <- 10; a <- 2; b <- 2
N <- 12000; burnin <- 2000
p <- numeric(N); p[1] <- 0.2
aceitos <- 0L
log_nucleo <- function(z) {
  if (z <= 0 || z >= 1) return(-Inf)
  (y+a-1)*log(z) + (n-y+b-1)*log1p(-z)
}
for (t in 2:N) {
  cand <- runif(1)
  log_alpha <- min(0, log_nucleo(cand)-log_nucleo(p[t-1]))
  if (log(runif(1)) < log_alpha) {
    p[t] <- cand
    aceitos <- aceitos+1L
  } else p[t] <- p[t-1]
}
amostra <- p[(burnin+1):N]
cat('Media exata:', (a+y)/(a+b+n), '\n')
cat('Media MCMC:', mean(amostra), '\n')
cat('IC 95%:', quantile(amostra, c(.025,.975)), '\n')
cat('Aceitacao:', aceitos/(N-1), '\n')
par(mfrow=c(1,2))
plot(p, type='l', main='Traco de p', xlab='Iteracao', ylab='p')
hist(amostra, probability=TRUE, breaks=35, main='Posteriori de p', xlab='p')
curve(dbeta(x, a+y, b+n-y), add=TRUE, col='red', lwd=2)
par(mfrow=c(1,1))
