# Gibbs Normal: Y_i ~ N(mu, 1/tau), mu ~ N(m0,s0^2), tau ~ Gama(a,b)
# Gama parametrizada com taxa (rate)
set.seed(20261015)
n <- 100; y <- rnorm(n, mean=5, sd=3)
m0 <- 0; s02 <- 100; a <- 0.5; b <- 0.5
N <- 14000; burnin <- 2000

gibbs_normal <- function(y, N, mu_inicial, tau_inicial) {
  stopifnot(tau_inicial > 0)
  n <- length(y)
  mu <- tau <- numeric(N)
  mu[1] <- mu_inicial; tau[1] <- tau_inicial
  for (t in 2:N) {
    V <- 1/(n*tau[t-1]+1/s02)
    M <- V*(tau[t-1]*sum(y)+m0/s02)
    mu[t] <- rnorm(1, mean=M, sd=sqrt(V))
    tau[t] <- rgamma(1, shape=a+n/2,
                    rate=b+sum((y-mu[t])^2)/2)
  }
  cbind(mu=mu, tau=tau, sigma2=1/tau)
}
set.seed(123)
cadeia1 <- gibbs_normal(y,N, mu_inicial=0, tau_inicial=0.1)
set.seed(456)
cadeia2 <- gibbs_normal(y,N, mu_inicial=10, tau_inicial=2)
pos1 <- cadeia1[(burnin+1):N, ,drop=FALSE]
pos2 <- cadeia2[(burnin+1):N, ,drop=FALSE]
pos <- rbind(pos1,pos2)
print(colMeans(pos))
print(apply(pos,2,quantile,probs=c(.025,.5,.975)))
cat('P(mu>0|y) estimada:', mean(pos[,'mu']>0),'\n')
par(mfrow=c(2,2))
matplot(cbind(cadeia1[,1],cadeia2[,1]),type='l', lty=1,
        main='Tracos de mu',xlab='Iteracao',ylab='mu')
matplot(cbind(cadeia1[,2],cadeia2[,2]),type='l',lty=1,
        main='Tracos de tau',xlab='Iteracao',ylab='tau')
acf(pos1[,'mu'],main='ACF mu')
acf(pos1[,'tau'],main='ACF tau')
par(mfrow=c(1,1))
if (requireNamespace('coda',quietly=TRUE)) {
  mcmc_list <- coda::mcmc.list(coda::mcmc(pos1[,c('mu','tau')]),
                              coda::mcmc(pos2[,c('mu','tau')]))
  print(coda::gelman.diag(mcmc_list))
  print(coda::effectiveSize(mcmc_list))
} else message('Pacote coda ausente: diagnosticos adicionais ignorados.')
