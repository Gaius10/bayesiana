# MH passeio aleatorio para N(0,1), proposta N(x_atual, s^2)
set.seed(20261015)
metropolis <- function(N,s,inicial=0) {
  x <- numeric(N); x[1] <- inicial; aceitos <- 0L
  for (t in 2:N) {
    candidato <- rnorm(1,mean=x[t-1],sd=s)
    logalpha <- min(0,-0.5*(candidato^2-x[t-1]^2))
    if (log(runif(1)) < logalpha) {
      x[t] <- candidato; aceitos <- aceitos + 1L
    } else x[t] <- x[t-1]
  }
  list(x=x, taxa=aceitos/(N-1))
}
valores <- c(0.2,1,5)
par(mfrow=c(3,2))
for (s in valores) {
  ans <- metropolis(10000,s)
  plot(ans$x[1:300],type='l',main=paste('s =',s,'aceitacao =',round(ans$taxa,2)),ylab='X')
  acf(ans$x[2001:10000],main=paste('ACF s =',s))
  cat('s =',s,'taxa =',ans$taxa,'media =',mean(ans$x[2001:10000]),'\n')
}
par(mfrow=c(1,1))
