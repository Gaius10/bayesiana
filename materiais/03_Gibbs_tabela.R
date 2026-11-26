# Gibbs para tabela 2 x 2
set.seed(20261015)
# linhas: X=0/1; colunas: Y=0/1
P <- matrix(c(.30,.10,.20,.40), nrow=2)
stopifnot(abs(sum(P)-1)<1e-12)
N <- 12000; burnin <- 2000
X <- Y <- integer(N)
for(t in 2:N) {
  y0 <- Y[t-1] + 1L
  px1 <- P[2,y0]/sum(P[,y0])
  X[t] <- rbinom(1,1,px1)
  x0 <- X[t] + 1L
  py1 <- P[x0,2]/sum(P[x0,])
  Y[t] <- rbinom(1,1,py1)
}
observada <- table(factor(X[(burnin+1):N],levels=0:1),
                  factor(Y[(burnin+1):N],levels=0:1))
cat('Conjunta teorica:\n'); print(P)
cat('Conjunta simulada:\n'); print(prop.table(observada))
cat('P(X=1) simulada:', mean(X[(burnin+1):N]), '\n')
cat('P(Y=1) simulada:', mean(Y[(burnin+1):N]), '\n')
par(mfrow=c(2,1))
plot(X[1:300], type='s', main='Cadeia de X', ylab='X')
plot(Y[1:300], type='s', main='Cadeia de Y', ylab='Y')
par(mfrow=c(1,1))
