# Metropolis-Hastings para uma distribuicao em {1, 2, 3}
set.seed(20261015)
pi_alvo <- c(0.2, 0.5, 0.3)
N <- 10000
x <- integer(N)
x[1] <- 1L
aceitos <- 0L
for (t in 2:N) {
  atual <- x[t-1]
  candidato <- sample(setdiff(1:3, atual), 1)
  alfa <- min(1, pi_alvo[candidato]/pi_alvo[atual])
  if (runif(1) < alfa) {
    x[t] <- candidato
    aceitos <- aceitos + 1L
  } else x[t] <- atual
}
cat('Distribuicao esperada:', pi_alvo, '\n')
cat('Frequencias simuladas:', prop.table(table(factor(x[1001:N], levels=1:3))), '\n')
cat('Taxa de aceitacao:', aceitos/(N-1), '\n')
par(mfrow=c(1,2))
plot(x[1:250], type='s', xlab='Iteração', ylab='Estado', main='Início da cadeia')
barplot(prop.table(table(factor(x[1001:N], levels=1:3))), ylim=c(0,0.6),
        names.arg=1:3, main='Frequencias após o aquecimento')
par(mfrow=c(1,1))
