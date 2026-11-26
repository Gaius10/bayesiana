

library(R2jags)
library(coda)

# 1. Funcao de probabilidade da distribuicao discreta
fdp <- function(x, beta1) {
  pgamma(x + 1, shape = 2, scale = beta1) -
    pgamma(x, shape = 2, scale = beta1)
}

# 2. Geracao de uma amostra da distribuicao discreta

sim_categ <- function(N, beta1) {
  y <- rgamma(N, shape = 2, scale = beta1)
  x <- floor(y)
  return(x)
}

# 3. Modelo JAGS (ones trick)
# P(X=x) = exp(-x/beta) * [(1+x/beta)-(1+(x+1)/beta)*exp(-1/beta)]
# Como P(X=x) <= 1, usamos C = 1.
cat('model {
  for (i in 1:N) {
    p[i] <- exp(-x[i]/beta1) *
            ((1 + x[i]/beta1) -
            (1 + (x[i]+1)/beta1)*exp(-1/beta1))
    ones[i] ~ dbern(p[i])
  }
  beta1 ~ dgamma(0.1, 0.01)
}', file = 'dmex.jags')

# 4. Configuracao da simulacao
set.seed(4000)
S <- 10                 # Fazer com pelo menos 500 replicas
valores_beta <- c(2, 4, 5)  # Testar com mais valores
tamanhos <- c(25, 50, 100)   #Fazer com tamanhos amostrais maiores

# Resultado final: uma linha para cada combinacao de beta e N
df <- data.frame()

# 5. Repeticao dos cenarios
for (beta0 in valores_beta) {
  for (N in tamanhos) {

    # Guardar uma estimativa e um intervalo por replica
    beta_c <- rep(NA_real_, S)
    beta_inf <- rep(NA_real_, S)
    beta_sup <- rep(NA_real_, S)
    dp_post <- rep(NA_real_, S)
    rhat <- rep(NA_real_, S)
    ess <- rep(NA_real_, S)

    for (j in 1:S) {
      # 5.1. Gerar dados usando o valor verdadeiro beta0
      x <- sim_categ(N, beta0)
      ones <- rep(1, N)
      dados <- list(N = N, x = x, ones = ones)

      # 5.2. Usar tres cadeias com iniciais diferentes
      # Sementes do JAGS definidas explicitamente para reproducibilidade
      semente <- 100000 + 10000*match(beta0, valores_beta) +
                 1000*match(N, tamanhos) + 3*j
      iniciais <- list(
        list(beta1 = 1, .RNG.name = 'base::Mersenne-Twister', .RNG.seed = semente),
        list(beta1 = 3, .RNG.name = 'base::Mersenne-Twister', .RNG.seed = semente+1),
        list(beta1 = 6, .RNG.name = 'base::Mersenne-Twister', .RNG.seed = semente+2)
      )

      # 5.3. Ajustar o modelo
      ajuste <- jags(data = dados,
                     inits = iniciais,
                     parameters.to.save = c('beta1'),
                     n.chains = 3,
                     n.iter = 20000,
                     n.burnin = 5000,
                     n.thin = 10,
                     model.file = 'dmex.jags',
                     DIC = FALSE)

      # 5.4. Resumir a amostra a posteriori
      cadeias <- as.mcmc(ajuste)
      amostras <- as.numeric(ajuste$BUGSoutput$sims.list$beta1)

      beta_c[j] <- mean(amostras)
      beta_inf[j] <- unname(quantile(amostras, 0.025))
      beta_sup[j] <- unname(quantile(amostras, 0.975))
      dp_post[j] <- sd(amostras)

      # Diagnosticos: Rhat (Gelman-Rubin) e ESS
      rhat[j] <- gelman.diag(cadeias, autoburnin = FALSE)$psrf[1, 1]
      ess[j] <- as.numeric(coda::effectiveSize(cadeias)[1])

      cat('beta =', beta0, '| N =', N, '| replica =', j, 'de', S, '\n')
    }

    # 6. Medidas de desempenho no cenario
    vies <- mean(beta_c - beta0)
    eqm <- mean((beta_c - beta0)^2)
    rmse <- sqrt(eqm)
    pc <- mean(beta_inf <= beta0 & beta0 <= beta_sup)
    comprimento <- mean(beta_sup - beta_inf)

    linha <- data.frame(beta = beta0, N = N,
                        media = mean(beta_c), vies = vies,
                        eqm = eqm, rmse = rmse, pc = pc,
                        dp_posterior = mean(dp_post),
                        comprimento = comprimento,
                        rhat_medio = mean(rhat),
                        prop_rhat_101 = mean(rhat < 1.01),
                        ess_medio = mean(ess))
    df <- rbind(df, linha)

    # Salva a tabela sempre que termina um cenario
    write.csv(df, 'resumo_simulacao_bayes.csv', row.names = FALSE)
    print(round(linha, 4))
  }
}

# 7. Tabela final
print(round(df, 4))
write.csv(df, 'resumo_simulacao_bayes.csv', row.names = FALSE)
