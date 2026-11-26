# ESTUDO DE SIMULACAO CLASSICO 
# X = floor(Y), com Y ~ Gama(forma = 2, escala = beta)
# Mesmos cenarios do estudo bayesiano: beta = 2,3,4,5; N = 25,50,100,200,400.

# 1. Funcao de probabilidade (para consulta)
fdp <- function(x, beta1) {
  pgamma(x + 1, shape = 2, scale = beta1) -
    pgamma(x, shape = 2, scale = beta1)
}

# 2. Geracao dos dados, sem limitar os valores a 0,...,99
sim_categ <- function(N, beta1) {
  y <- rgamma(N, shape = 2, scale = beta1)
  x <- floor(y)
  return(x)
}

# 3. Log-verossimilhanca
# A otimizacao usa eta = log(beta1), portanto beta1 = exp(eta) > 0.
# A formula usa expm1() para evitar a subtracao de numeros muito proximos.
ll <- function(eta, x) {
  beta1 <- exp(eta)
  u <- 1 / beta1
  # P(X=x) = exp(-x*u) * [P(Gama(2,1) <= u) + x*u*(1-exp(-u))]
  parte <- pgamma(u, shape = 2, rate = 1) + x*u*(-expm1(-u))
  if (!is.finite(beta1) || any(!is.finite(parte)) || any(parte <= 0)) {
    return(-1e100)
  }
  return(sum(-x*u + log(parte)))
}

# 4. Configuracao do estudo
set.seed(4000)
S <- 10   # Para testar: 10; estudo definitivo: 1000
valores_beta <- c(2, 3, 4, 5)
tamanhos <- c(25, 50, 100, 200, 400)
z <- qnorm(0.975)

df <- data.frame()

# 5. Repetir os cenarios
for (beta0 in valores_beta) {
  for (N in tamanhos) {

    beta_c <- rep(NA_real_, S)
    beta_inf <- rep(NA_real_, S)
    beta_sup <- rep(NA_real_, S)
    erro_padrao <- rep(NA_real_, S)
    conv <- rep(FALSE, S)

    for (j in 1:S) {
      # 5.1. Gerar a amostra
      x <- sim_categ(N, beta0)

      # 5.2. Encontrar a estimativa de maxima verossimilhanca
      ajuste <- optim(par = log(max(mean(x)/2, 0.1)),
                      fn = ll, x = x,
                      method = 'BFGS',
                      control = list(fnscale = -1),
                      hessian = TRUE)

      beta_c[j] <- exp(ajuste$par)

      # 5.3. Erro-padrao e intervalo de confianca
      # H_eta e a segunda derivada da log-verossimilhanca em eta.
      # Var(eta_hat) ~ 1/(-H_eta); pelo metodo delta,
      # EP(beta_hat) ~ beta_hat / sqrt(-H_eta).
      H_eta <- as.numeric(ajuste$hessian[1, 1])
      valido <- (ajuste$convergence == 0 && is.finite(H_eta) && H_eta < 0)
      conv[j] <- valido

      if (valido) {
        erro_padrao[j] <- beta_c[j] / sqrt(-H_eta)
        # Intervalo de Wald na escala log(beta): sempre positivo
        ep_eta <- 1 / sqrt(-H_eta)
        beta_inf[j] <- exp(ajuste$par - z*ep_eta)
        beta_sup[j] <- exp(ajuste$par + z*ep_eta)
      }

      cat('beta =', beta0, '| N =', N,
          '| replica =', j, 'de', S, '\n')
    }

    # 6. Medidas de desempenho
    # Viés e EQM são calculados para todas as estimativas finitas.
    boas <- is.finite(beta_c)
    int_validos <- is.finite(beta_inf) & is.finite(beta_sup)
    vies <- if (any(boas)) mean(beta_c[boas] - beta0) else NA_real_
    eqm <- if (any(boas)) mean((beta_c[boas] - beta0)^2) else NA_real_
    rmse <- sqrt(eqm)
    # Cobertura entre intervalos validos; informar também a proporção valida.
    pc <- if (any(int_validos)) {
      mean(beta_inf[int_validos] <= beta0 & beta0 <= beta_sup[int_validos])
    } else NA_real_
    comprimento <- if (any(int_validos)) mean(beta_sup[int_validos] - beta_inf[int_validos]) else NA_real_
    ep_medio <- if (any(is.finite(erro_padrao))) mean(erro_padrao[is.finite(erro_padrao)]) else NA_real_

    linha <- data.frame(beta = beta0, N = N,
                        media = if (any(boas)) mean(beta_c[boas]) else NA_real_,
                        vies = vies, eqm = eqm, rmse = rmse,
                        pc = pc, ep_medio = ep_medio,
                        comprimento = comprimento,
                        prop_convergencia = mean(conv),
                        prop_ic_valido = mean(int_validos))
    df <- rbind(df, linha)

    # Guardar resultados apos cada cenario
    write.csv(df, 'resumo_simulacao_classica.csv', row.names = FALSE)
    print(round(linha, 4))
  }
}

# 7. Tabela final
print(round(df, 4))
write.csv(df, 'resumo_simulacao_classica.csv', row.names = FALSE)
