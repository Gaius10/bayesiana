

# Dados da China
dados_china <- c(3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9, 10, 11, 11, 13,
                 3, 14, 15, 16, 17, 22, 22, 24, 26, 26, 27, 28, 29, 30,
                 31, 31, 35, 38, 38, 42, 43, 44, 45, 46, 47, 52, 57, 64,
                 65, 71, 73, 73, 86, 89, 97, 97, 97, 98, 105, 108, 109,
                 114, 118, 121, 136, 142, 143, 146, 150)

# Dados da Europa
dados_europa <- c(6, 18, 28, 29, 44, 47, 55, 116, 118, 129, 150,
                  184, 219, 236, 237, 336, 421, 434, 612, 648, 706,
                  838, 1129, 1393, 1540, 1941, 2175, 2278, 2667,
                  2803, 2824)



# Se T ~ Gama(2, escala=beta) e X = floor(T), entao
# P(X=x) = P(x <= T < x+1), x=0,1,2,...

# Logaritmo de P(X=x). Esta forma evita subtracoes instaveis.
log_pmf <- function(x, beta) {
  u <- 1 / beta
  a <- pgamma(u, shape = 2, scale = 1)
  b <- -expm1(-u)
  -x * u + log(a + x * u * b)
}

# O parametro otimizado e eta=log(beta), garantindo beta>0.
log_verossimilhanca <- function(eta, x) {
  beta <- exp(eta)
  sum(log_pmf(x, beta))
}

ajustar_proposto <- function(x, regiao) {
  # Valor inicial baseado aproximadamente em E(X) = 2*beta - 0.5.
  beta_inicial <- max((mean(x) + 0.5) / 2, 0.01)

  ajuste <- optim(par = log(beta_inicial), fn = log_verossimilhanca,
                  x = x, method = "BFGS",
                  control = list(fnscale = -1, maxit = 3000),
                  hessian = TRUE)

  eta_hat <- ajuste$par
  beta_hat <- exp(eta_hat)
  loglik <- log_verossimilhanca(eta_hat, x)

  # Informação observada na escala log(beta).
  info_eta <- -as.numeric(ajuste$hessian[1, 1])
  se_eta <- if (is.finite(info_eta) && info_eta > 0) 1/sqrt(info_eta) else NA_real_
  ep_beta <- beta_hat * se_eta  # Metodo delta
  ic <- exp(eta_hat + qnorm(c(0.025, 0.975)) * se_eta)

  resumo <- data.frame(regiao = regiao, n = length(x), beta_hat = beta_hat,
                       ep_beta = ep_beta, ic95_inf = ic[1], ic95_sup = ic[2],
                       logLik = loglik, AIC = -2*loglik + 2,
                       BIC = -2*loglik + log(length(x)),
                       convergencia = ajuste$convergence)

  list(resumo = resumo, beta_hat = beta_hat, dados = x)
}

# Ajustar China e Europa separadamente.
resultado_proposto <- list()
resultado_proposto$China <- ajustar_proposto(dados_china, "China")
resultado_proposto$Europa <- ajustar_proposto(dados_europa, "Europa")

resumo_proposto <- do.call(rbind, lapply(resultado_proposto, function(z) z$resumo))
print(resumo_proposto, row.names = FALSE, digits = 5)
write.csv(resumo_proposto, file.path(pasta, "resumo_proposto.csv"), row.names = FALSE)


# ------------------------------------------------------------

plotar_diagnosticos <- function(resultado, regiao) {
  x <- resultado$dados
  beta <- resultado$beta_hat
  n <- length(x)
  k <- 0:max(x)

  # 1. FDA empirica e FDA ajustada
  plot(ecdf(x), verticals = TRUE, do.points = FALSE,
       main = paste("FDA -", regiao), xlab = "x", ylab = "F(x)",
       xlim = c(0, max(x)), col = "black", lwd = 2)
  lines(k, pgamma(k+1, shape=2, scale=beta), type="s", col="blue", lwd=2)
  legend("bottomright", c("Empirica", "Ajustada"),
         col=c("black", "blue"), lty=1, lwd=2, bty="n")

  # 2. Quantis observados versus quantis ajustados
  probs <- ppoints(n)
  q_obs <- sort(x)
  q_modelo <- floor(qgamma(probs, shape=2, scale=beta))
  plot(q_modelo, q_obs, pch=16, cex=0.6,
       main=paste("Q-Q discreto -", regiao),
       xlab="Quantis ajustados", ylab="Quantis observados")
  abline(0, 1, col="red", lty=2)

  # 3. Log-verossimilhanca em torno do maximo
  grade_beta <- seq(0.55*beta, 1.5*beta, length.out=150)
  valores <- sapply(grade_beta, function(b) sum(log_pmf(x,b)))
  plot(grade_beta, valores, type="l", lwd=2, col="blue",
       main=paste("Log-verossimilhanca -", regiao),
       xlab=expression(beta), ylab="log L(beta)")
  abline(v=beta, col="red", lty=2)
}

# Graficos 
op <- par(mfrow=c(1,3))
plotar_diagnosticos(resultado_proposto$China, "China")
par(op)
op <- par(mfrow=c(1,3))
plotar_diagnosticos(resultado_proposto$Europa, "Europa")
par(op)



# ------------------------------------------------------------
# ETAPA 2 - CONCORRENTES 
# ------------------------------------------------------------

  log_sub_exp <- function(a, b) {
    a + log(-expm1(b - a))
  }

  # Parametrizacoes internas: par = log(parametros positivos).
  modelos <- list(
    Poisson = list(
      k = 1L, nomes = "lambda",
      inicio = function(x) log(max(mean(x), 0.01)),
      logp = function(x, p) dpois(x, lambda=exp(p[1]), log=TRUE),
      F = function(k, p) ppois(k, lambda=exp(p[1]))
    ),
    Geometrica = list(
      k = 1L, nomes = "prob",
      inicio = function(x) log(max(mean(x), 0.01)),
      # p=1/(1+media); parametro interno = log(media)
      logp = function(x, p) dgeom(x, prob=1/(1+exp(p[1])), log=TRUE),
      F = function(k, p) pgeom(k, prob=1/(1+exp(p[1])))
    ),
    `Binomial negativa` = list(
      k = 2L, nomes = c("size", "mu"),
      inicio = function(x) {
        m <- mean(x); v <- var(x)
        log(c(if (v > m) max(m^2/(v-m), 0.1) else 100, max(m, 0.01)))
      },
      logp = function(x, p) dnbinom(x, size=exp(p[1]), mu=exp(p[2]), log=TRUE),
      F = function(k, p) pnbinom(k, size=exp(p[1]), mu=exp(p[2]))
    ),
    `Weibull discretizada` = list(
      k = 2L, nomes = c("shape", "scale"),
      inicio = function(x) log(c(1.5, max(mean(x), 0.1))),
      logp = function(x, p) {
        shape <- exp(p[1]); scale <- exp(p[2]);
        a <- (x/scale)^shape
        b <- ((x+1)/scale)^shape
        # P(X=x)=exp(-a)-exp(-b)
        log_sub_exp(-a, -b)
      },
      F = function(k, p) pweibull(k+1, shape=exp(p[1]), scale=exp(p[2]))
    )
  )

  ajustar_concorrente <- function(x, nome, mod) {
    objetivo <- function(par) {
      lp <- tryCatch(mod$logp(x, par), error=function(e) NA_real_)
      if (length(lp)!=length(x) || any(!is.finite(lp))) return(1e100)
      -sum(lp)
    }
    inicial <- mod$inicio(x)
    # Tentativas alternativas contribuem para evitar maximos locais/maus inicios.
    pontos <- list(inicial, inicial + log(1.2), inicial - log(1.2))
    ajustes <- lapply(pontos, function(p0) {
      tryCatch(optim(p0, objetivo, method="BFGS",
                     control=list(reltol=1e-11, maxit=3000)),
               error=function(e) NULL)
    })
    bons <- Filter(function(z) !is.null(z) && is.finite(z$value) &&
                     z$value < 1e99, ajustes)
    if (!length(bons)) stop("Nao foi possivel ajustar: ", nome)
    melhor <- bons[[which.min(vapply(bons, function(z) z$value, numeric(1)))]]
    if (melhor$convergence != 0)
      warning(sprintf("Verifique a convergencia do modelo %s (codigo %s)",
                      nome, melhor$convergence))
    est <- exp(melhor$par)
    if (nome=="Geometrica") est <- 1/(1+est)
    names(est) <- mod$nomes
    ll <- -melhor$value
    list(nome=nome, par=melhor$par, estimativas=est, logLik=ll,
         AIC=-2*ll+2*mod$k, BIC=-2*ll+log(length(x))*mod$k,
         convergencia=melhor$convergence, F=mod$F)
  }

  tabela_concorrentes <- list()
  fda_concorrentes <- list()
  for (regiao in names(resultado_proposto)) {
    x <- resultado_proposto[[regiao]]$dados
    ajustes <- lapply(names(modelos), function(nome)
      ajustar_concorrente(x, nome, modelos[[nome]]))
    names(ajustes) <- names(modelos)


    r <- resultado_proposto[[regiao]]$resumo
    tab <- data.frame(conjunto=regiao, modelo="Gama(2) discretizada",
                      parametros=paste0("beta=", round(r$beta_hat,5)),
                      k=1, logLik=r$logLik, AIC=r$AIC, BIC=r$BIC,
                      convergencia=r$convergencia)
    for (nome in names(ajustes)) {
      a <- ajustes[[nome]]
      tab <- rbind(tab, data.frame(conjunto=regiao, modelo=nome,
                   parametros=paste(sprintf("%s=%.5f", names(a$estimativas),
                                             a$estimativas), collapse="; "),
                   k=length(a$estimativas), logLik=a$logLik,
                   AIC=a$AIC, BIC=a$BIC, convergencia=a$convergencia))
    }
    tab <- tab[order(tab$AIC), ]
    rownames(tab) <- NULL
    tabela_concorrentes[[regiao]] <- tab
    fda_concorrentes[[regiao]] <- ajustes
    print(tab, row.names=FALSE, digits=5)
  }



  # FDA empirica e modelos concorrentes 
  for (regiao in names(resultado_proposto)) {
    x <- resultado_proposto[[regiao]]$dados
    k <- 0:max(x)
    plot(ecdf(x), verticals=TRUE, do.points=FALSE, col="black", lwd=2,
         main=paste("Comparacao -",regiao), xlab="x", ylab="F(x)")
    lines(k, pgamma(k+1, shape=2,
                    scale=resultado_proposto[[regiao]]$beta_hat),
          type="s", col=2, lwd=2)
    concorrentes <- fda_concorrentes[[regiao]]
    for (j in seq_along(concorrentes)) {
      a <- concorrentes[[j]]
      lines(k, a$F(k,a$par), type="s", col=j+2, lwd=2)
    }
    legend("bottomright", c("Empirica", "Gama(2) discretizada", names(concorrentes)),
           col=seq_len(length(concorrentes)+2), lty=1, lwd=2, cex=0.7, bty="n")
  }


