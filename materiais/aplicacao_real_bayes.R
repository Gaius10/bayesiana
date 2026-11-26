# ============================================================================
# ANALISE BAYESIANA = GAMA(2) DISCRETIZADA
# ============================================================================

library(rjags)
library(coda)
set.seed(20261028)

comparar_modelos <- FALSE

# Para uma verificacao rapida, diminua n_iter e n_burn.
n_cadeias <- 3
n_adapt <- 1000
n_burn <- 5000
n_iter <- 15000  # iteracoes apos o burn-in, por cadeia
thin <- 5


# 1. DADOS -------------------------------------------------------------------
china <- c(3,3,4,5,5,6,6,7,7,7,8,8,9,10,11,11,13,
           3,14,15,16,17,22,22,24,26,26,27,28,29,30,31,31,35,
           38,38,42,43,44,45,46,47,52,57,64,65,71,73,73,86,89,
           97,97,97,98,105,108,109,114,118,121,136,142,143,146,150)
europa <- c(6,18,28,29,44,47,55,116,118,129,150,
            184,219,236,237,336,421,434,612,648,706,838,
            1129,1393,1540,1941,2175,2278,2667,2803,2824)
conjuntos <- list(China = china, Europa = europa)

# 2. MODELO PROPOSTO ---------------------------------------------------------
# X = floor(Y), Y ~ Gamma(forma=2, escala=beta1).
# P(X=x) = S_Y(x) - S_Y(x+1), x=0,1,2,...
# Priori original: beta1 ~ Gamma(0.1, 0.01) [forma, taxa].
# No ones trick: ones[i]=1 e P(ones[i]=1|beta1)=P(X=x[i]|beta1).
# Como P(X=x|beta1) <= 1, usamos C=1.
modelo_proposto <- '
model {
  for (i in 1:N) {
    s0[i] <- exp(-x[i]/beta1) * (1 + x[i]/beta1)
    s1[i] <- exp(-(x[i]+1)/beta1) * (1 + (x[i]+1)/beta1)
    prob[i] <- s0[i] - s1[i]
    ones[i] ~ dbern(prob[i])
  }
  beta1 ~ dgamma(0.1, 0.01)
}'

# 3. FUNCOES SIMPLES ---------------------------------------------------------
# Funcao de log-probabilidade do modelo proposto; versao estavel para usar em R.
logp_gama <- function(x, beta1) {
  u <- 1/beta1
  -x*u + log(pgamma(u, shape=2, rate=1) + x*u*(-expm1(-u)))
}

# Log(media(exp(z))) sem problemas de underflow/overflow.
log_media_exp <- function(z) {
  a <- max(z)
  a + log(mean(exp(z-a)))
}


calcular_criterios <- function(ll, ll_na_media) {
  dev <- -2 * rowSums(ll)
  Dbarra <- mean(dev)
  Dtheta <- -2 * sum(ll_na_media)
  pD <- Dbarra - Dtheta
  DIC <- Dbarra + pD
  logCPO <- -apply(-ll, 2, log_media_exp)
  LPML <- sum(logCPO)
  lppd <- sum(apply(ll, 2, log_media_exp))
  pWAIC <- sum(apply(ll, 2, var))
  WAIC <- -2 * (lppd - pWAIC)
  list(DIC=DIC, LPML=LPML, WAIC=WAIC, logCPO=logCPO)
}


resultado_proposto <- list()
resumo_proposto <- data.frame()
criterios_proposto <- data.frame()

for (regiao in names(conjuntos)) {
  cat('\n===== MODELO PROPOSTO:', regiao, '=====\n')
  x <- conjuntos[[regiao]]
  N <- length(x)
  ones <- rep(1, N)
  dados <- list(N=N, x=x, ones=ones)

  # Tres cadeias, com valores iniciais distintos.
  inicio <- max(0.1, (mean(x)+0.5)/2)
  iniciais <- list(list(beta1=0.8*inicio),
                   list(beta1=1.0*inicio),
                   list(beta1=1.2*inicio))

  con <- textConnection(modelo_proposto)
  jm <- jags.model(con, data=dados, inits=iniciais,
                   n.chains=n_cadeias, n.adapt=n_adapt)
  close(con)

  update(jm, n_burn)
  cadeias <- coda.samples(jm, variable.names='beta1',
                          n.iter=n_iter, thin=thin)
  amostras <- as.numeric(as.matrix(cadeias)[, 1])

  # Resumo posterior
  media_beta <- mean(amostras)
  dp_beta <- sd(amostras)
  IC <- quantile(amostras, c(0.025, 0.975))
  rhat <- gelman.diag(cadeias, autoburnin=FALSE)$psrf[1,1]
  ess <- as.numeric(effectiveSize(cadeias)[1])
  cat('Media posterior:', round(media_beta,4), '\n')
  cat('DP posterior:', round(dp_beta,4), '\n')
  cat('IC 95%:', round(IC[1],4), 'a', round(IC[2],4), '\n')
  cat('Rhat:', round(rhat,4), ' | ESS:', round(ess,1), '\n')

  resumo_proposto <- rbind(resumo_proposto,
      data.frame(regiao=regiao, N=N, media=media_beta, dp=dp_beta,
                 IC_inf=IC[1], IC_sup=IC[2], Rhat=rhat, ESS=ess,
                 row.names=NULL))

  # GRAFICOS DE DIAGNOSTICO: seis paginas em um unico arquivo por regiao.

  traceplot(cadeias, main=paste('Traceplot -',regiao))
  densplot(cadeias, main=paste('Densidade posterior -',regiao))
  acfplot(cadeias, main=paste('Autocorrelacao (ACF) -',regiao))
  plot(amostras, type='l', col='steelblue',
       main=paste('Amostras combinadas (apenas visualizacao) -',regiao),
       xlab='Indice da amostra',ylab='beta1')
  hist(amostras, breaks=30, freq=FALSE, col='gray85',border='white',
       main=paste('Histograma posterior -',regiao),xlab='beta1')
  abline(v=media_beta, col='red', lwd=2)
  dev.off()

  # AJUSTE AOS DADOS: F(x)=P(X<=x)=P(Y<x+1)
  xx <- 0:max(x)
  F_emp <- ecdf(x)(xx)
  indice <- unique(round(seq(1,length(amostras),length.out=min(1000,length(amostras)))))
  F_mcmc <- vapply(indice, function(k)
    pgamma(xx+1, shape=2, scale=amostras[k]), numeric(length(xx)))
  F_media <- rowMeans(F_mcmc)
  F_inf <- apply(F_mcmc, 1, quantile, probs=0.025)
  F_sup <- apply(F_mcmc, 1, quantile, probs=0.975)


  plot(xx,F_emp,type='s',lwd=2,col='black',ylim=c(0,1),
       main=paste('FDA empirica e preditiva posterior -',regiao),
       xlab='x',ylab='F(x)')
  lines(xx,F_media,type='s',col='blue',lwd=2)
  lines(xx,F_inf,type='s',col='blue',lty=2)
  lines(xx,F_sup,type='s',col='blue',lty=2)
  legend('bottomright',c('Empirica','Preditiva posterior media','IC pontual 95%'),
         col=c('black','blue','blue'),lty=c(1,1,2),bty='n')
  dev.off()

  # Criterios
  ll <- vapply(amostras, function(b) logp_gama(x,b), numeric(N))
  ll <- t(ll)
  crit <- calcular_criterios(ll, logp_gama(x,media_beta))
  criterios_proposto <- rbind(criterios_proposto,
     data.frame(regiao=regiao, modelo='Gama2_discreta',
                DIC=crit$DIC, LPML=crit$LPML, WAIC=crit$WAIC))

  resultado_proposto[[regiao]] <- list(cadeias=cadeias,
                                      amostras=amostras,
                                      fd=data.frame(x=xx,F_emp=F_emp,F_media=F_media))
}

cat('\n=== RESUMO DO MODELO PROPOSTO ===\n')
print(round(resumo_proposto[, -1], 4))
print(criterios_proposto)

