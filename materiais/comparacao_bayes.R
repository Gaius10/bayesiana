# ================================================================
# COMPARACAO BAYESIANA 
# ================================================================

library(rjags)
library(coda)
set.seed(20261028)


china <- c(3,3,4,5,5,6,6,7,7,7,8,8,9,10,11,11,13,
           3,14,15,16,17,22,22,24,26,26,27,28,29,30,31,31,35,
           38,38,42,43,44,45,46,47,52,57,64,65,71,73,73,86,89,
           97,97,97,98,105,108,109,114,118,121,136,142,143,146,150)
europa <- c(6,18,28,29,44,47,55,116,118,129,150,
            184,219,236,237,336,421,434,612,648,706,838,
            1129,1393,1540,1941,2175,2278,2667,2803,2824)
conjuntos <- list(China=china, Europa=europa)

 PARAMETROS DA SIMULACAO MCMC
n_cadeias <- 3
n_adapt <- 1000
n_burn <- 5000
n_iter <- 15000       # por cadeia, apos o burn-in
thin <- 5

#  MODELOS CONCORRENTES
# Todos os modelos aceitam x = 0, 1, 2, ... (SEM truncamento em zero).
# As prioris abaixo sao as mesmas da versao anterior da comparacao.
modelos <- list(
  Poisson=list(
    codigo='model {
      for (i in 1:N) { x[i] ~ dpois(lambda) }
      lambda ~ dgamma(0.1, 0.01)
    }',
    parametros='lambda',
    inicio=function(x,k) list(lambda=max(0.1,mean(x))*c(0.8,1,1.2)[k]),
    logp=function(x,p) dpois(x,lambda=p[['lambda']],log=TRUE),
    fda=function(xx,p) ppois(xx,lambda=p[['lambda']])
  ),
  Geometrica=list(
    codigo='model {
      for (i in 1:N) { x[i] ~ dnegbin(p, 1) }
      p ~ dbeta(1,1)
    }',
    parametros='p',
    inicio=function(x,k) list(p=min(0.99,max(0.001,c(0.8,1,1.2)[k]/(1+mean(x))))),
    logp=function(x,p) dgeom(x,prob=p[['p']],log=TRUE),
    fda=function(xx,p) pgeom(xx,prob=p[['p']])
  ),
  Binomial_negativa=list(
    codigo='model {
      p <- r/(r+mu)
      for (i in 1:N) { x[i] ~ dnegbin(p,r) }
      mu ~ dgamma(0.1,0.01)
      r ~ dgamma(2,0.1)
    }',
    parametros=c('mu','r'),
    inicio=function(x,k) list(mu=max(0.1,mean(x))*c(0.8,1,1.2)[k],
                              r=c(0.7,2,8)[k]),
    logp=function(x,p) dnbinom(x,mu=p[['mu']],size=p[['r']],log=TRUE),
    fda=function(xx,p) pnbinom(xx,mu=p[['mu']],size=p[['r']])
  ),
  Weibull_discreta=list(
    codigo='model {
      for (i in 1:N) {
        s0[i] <- exp(-pow(x[i]/b,a))
        s1[i] <- exp(-pow((x[i]+1)/b,a))
        prob[i] <- s0[i]-s1[i]
        ones[i] ~ dbern(prob[i])
      }
      a ~ dgamma(2,1)
      b ~ dgamma(0.1,0.01)
    }',
    parametros=c('a','b'),
    inicio=function(x,k) list(a=c(0.7,1,1.5)[k],
                              b=max(1,mean(x))*c(0.8,1,1.2)[k]),
    logp=function(x,p) {
      a <- p[['a']]; b <- p[['b']]
      v <- (x/b)^a
      dif <- ((x+1)/b)^a-v
      -v+log(-expm1(-dif))
    },
    fda=function(xx,p) pweibull(xx+1,shape=p[['a']],scale=p[['b']])
  )
)

#  FUNCOES AUXILIARES PARA OS CRITERIOS
log_media_exp <- function(z) {
  m <- max(z)
  m + log(mean(exp(z-m)))
}

calcular_criterios <- function(ll, ll_media) {
  Dbarra <- mean(-2*rowSums(ll))
  Dhat <- -2*sum(ll_media)
  DIC <- 2*Dbarra-Dhat
  LPML <- sum(-apply(-ll, 2, log_media_exp))
  lppd <- sum(apply(ll, 2, log_media_exp))
  pWAIC <- sum(apply(ll, 2, var))
  WAIC <- -2*(lppd-pWAIC)
  data.frame(DIC=DIC, LPML=LPML, WAIC=WAIC)
}

#  MODELOS ALTERNATIVOS
comparacao <- criterios_proposto
resumo_alternativos <- data.frame()
fda_ajustada <- list()

for (regiao in names(conjuntos)) {
  x <- conjuntos[[regiao]]
  N <- length(x)
  xx <- sort(unique(c(0,x)))

  for (nome in names(modelos)) {
    cat('\nAjustando', nome, '-', regiao, '\n')
    mod <- modelos[[nome]]
    dados <- list(N=N, x=as.integer(x))
    # Apenas a Weibull discreta usa o ones trick nesta comparacao.
    if (nome == 'Weibull_discreta') {
      dados$ones <- rep(1L, N)
    }
    iniciais <- lapply(1:n_cadeias,function(k) mod$inicio(x,k))

    con <- textConnection(mod$codigo)
    jm <- tryCatch(jags.model(con,data=dados,inits=iniciais,
                              n.chains=n_cadeias,n.adapt=n_adapt,quiet=TRUE),
                   finally=close(con))
    update(jm,n_burn,progress.bar='none')
    cadeias <- coda.samples(jm, variable.names=mod$parametros,
                            n.iter=n_iter,thin=thin,progress.bar='none')
    mat <- as.matrix(cadeias)
    medias <- as.list(colMeans(mat))

    # Resumos posteriores e convergencia para cada parametro
    rhat <- gelman.diag(cadeias,autoburnin=FALSE,multivariate=FALSE)$psrf[,1]
    ess <- effectiveSize(cadeias)
    for (param in mod$parametros) {
      valores <- mat[,param]
      resumo_alternativos <- rbind(resumo_alternativos,
        data.frame(regiao=regiao,modelo=nome,parametro=param,
                   media=mean(valores),dp=sd(valores),
                   IC_inf=unname(quantile(valores,0.025)),
                   IC_sup=unname(quantile(valores,0.975)),
                   Rhat=unname(rhat[param]),ESS=unname(ess[param])))
    }
}}
 
