# Funcoes do modelo DMEx (moment exponential discreta)
# X = floor(Y), Y ~ Gama(forma = 2, escala = beta)
# S(y) = P(Y > y) = (1 + y/beta) exp(-y/beta)
# P(X = x) = S(x) - S(x+1), x = 0, 1, 2, ...

library(rjags)
library(coda)

# Parte estavel da fmp: P(X=x) = exp(-x/beta) * A(x, beta)
# A = P(Gama(2,1) <= 1/beta) + (x/beta) (1 - exp(-1/beta))
A_dmex <- function(x, beta) {
  u <- 1/beta
  pgamma(u, shape = 2) - x*u*expm1(-u)
}

logp_dmex <- function(x, beta) -x/beta + log(A_dmex(x, beta))
ddmex <- function(x, beta) exp(logp_dmex(x, beta))
pdmex <- function(x, beta) pgamma(x + 1, shape = 2, scale = beta)   # F(x) = P(Y < x+1)
rdmex <- function(n, beta) floor(rgamma(n, shape = 2, scale = beta))

# Momentos em forma fechada: E(X) = sum_{k>=1} S(k), E(X^2) = sum_{k>=1} (2k-1) S(k)
media_dmex <- function(beta) {
  q <- exp(-1/beta)
  q/(1 - q) + q/(beta*(1 - q)^2)
}
var_dmex <- function(beta) {
  q <- exp(-1/beta)
  m <- media_dmex(beta)
  2*(q/(1 - q)^2 + q*(1 + q)/(beta*(1 - q)^3)) - m - m^2
}

# Log-verossimilhanca em eta = log(beta) e escore analitico em beta
loglik <- function(eta, x) sum(logp_dmex(x, exp(eta)))
escore <- function(beta, x) sum((x^2 - (x + 1)^2*exp(-1/beta))/A_dmex(x, beta))/beta^3

# EMV via BFGS em eta (gradiente = beta * U(beta)); IC de Wald em eta
emv <- function(x) {
  aj <- optim(log((mean(x) + 0.5)/2), loglik, function(eta, x) exp(eta)*escore(exp(eta), x),
              x = x, method = "BFGS", control = list(fnscale = -1), hessian = TRUE)
  ep_eta <- 1/sqrt(-aj$hessian[1, 1])
  b <- exp(aj$par)
  c(beta = b, ep = b*ep_eta, inf = exp(aj$par - 1.96*ep_eta), sup = exp(aj$par + 1.96*ep_eta),
    logLik = aj$value, conv = aj$convergence == 0 & aj$hessian[1, 1] < 0)
}

# Modelos JAGS (ones trick, P(X=x) <= 1, logo C = 1), na forma estavel de P(X=x)
verossim_jags <- '
  for (i in 1:N) {
    p[i] <- exp(-x[i]/beta) * (1 - exp(-1/beta)*(1 + 1/beta) + x[i]/beta*(1 - exp(-1/beta)))
    ones[i] ~ dbern(p[i])
  }'
modelo_gama <- paste('model {', verossim_jags, '
  beta ~ dgamma(0.1, 0.01)
}')
modelo_jeffreys <- paste('model {', verossim_jags, '
  eta ~ dunif(-10, 10)
  beta <- exp(eta)
}')

# Iniciais da DMEx para as 3 cadeias: espalhadas em torno de (media + 0.5)/2
ini_dmex <- function(x, modelo) {
  b0 <- (mean(x) + 0.5)/2*c(0.5, 1, 2)
  lapply(b0, function(b) if (grepl("dunif", modelo)) list(eta = log(b)) else list(beta = b))
}

# Ajuste MCMC com 3 cadeias e sementes fixas (o ones so entra se o modelo usar)
ajusta_jags <- function(x, modelo, semente, n_burn, n_iter, thin, ini = ini_dmex(x, modelo),
                        vars = "beta", n_adapt = 500) {
  ini <- lapply(1:3, function(k) c(ini[[k]], .RNG.name = "base::Mersenne-Twister", .RNG.seed = semente + k))
  dados <- list(N = length(x), x = x)
  if (grepl("ones", modelo)) dados$ones <- rep(1, length(x))
  jm <- jags.model(textConnection(modelo), dados, ini, n.chains = 3, n.adapt = n_adapt, quiet = TRUE)
  update(jm, n_burn, progress.bar = "none")
  coda.samples(jm, vars, n_iter, thin = thin, progress.bar = "none")
}

# Resumo a posteriori de uma mcmc.list de beta
resumo_post <- function(cad) {
  b <- unlist(cad)
  c(media = mean(b), mediana = median(b), dp = sd(b), quantile(b, c(0.025, 0.975)),
    rhat = gelman.diag(cad, autoburnin = FALSE)$psrf[[1]], ess = effectiveSize(cad)[[1]])
}

# Criterios bayesianos a partir da matriz ll (amostras x observacoes)
log_media_exp <- function(z) max(z) + log(mean(exp(z - max(z))))
criterios <- function(ll, ll_media) {
  D <- -2*rowSums(ll)
  c(DIC = 2*mean(D) + 2*sum(ll_media), WAIC = -2*(sum(apply(ll, 2, log_media_exp)) - sum(apply(ll, 2, var))),
    LPML = -sum(apply(-ll, 2, log_media_exp)))
}
