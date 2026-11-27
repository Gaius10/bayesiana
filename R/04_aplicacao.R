# Aplicacao a dados reais: mortes diarias por COVID-19 (Afify et al., 2022)
# Base: materiais/trabalho/aplicacao_real_classica.R, aplicacao_real_bayes.R e comparacao_bayes.R
source("R/modelo.R")
set.seed(20261028)
dir.create("figuras", showWarnings = FALSE)

# China: 23/01 a 28/03/2020 (o "3" apos o 13 esta assim no artigo; ver answers/progresso.markdown)
china <- c(3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9, 10, 11, 11, 13, 3, 14, 15, 16, 17, 22, 22, 24, 26, 26, 27, 28,
           29, 30, 31, 31, 35, 38, 38, 42, 43, 44, 45, 46, 47, 52, 57, 64, 65, 71, 73, 73, 86, 89, 97, 97, 97, 98,
           105, 108, 109, 114, 118, 121, 136, 142, 143, 146, 150)
# Europa: 01/03 a 31/03/2020
europa <- c(6, 18, 28, 29, 44, 47, 55, 116, 118, 129, 150, 184, 219, 236, 237, 336, 421, 434, 612, 648, 706, 838,
            1129, 1393, 1540, 1941, 2175, 2278, 2667, 2803, 2824)
dados <- list(China = china, Europa = europa)

# ---------------------------------------------------------------------------
# 1. CLASSICO: EMV, EP, IC 95%, AIC e BIC do proposto e de 4 concorrentes
# ---------------------------------------------------------------------------
# tr leva o parametro irrestrito (otimizado) para a escala natural
log_sub_exp <- function(a, b) a + log(-expm1(b - a))
concorrentes <- list(
  Poisson = list(tr = exp, ini = function(x) log(mean(x)), nomes = "lambda",
                 logp = function(x, t) dpois(x, t[1], log = TRUE), F = function(k, t) ppois(k, t[1])),
  Geometrica = list(tr = plogis, ini = function(x) qlogis(1/(1 + mean(x))), nomes = "p",
                    logp = function(x, t) dgeom(x, t[1], log = TRUE), F = function(k, t) pgeom(k, t[1])),
  `Binomial negativa` = list(tr = exp, ini = function(x) log(c(mean(x)^2/(var(x) - mean(x)), mean(x))),
                             nomes = c("r", "mu"),
                             logp = function(x, t) dnbinom(x, size = t[1], mu = t[2], log = TRUE),
                             F = function(k, t) pnbinom(k, size = t[1], mu = t[2])),
  `Weibull discreta` = list(tr = exp, ini = function(x) log(c(1.5, mean(x))), nomes = c("forma", "escala"),
                            logp = function(x, t) log_sub_exp(-(x/t[2])^t[1], -((x + 1)/t[2])^t[1]),
                            F = function(k, t) pweibull(k + 1, t[1], t[2])))

# EP pela hessiana na escala natural; IC de Wald na escala irrestrita (respeita o suporte)
ajusta_conc <- function(x, m) {
  aj <- optim(m$ini(x), function(p) -sum(m$logp(x, m$tr(p))), method = "BFGS", hessian = TRUE,
              control = list(maxit = 1000, reltol = 1e-12))
  t <- m$tr(aj$par)
  ep_p <- sqrt(diag(solve(aj$hessian)))
  data.frame(parametro = m$nomes, estimativa = t,
             ep = sqrt(diag(solve(optimHess(t, function(t) -sum(m$logp(x, t)), control = list(parscale = t))))),
             inf = m$tr(aj$par - 1.96*ep_p), sup = m$tr(aj$par + 1.96*ep_p), logLik = -aj$value,
             k = length(t), conv = aj$convergence == 0)
}

classico <- do.call(rbind, lapply(names(dados), function(r) {
  x <- dados[[r]]
  e <- emv(x)
  tab <- rbind(data.frame(modelo = "DMEx", parametro = "beta", estimativa = e[["beta"]], ep = e[["ep"]],
                          inf = e[["inf"]], sup = e[["sup"]], logLik = e[["logLik"]], k = 1, conv = e[["conv"]] == 1),
               do.call(rbind, Map(function(nm, m) cbind(modelo = nm, ajusta_conc(x, m)), names(concorrentes), concorrentes)))
  cbind(regiao = r, tab, AIC = -2*tab$logLik + 2*tab$k, BIC = -2*tab$logLik + log(length(x))*tab$k)
}))
write.csv(classico, "resultados/aplic_classica.csv", row.names = FALSE)

# FDA empirica x ajustadas e Q-Q discreto da DMEx
pdf("figuras/ajuste_classico.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (r in names(dados)) {
  x <- dados[[r]]; k <- 0:max(x); cl <- classico[classico$regiao == r, ]
  b <- cl$estimativa[cl$modelo == "DMEx"]
  plot(ecdf(x), verticals = TRUE, do.points = FALSE, lwd = 2, main = paste("FDA -", r), xlab = "x", ylab = "F(x)")
  lines(k, pdmex(k, b), type = "s", col = 2, lwd = 2)
  for (j in seq_along(concorrentes))
    lines(k, concorrentes[[j]]$F(k, cl$estimativa[cl$modelo == names(concorrentes)[j]]), type = "s", col = j + 2)
  legend("bottomright", c("Empirica", "DMEx", names(concorrentes)), col = 1:6, lty = 1, bty = "n", cex = 0.8)
  plot(floor(qgamma(ppoints(length(x)), 2, scale = b)), sort(x), pch = 16, main = paste("Q-Q DMEx -", r),
       xlab = "Quantis ajustados", ylab = "Quantis observados")
  abline(0, 1, col = 2, lty = 2)
}
invisible(dev.off())

# ---------------------------------------------------------------------------
# 2. BAYESIANO: DMEx com as duas prioris e os 4 concorrentes (JAGS)
# ---------------------------------------------------------------------------
n_burn <- 5000; n_iter <- 15000; thin <- 5   # 3 cadeias x 3000 = 9000 amostras

modelos_bayes <- list(
  `DMEx - Gama` = list(codigo = modelo_gama, vars = "beta", ini = function(x) ini_dmex(x, modelo_gama),
                       logp = function(x, p) logp_dmex(x, p[["beta"]])),
  `DMEx - 1/beta` = list(codigo = modelo_jeffreys, vars = "beta", ini = function(x) ini_dmex(x, modelo_jeffreys),
                         logp = function(x, p) logp_dmex(x, p[["beta"]])),
  Poisson = list(codigo = 'model {
      for (i in 1:N) { x[i] ~ dpois(lambda) }
      lambda ~ dgamma(0.1, 0.01)
    }', vars = "lambda", ini = function(x) lapply(c(0.8, 1, 1.2), function(f) list(lambda = f*mean(x))),
    logp = function(x, p) dpois(x, p[["lambda"]], log = TRUE)),
  Geometrica = list(codigo = 'model {
      for (i in 1:N) { x[i] ~ dnegbin(p, 1) }
      p ~ dbeta(1, 1)
    }', vars = "p", ini = function(x) lapply(c(0.8, 1, 1.2), function(f) list(p = f/(1 + mean(x)))),
    logp = function(x, p) dgeom(x, p[["p"]], log = TRUE)),
  `Binomial negativa` = list(codigo = 'model {
      for (i in 1:N) { x[i] ~ dnegbin(r/(r + mu), r) }
      mu ~ dgamma(0.1, 0.01)
      r ~ dgamma(2, 0.1)
    }', vars = c("mu", "r"), ini = function(x) Map(function(f, r) list(mu = f*mean(x), r = r), c(0.8, 1, 1.2), c(0.7, 2, 8)),
    logp = function(x, p) dnbinom(x, size = p[["r"]], mu = p[["mu"]], log = TRUE)),
  `Weibull discreta` = list(codigo = 'model {
      for (i in 1:N) {
        prob[i] <- exp(-pow(x[i]/b, a)) - exp(-pow((x[i] + 1)/b, a))
        ones[i] ~ dbern(prob[i])
      }
      a ~ dgamma(2, 1)
      b ~ dgamma(0.1, 0.01)
    }', vars = c("a", "b"), ini = function(x) Map(function(a, f) list(a = a, b = f*mean(x)), c(0.7, 1, 1.5), c(0.8, 1, 1.2)),
    logp = function(x, p) log_sub_exp(-(x/p[["b"]])^p[["a"]], -((x + 1)/p[["b"]])^p[["a"]])))

ajustes <- lapply(dados, function(x) lapply(modelos_bayes, function(m)
  ajusta_jags(x, m$codigo, 2026, n_burn, n_iter, thin, ini = m$ini(x), vars = m$vars, n_adapt = 1000)))
saveRDS(ajustes, "resultados/aplic_cadeias.rds")

# Resumo a posteriori, R-hat e ESS de cada parametro
resumo <- do.call(rbind, lapply(names(dados), function(r) do.call(rbind, lapply(names(modelos_bayes), function(nm) {
  cad <- ajustes[[r]][[nm]]
  do.call(rbind, lapply(varnames(cad), function(v) cbind(regiao = r, modelo = nm, parametro = v,
                                                           as.data.frame(t(resumo_post(cad[, v]))))))
}))))
write.csv(resumo, "resultados/aplic_bayes_resumo.csv", row.names = FALSE)

# DIC, WAIC e LPML (log-verossimilhanca pontual em cada amostra a posteriori)
crit <- do.call(rbind, lapply(names(dados), function(r) do.call(rbind, lapply(names(modelos_bayes), function(nm) {
  x <- dados[[r]]; m <- modelos_bayes[[nm]]; a <- as.matrix(ajustes[[r]][[nm]])
  ll <- t(apply(a, 1, function(p) m$logp(x, p)))
  data.frame(regiao = r, modelo = nm, t(criterios(ll, m$logp(x, colMeans(a)))))
}))))
write.csv(crit, "resultados/aplic_bayes_criterios.csv", row.names = FALSE)

# Diagnosticos da DMEx: traco, ACF e densidade (uma linha por regiao x priori)
pdf("figuras/diag_mcmc.pdf", width = 10, height = 10)
par(mfrow = c(4, 3), mar = c(4, 4, 2.5, 1))
for (r in names(dados)) for (nm in names(modelos_bayes)[1:2]) {
  cad <- ajustes[[r]][[nm]]
  traceplot(cad, main = paste("Traco -", r, "-", nm), ylab = expression(beta))
  acf(cad[[1]][, 1], main = paste("ACF (cadeia 1) -", r))
  densplot(cad, main = paste("Densidade -", r, "-", nm), xlab = expression(beta))
}
invisible(dev.off())

# Posterioris de beta (duas prioris) e FDA preditiva a posteriori
pdf("figuras/bayes_aplicacao.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (r in names(dados)) {
  b <- lapply(ajustes[[r]][1:2], function(cad) unlist(cad[, "beta"]))
  d <- lapply(b, density)
  plot(d[[1]], lwd = 2, main = paste("Posteriori de beta -", r), xlab = expression(beta),
       ylim = c(0, max(sapply(d, function(z) max(z$y)))))
  lines(d[[2]], col = 4, lwd = 2, lty = 2)
  abline(v = classico$estimativa[classico$regiao == r & classico$modelo == "DMEx"], col = 2, lty = 3)
  legend("topright", c("Priori Gama(0,1; 0,01)", "Priori 1/beta", "EMV"), col = c(1, 4, 2), lty = c(1, 2, 3), bty = "n")
}
for (r in names(dados)) {
  x <- dados[[r]]; k <- 0:max(x)
  Fp <- sapply(unlist(ajustes[[r]][[1]][, "beta"])[seq(1, 9000, by = 9)], function(bb) pdmex(k, bb))
  plot(ecdf(x), verticals = TRUE, do.points = FALSE, lwd = 2, main = paste("FDA preditiva -", r), xlab = "x", ylab = "F(x)")
  lines(k, rowMeans(Fp), type = "s", col = 4, lwd = 2)
  matlines(k, t(apply(Fp, 1, quantile, c(0.025, 0.975))), type = "s", col = 4, lty = 2)
  legend("bottomright", c("Empirica", "Preditiva (media)", "Faixa 95%"), col = c(1, 4, 4), lty = c(1, 1, 2), bty = "n")
}
invisible(dev.off())

# Tabela final: criterios classicos (AIC, BIC) e bayesianos (DIC, WAIC, LPML) lado a lado
crit$base <- sub(" - .*", "", crit$modelo)
final <- merge(crit, unique(classico[, c("regiao", "modelo", "logLik", "AIC", "BIC")]),
               by.x = c("regiao", "base"), by.y = c("regiao", "modelo"))
final <- final[order(final$regiao, final$DIC), c("regiao", "modelo", "logLik", "AIC", "BIC", "DIC", "WAIC", "LPML")]
write.csv(final, "resultados/aplic_tabela_final.csv", row.names = FALSE)
print(final, digits = 6)
