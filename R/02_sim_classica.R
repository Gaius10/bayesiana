# Estudo de simulacao classico (EMV) sobre as amostras de R/01_amostras.R
source("R/modelo.R")
sim <- readRDS("resultados/amostras.rds")

# EMV de cada replica: matriz S x 6 (beta, ep, inf, sup, logLik, conv) por cenario
est <- lapply(sim$amostras, lapply, function(m) t(apply(m, 1, emv)))
saveRDS(est, "resultados/emv_replicas.rds")

medidas <- function(e, b0) {
  ok <- e[, "conv"] == 1
  c(media = mean(e[ok, "beta"]), vies = mean(e[ok, "beta"] - b0), eqm = mean((e[ok, "beta"] - b0)^2),
    cobertura = mean(e[ok, "inf"] <= b0 & b0 <= e[ok, "sup"]), comprimento = mean(e[ok, "sup"] - e[ok, "inf"]),
    ep_medio = mean(e[ok, "ep"]), prop_conv = mean(ok))
}
tab <- do.call(rbind, lapply(seq_along(sim$betas), function(i)
  cbind(beta = sim$betas[i], n = sim$tamanhos, t(sapply(est[[i]], medidas, b0 = sim$betas[i])))))
write.csv(tab, "resultados/sim_classica.csv", row.names = FALSE)
print(round(tab, 4))
