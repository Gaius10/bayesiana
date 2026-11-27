# Estudo de simulacao bayesiano (Priori 1: Gama(0.1, 0.01)) sobre as mesmas amostras
# Uso: Rscript R/03_sim_bayes.R [S]   (S < 1000 roda um piloto em resultados/piloto)
# Cada cenario e salvo em <pasta>/cen_b?_n?.rds; cenarios ja salvos sao pulados.
library(parallel)
source("R/modelo.R")
sim <- readRDS("resultados/amostras.rds")
S <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else sim$S
pasta <- if (S < sim$S) "resultados/piloto" else "resultados/bayes"
dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
log_txt <- if (S < sim$S) file.path(pasta, "log.txt") else "resultados/log_sim_bayes.txt"
registra <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", file = log_txt, append = TRUE)

n_burn <- 1500; n_iter <- 3500; thin <- 1   # ajustado pelo piloto (S = 10)

cl <- makeCluster(detectCores() - 1)
invisible(clusterEvalQ(cl, source("R/modelo.R")))
clusterExport(cl, c("n_burn", "n_iter", "thin"))
registra("inicio: S =", S, "| burn =", n_burn, "| iter =", n_iter, "| thin =", thin)

for (ib in seq_along(sim$betas)) for (iN in seq_along(sim$tamanhos)) {
  b0 <- sim$betas[ib]; n <- sim$tamanhos[iN]
  arq <- file.path(pasta, sprintf("cen_b%d_n%d.rds", b0, n))
  if (file.exists(arq)) next
  t0 <- Sys.time()
  m <- sim$amostras[[ib]][[iN]][1:S, , drop = FALSE]
  res <- parLapply(cl, 1:S, function(j, m, sem) {
    cad <- ajusta_jags(m[j, ], modelo_gama, sem + 3*j, n_burn, n_iter, thin)
    list(resumo = resumo_post(cad), cadeias = if (j == 1) cad)
  }, m = m, sem = 1e5*ib + 1e4*iN)
  if (n %in% c(25, 400)) saveRDS(res[[1]]$cadeias, file.path(pasta, sprintf("cadeias_b%d_n%d.rds", b0, n)))
  saveRDS(t(sapply(res, `[[`, "resumo")), arq)
  registra(sprintf("beta = %d | n = %3d | %.1f min", b0, n, difftime(Sys.time(), t0, units = "mins")))
}
stopCluster(cl)

# Tabela final: estimador = media a posteriori; IC de credibilidade 95% por quantis
tab <- do.call(rbind, lapply(sim$betas, function(b0) t(sapply(sim$tamanhos, function(n) {
  r <- readRDS(file.path(pasta, sprintf("cen_b%d_n%d.rds", b0, n)))
  c(beta = b0, n = n, media = mean(r[, "media"]), vies = mean(r[, "media"] - b0), eqm = mean((r[, "media"] - b0)^2),
    cobertura = mean(r[, "2.5%"] <= b0 & b0 <= r[, "97.5%"]), comprimento = mean(r[, "97.5%"] - r[, "2.5%"]),
    dp_post = mean(r[, "dp"]), rhat_medio = mean(r[, "rhat"]), rhat_max = max(r[, "rhat"]),
    prop_rhat_101 = mean(r[, "rhat"] < 1.01), ess_medio = mean(r[, "ess"]), ess_min = min(r[, "ess"]))
}))))
write.csv(tab, if (S < sim$S) file.path(pasta, "sim_bayes_piloto.csv") else "resultados/sim_bayes.csv", row.names = FALSE)
registra("fim")
