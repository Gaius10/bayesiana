# Gera uma unica vez todas as amostras do estudo de simulacao
# amostras[["b2"]][["n25"]] e uma matriz S x n (uma replica por linha)
source("R/modelo.R")
set.seed(20261127)
betas <- c(2, 4, 5)
tamanhos <- c(25, 50, 100, 200, 400)
S <- 1000

amostras <- lapply(setNames(betas, paste0("b", betas)), function(b)
  lapply(setNames(tamanhos, paste0("n", tamanhos)), function(n) matrix(rdmex(S*n, b), S, n)))

dir.create("resultados", showWarnings = FALSE)
saveRDS(list(betas = betas, tamanhos = tamanhos, S = S, amostras = amostras), "resultados/amostras.rds")
