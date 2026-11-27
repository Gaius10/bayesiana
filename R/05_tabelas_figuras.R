# Tabelas LaTeX (relatorio/tabelas/*.tex) e figuras do estudo de simulacao
# A parte bayesiana so e gerada quando resultados/sim_bayes.csv existir
source("R/modelo.R")
dir.create("relatorio/tabelas", recursive = TRUE, showWarnings = FALSE)

# Escreve as linhas de um tabular (virgula decimal, sinal de menos, inteiros sem casas)
tabular <- function(df, arq, dig = 4, fmt = "f") {
  f <- function(v) {
    if (!is.numeric(v)) return(v)
    s <- if (all(v == round(v))) as.character(v) else formatC(v, digits = dig, format = fmt, flag = if (fmt == "fg") "#" else "")
    sub("^-", "$-$", gsub(".", ",", sub("\\.$", "", trimws(s)), fixed = TRUE))
  }
  writeLines(paste(do.call(paste, c(lapply(df, f), sep = " & ")), "\\\\"), file.path("relatorio/tabelas", arq))
}

# Simulacao
cl <- read.csv("resultados/sim_classica.csv")
tabular(cl[, c("beta", "n", "media", "vies", "eqm", "cobertura", "comprimento")], "sim_classica.tex")
tem_bayes <- file.exists("resultados/sim_bayes.csv")
if (tem_bayes) {
  ba <- read.csv("resultados/sim_bayes.csv")
  tabular(ba[, c("beta", "n", "media", "vies", "eqm", "cobertura", "comprimento")], "sim_bayes.tex")
  tabular(ba[, c("beta", "n", "rhat_medio", "rhat_max", "prop_rhat_101", "ess_medio", "ess_min")], "sim_bayes_conv.tex", 3)
  tabular(ba[, c("beta", "n", "rhat_max", "prop_rhat_101", "ess_min")], "slide_conv.tex", 3)
}

# Aplicacao
ac <- read.csv("resultados/aplic_classica.csv")
tabular(ac[, c("regiao", "modelo", "parametro", "estimativa", "ep", "inf", "sup")], "aplic_classica.tex", fmt = "fg")
ab <- read.csv("resultados/aplic_bayes_resumo.csv", check.names = FALSE)
tabular(ab[, c("regiao", "modelo", "parametro", "media", "dp", "2.5%", "97.5%", "rhat")], "aplic_bayes.tex", fmt = "fg")
af <- read.csv("resultados/aplic_tabela_final.csv")
tabular(af[, c("regiao", "modelo", "AIC", "BIC", "DIC", "WAIC", "LPML")], "aplic_final.tex", 2)

# Figura: fmp e fda da DMEx para os valores de beta da simulacao
pdf("figuras/fmp_dmex.pdf", width = 10, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
k <- 0:30
matplot(k, sapply(c(2, 4, 5), ddmex, x = k), type = "h", lwd = 2, lty = 1, col = c(1, 2, 4), xlab = "x",
        ylab = "P(X = x)", main = "Funcao de probabilidade")
legend("topright", paste("beta =", c(2, 4, 5)), col = c(1, 2, 4), lwd = 2, bty = "n")
matplot(k, sapply(c(2, 4, 5), pdmex, x = k), type = "s", lwd = 2, lty = 1, col = c(1, 2, 4), xlab = "x",
        ylab = "F(x)", main = "Funcao de distribuicao")
invisible(dev.off())

# Figura: vies, EQM e cobertura por n (EMV x media a posteriori)
pdf("figuras/sim_comparacao.pdf", width = 10, height = 8)
par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))
for (med in c("vies", "eqm", "cobertura")) for (b in unique(cl$beta)) {
  y <- cbind(cl[cl$beta == b, med], if (tem_bayes) ba[ba$beta == b, med])
  matplot(cl$n[cl$beta == b], y, type = "b", pch = 16:17, lty = 1, col = c(1, 4), log = "x", xlab = "n", ylab = med,
          main = paste0(med, " (beta = ", b, ")"), xaxt = "n")
  axis(1, at = cl$n[cl$beta == b])
  abline(h = c(vies = 0, eqm = 0, cobertura = 0.95)[med], lty = 3, col = 2)
  if (b == 2 && med == "vies") legend("topright", c("EMV", "Bayes (Gama)"), col = c(1, 4), pch = 16:17, bty = "n")
}
invisible(dev.off())

# Figura: convergencia das cadeias guardadas (1 replica com n = 25 e n = 400 por beta)
if (tem_bayes) {
  pdf("figuras/sim_convergencia.pdf", width = 10, height = 8)
  par(mfrow = c(3, 4), mar = c(4, 4, 2.5, 1))
  for (b in c(2, 4, 5)) for (n in c(25, 400)) {
    cad <- readRDS(sprintf("resultados/bayes/cadeias_b%d_n%d.rds", b, n))
    traceplot(cad, main = paste0("beta = ", b, ", n = ", n), ylab = "beta")
    abline(h = b, col = 2, lty = 2)
    acf(cad[[1]][, 1], main = "ACF (cadeia 1)")
  }
  invisible(dev.off())
}
