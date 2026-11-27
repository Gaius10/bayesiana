# Checagens numericas das funcoes de R/modelo.R
source("R/modelo.R")
set.seed(2026)
betas <- c(0.5, 2, 4, 5, 25, 400)

# 1. fmp soma 1 e coincide com S(x) - S(x+1)
S <- function(y, beta) (1 + y/beta)*exp(-y/beta)
soma <- sapply(betas, function(b) sum(ddmex(0:(200*b), b)))
dif <- sapply(betas, function(b) max(abs(ddmex(0:100, b) - (S(0:100, b) - S(1:101, b)))))

# 2. escore analitico x derivada numerica (diferenca central)
x <- rdmex(50, 4)
h <- 1e-5
num <- sapply(betas, function(b) (sum(logp_dmex(x, b + h)) - sum(logp_dmex(x, b - h)))/(2*h))
ana <- sapply(betas, escore, x = x)

# 3. E(X) e Var(X) em forma fechada x simulacao (10^6 valores)
sim <- sapply(betas, function(b) {z <- rdmex(1e6, b); c(mean(z), var(z))})

print(round(data.frame(beta = betas, soma_fmp = soma, max_dif_S = dif, U_analitico = ana, U_numerico = num,
                       EX = media_dmex(betas), EX_sim = sim[1, ], VarX = var_dmex(betas), VarX_sim = sim[2, ]), 6))
