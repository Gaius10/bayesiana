---
name: r-enxuto
description: Use sempre que escrever ou alterar código R neste projeto
---

# Código R enxuto (trabalho SME0809)

Siga o estilo dos códigos do professor em `materiais/codigos/`: scripts diretos, de cima para baixo.

## Regras

- **Pacotes:** base R + `rjags`/`R2jags`/`coda` (+ `parallel` na simulação). Sem tidyverse, sem frameworks.
- **Funções do modelo só em `R/modelo.R`.** Os scripts fazem `source("R/modelo.R")`. Nada duplicado: antes de escrever uma função, veja se já existe em `R/modelo.R` ou em `materiais/trabalho/` e reaproveite.
- **Vetorize.** Use `sapply`/`vapply`/operações em vetor no lugar de laços quando ficar mais claro.
- **Sem defesas desnecessárias:** nada de `tryCatch`, `stopifnot`, checagens de `is.finite` "por garantia", `print`/`cat` de depuração ou código morto. Mantenha só o que o resultado usa.
- **Comentários curtos em português**, sem acentos se o arquivo for lido no Windows com outra codificação.
- **Reprodutível:** `set.seed()` fixo no topo; sementes do JAGS (`.RNG.name`, `.RNG.seed`) fixas.
- **Saídas:** tabelas em `resultados/*.csv` ou `resultados/*.rds`, figuras em `figuras/*.pdf` (sempre `pdf()` ... `dev.off()`; gráficos lattice com `print()`), texto em `relatorio/`.
- **Caminhos relativos à raiz do projeto** (rodar com `Rscript R/xx.R` a partir da raiz).
- **Nunca altere nada em `materiais/`.**

## Ao terminar cada script

Releia e reduza o que der sem perder clareza: junte linhas repetidas, remova variáveis usadas uma vez só, elimine o que não vai para uma tabela ou figura.
