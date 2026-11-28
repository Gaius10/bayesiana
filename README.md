# Trabalho de SME0809 – Inferência Bayesiana

ICMC-USP, 2º semestre de 2026. Prof. Adriano K. Suzuki.

O ponto de partida é o artigo de Afify et al. (2022), que propõe a **moment exponential discreta (DMEx)** e a ajusta a dados de mortes diárias por COVID-19 na China e na Europa. O artigo só usa métodos clássicos (máxima verossimilhança e mais seis estimadores). Neste trabalho refizemos a parte clássica e acrescentamos a parte bayesiana, que não existia: duas prioris, estudo de simulação comparando EMV e média a posteriori, e a aplicação com DIC, WAIC e LPML.

O modelo é bem simples: se Y ~ Gama(forma 2, escala β), então X = ⌊Y⌋ segue a DMEx(β), com

```
P(X = x) = S(x) - S(x+1),   S(y) = (1 + y/β) exp(-y/β),   x = 0, 1, 2, ...
```

## Grupo

| Nome | Nº USP | Parte na apresentação |
|---|---|---|
| Francisco Chiarello Ferreira | 16228894 | Abertura e considerações finais |
| Caique Figueiredo Freitas | 11836121 | Modelo |
| Caio Corrêa Chaves | 15444406 | Inferência clássica |
| Paulo Eduardo De Vincenzi | 12743130 | Inferência bayesiana |
| Pedro Paulo Silva Santos | 16993206 | Simulação clássica |
| João Vitor Costa Morais Lemes | 16877300 | Simulação bayesiana e comparação |
| Antônio Carlos Cavalcante Gonçalves | 16854872 | Aplicação a dados reais |

## Organização

```
R/                  scripts (rodar da raiz do projeto, na ordem)
  modelo.R            funcoes do modelo: fmp, fda, gerador, momentos,
                      verossimilhanca, escore, EMV e modelos JAGS
  00_checagens.R      confere fmp, escore e momentos numericamente
  01_amostras.R       gera todas as amostras da simulacao (uma vez so)
  02_sim_classica.R   simulacao com o EMV
  03_sim_bayes.R      simulacao bayesiana em paralelo (salva por cenario)
  04_aplicacao.R      dados da China e da Europa, classico e bayesiano
  05_tabelas_figuras.R  tabelas em LaTeX e figuras da simulacao
resultados/         csv e rds gerados pelos scripts
figuras/            figuras em pdf
relatorio/          relatorio.tex e slides.tex (+ pdfs compilados)
answers/            anotacoes nossas (derivacoes, andamento)
materiais/          aulas, codigos e roteiro passados pelo professor
```

## Como rodar

Precisa de R (usamos o 4.6), JAGS 5.x e dos pacotes `rjags`, `R2jags` e `coda`. O `rjags` atual do CRAN só funciona com a JAGS 5. Se ela não estiver instalada no lugar padrão, defina a variável `JAGS_HOME` apontando para a pasta da JAGS.

```
Rscript R/00_checagens.R
Rscript R/01_amostras.R
Rscript R/02_sim_classica.R      # alguns segundos
Rscript R/03_sim_bayes.R         # ~1 h com 6 nucleos
Rscript R/04_aplicacao.R         # ~2 min
Rscript R/05_tabelas_figuras.R
```

Para acompanhar a simulação bayesiana, olhe `resultados/log_sim_bayes.txt`. Se ela parar no meio, é só rodar de novo: os cenários que já foram salvos em `resultados/bayes/` são pulados. Com `Rscript R/03_sim_bayes.R 10` o script roda um piloto com 10 réplicas por cenário.

Depois é só compilar o relatório e os slides (duas vezes, por causa das referências):

```
cd relatorio
pdflatex relatorio.tex && pdflatex relatorio.tex
pdflatex slides.tex && pdflatex slides.tex
```

## Escolhas principais

- Cenários: β ∈ {2, 4, 5}, n ∈ {25, 50, 100, 200, 400}, S = 1000 réplicas. As amostras são as mesmas para os dois métodos.
- Priori 1: β ~ Gama(0,1; 0,01). Priori 2: π(β) ∝ 1/β, que no JAGS vira log β ~ U(−10, 10). A segunda vem do resultado de Lindley (1958) para parâmetros de escala. A simulação usa só a Priori 1; a aplicação usa as duas.
- A verossimilhança entra no JAGS pelo ones trick. Usamos 3 cadeias e sementes fixas, e avaliamos a convergência pelo R̂ e pelo ESS.
- Concorrentes na aplicação: Poisson, geométrica, binomial negativa e Weibull discreta (os mesmos dos códigos do professor).

## Resultados em resumo

- **Simulação:**
  - O EMV praticamente não tem viés, o EQM fica perto de β²/(2n) e a cobertura fica perto de 95%.
  - A média a posteriori tem um viés positivo pequeno quando n é pequeno (até +0,10 com n = 25), por causa da assimetria da posteriori. A partir de n ≈ 100 as duas abordagens dão praticamente o mesmo resultado.
  - Todas as 15 mil réplicas bayesianas convergiram (R̂ < 1,01).
- **Aplicação:**
  - Na China, β̂ = 25,0 e as duas prioris dão a mesma posteriori.
  - Na Europa, β̂ = 393. Lá, a Gama(0,1; 0,01) deixa de ser vaga e puxa β para baixo.
  - Pelos critérios clássicos e bayesianos, a DMEx fica atrás da geométrica, da binomial negativa e da Weibull discreta. O motivo é que o coeficiente de variação dela é praticamente fixo, em torno de 0,71.

## Observações sobre os dados

- Na lista da China que aparece no artigo há um "3" fora de ordem, logo depois do 13. Mantivemos os dados como estão publicados. Trocando esse valor por 13, o EMV bate exatamente com o do artigo (β̂ = 25,122), então provavelmente é erro de digitação.
- Para a Europa, não conseguimos reproduzir a estimativa do artigo (409,4). Com os mesmos dados obtemos 393,2, com log-verossimilhança maior. Os AIC da tabela do artigo também não conferem com −2 logL + 2k.

## Referência principal

Afify, A. Z.; Ahsan-ul-Haq, M.; Aljohani, H. M.; Alghamdi, A. S.; Babar, A.; Gómez, H. W. (2022). A new one-parameter discrete exponential distribution: properties, inference, and applications to COVID-19 data. *Journal of King Saud University – Science*, 34(6), 102199.

As demais referências estão no relatório.
