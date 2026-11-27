# Progresso do trabalho (SME0809)

Atualizado em 2026-11-27.

## Referência do modelo (D6)

- **Artigo base (só inferência clássica):** Afify, A. Z.; Ahsan-ul-Haq, M.; Aljohani, H. M.; Alghamdi, A. S.; Babar, A.; Gómez, H. W. (2022). *A new one-parameter discrete exponential distribution: Properties, inference, and applications to COVID-19 data*. Journal of King Saud University – Science, 34(6), 102199. DOI: 10.1016/j.jksus.2022.102199.
- **Modelo:** discrete moment exponential (DMEx): X = floor(Y), Y ~ Gama(forma 2, escala β), que é a moment exponential de Dara e Ahmad (2012). Discretização pela sobrevivência: P(X = x) = S(x) − S(x+1), com S(x) = (1 + x/β) e^{−x/β}.
- **O que o artigo faz:** 7 estimadores clássicos (EMV, Anderson-Darling, Anderson-Darling da cauda direita, Cramér-von Mises, mínimos quadrados, mínimos quadrados ponderados, percentis), simulação e duas aplicações. **Não tem parte bayesiana.** É ela que o trabalho acrescenta.
- **Dados:** mortes diárias por COVID-19.
  - China: 23/01 a 28/03/2020, worldometers, n = 66.
  - Europa: 01/03 a 31/03/2020, OMS, n = 31.
- **Ligação com Lindley (1958), que está na pasta:** na Gama contínua, β é parâmetro de escala e vale a condição (A) com τ = log β. A fiducial coincide com a posteriori sob π(β) ∝ 1/β, e isso motiva a Priori 2.

## Etapas

- [x] 0. Skill `.claude/skills/r-enxuto/SKILL.md`.
- [x] 1. Ambiente.
  - **R:** 4.6.0, em `C:\Program Files\R\R-4.6.0\bin\Rscript.exe` (fora do PATH).
  - **JAGS:** o rjags 5-5 exige a JAGS 5.
    - O instalador pede administrador. Foi extraído com 7-Zip em `%LOCALAPPDATA%\Programs\JAGS\JAGS-5.0.2`.
    - `JAGS_HOME` foi gravado em `C:\Users\paulo\.Renviron` **e** como variável de ambiente do usuário. Um R lançado pelo PowerShell não lê esse `.Renviron`.
    - A JAGS 4.3.2, em Program Files, não é usada.
  - **Pacotes:** rjags 5-5, R2jags 0.8-9 e coda 0.19-4.1.
  - **LaTeX:** MiKTeX 25.12, por usuário, em `%LOCALAPPDATA%\Programs\MiKTeX\miktex\bin\x64`, com AutoInstall ligado.
- [x] 2. `R/modelo.R`. `R/00_checagens.R` confirma:
  - a fmp soma 1;
  - o escore analítico bate com a derivada numérica;
  - E(X) e Var(X) em forma fechada batem com 10^6 simulações.
- [x] 3. `answers/derivacoes.markdown`.
- [x] 4. `R/01_amostras.R` e `R/02_sim_classica.R`.
  - Saídas: `resultados/amostras.rds` e `resultados/sim_classica.csv`.
  - Rodados completos em 4 s, com 100% de convergência.
- [ ] 5. `R/03_sim_bayes.R`.
  - **Piloto (S = 10):** 38 s com 11 workers, R̂ máx. 1,007 e ESS mín. cerca de 2800.
  - **Rodada completa:** 3 cadeias, burn 1000, iter 2500, thin 1.
    - Usa `min(detectCores() - 1, 6)` workers porque a máquina tem só cerca de 0,8 GB livres. Com 11 workers deu "cannot allocate vector".
    - Iniciada às 15:04:45. Ritmo de cerca de 1,3 + 0,016·n min por cenário, com término previsto por volta das 16h10.
  - **Acompanhar:** `Get-Content resultados\log_sim_bayes.txt -Wait`.
  - **Se cair:** rodar de novo `Rscript R/03_sim_bayes.R`. Os cenários já salvos em `resultados/bayes/` são pulados.
- [x] 6. `R/04_aplicacao.R`. Gera:
  - `resultados/aplic_*.csv`;
  - `figuras/ajuste_classico.pdf`, `diag_mcmc.pdf` e `bayes_aplicacao.pdf`.
- [ ] 7. `relatorio/relatorio.tex` e `relatorio/slides.tex`.

## Achados da aplicação

- **China:** β̂ = 25,05 (EP 2,18).
  - **A DMEx não é o melhor modelo** contra Poisson, geométrica, binomial negativa e Weibull discreta. A geométrica tem AIC 650,6 contra 664,9 da DMEx, e DIC, WAIC e LPML concordam.
- **Europa:** β̂ = 393,2 (EP 49,9). A DMEx perde para a binomial negativa, a Weibull discreta e a geométrica, com diferença de cerca de 36 de AIC.
- **Comparação com o artigo:** ele só usou outros concorrentes (DBur, DBH, DR, DITL, DPr, DIR, Poisson).
- **Prioris:** na China, as posterioris são praticamente iguais. Na Europa, a Gama(0,1; 0,01) puxa β para baixo: média 376,8 contra 399,6 com 1/β, porque e^{−0,01β} pesa em β ≈ 400. A priori "vaga" não é vaga nessa escala.
- **Convergência:** todas as cadeias da aplicação têm R̂ ≤ 1,001 e ESS ≥ 8300.

## Pontos a confirmar

- [A CONFIRMAR] O "3" fora de ordem nos dados da China, logo depois do 13, também aparece no artigo. Mantido como está, mas há forte evidência de que seja 13:
  - com 13, o EMV dá β̂ = 25,122 e logL = −330,515, exatamente como no artigo; com 3 dá 25,046 e −331,467;
  - a Poisson do artigo (λ = 49,742) implica soma 10 unidades maior que a nossa.
- [A CONFIRMAR] Os números do artigo para a Europa são inconsistentes:
  - o artigo dá β̂ = 409,36 e logL = −257,12; nós obtemos 393,17 e −255,76, e o β do artigo aplicado aos nossos dados dá −255,81;
  - a Poisson do artigo implica soma 24357, e a nossa é 24361;
  - os AIC dos concorrentes na Tabela 6 do artigo não conferem com −2logL + 2k.
- [A CONFIRMAR] A referência completa de Dara e Ahmad (2012).
- [A CONFIRMAR] O push para `origin` (Gaius10/bayesiana) foi recusado com 403: a conta PauloVincenzi não tem permissão de escrita. Os commits estão só locais.
