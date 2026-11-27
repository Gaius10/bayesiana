# Progresso do trabalho (SME0809)

Atualizado em 2026-11-27, 16h. **Todas as etapas foram concluídas.** Falta só o que está em [A CONFIRMAR], no fim.

## Referência do modelo (D6)

- **Artigo base (só inferência clássica):** Afify, A. Z.; Ahsan-ul-Haq, M.; Aljohani, H. M.; Alghamdi, A. S.; Babar, A.; Gómez, H. W. (2022). *A new one-parameter discrete exponential distribution: Properties, inference, and applications to COVID-19 data*. Journal of King Saud University – Science, 34(6), 102199. DOI: 10.1016/j.jksus.2022.102199.
- **Modelo:** DMEx: X = floor(Y), Y ~ Gama(forma 2, escala β), que é a moment exponential de Dara e Ahmad (2012).
- **O que o artigo faz:** 7 estimadores clássicos, simulação e aplicação à COVID-19 (China, n = 66; Europa, n = 31). Não tem parte bayesiana.
- **Ligação com Lindley (1958):** na Gama contínua, log β é parâmetro de locação, então a fiducial coincide com a posteriori sob 1/β. Isso motiva a Priori 2.

## Etapas

- [x] 0. Skill `.claude/skills/r-enxuto/SKILL.md`.
- [x] 1. Ambiente.
  - **R:** 4.6.0, em `C:\Program Files\R\R-4.6.0\bin\Rscript.exe` (fora do PATH).
  - **JAGS:** 5.0.2, extraída sem administrador em `%LOCALAPPDATA%\Programs\JAGS\JAGS-5.0.2`. `JAGS_HOME` está no `.Renviron` e na variável de ambiente do usuário.
  - **Pacotes:** rjags 5-5, R2jags e coda.
  - **LaTeX:** MiKTeX por usuário em `%LOCALAPPDATA%\Programs\MiKTeX\miktex\bin\x64`.
- [x] 2. `R/modelo.R`. `R/00_checagens.R` confirma:
  - a fmp soma 1;
  - o escore analítico bate com a derivada numérica;
  - E(X) e Var(X) batem com 10^6 simulações.
- [x] 3. `answers/derivacoes.markdown`.
- [x] 4. `R/01_amostras.R` e `R/02_sim_classica.R`, que gera `resultados/sim_classica.csv`.
- [x] 5. `R/03_sim_bayes.R`, que gera `resultados/sim_bayes.csv`, `resultados/bayes/` e `resultados/log_sim_bayes.txt`.
  - **Piloto (S = 10):** R̂ máx. 1,007 e ESS mín. cerca de 2800.
  - **Rodada completa:** 3 cadeias, adaptação 500, burn 1000, iter 2500, thin 1. Usa 6 workers, porque com 11 faltou memória (só cerca de 0,8 GB livres).
  - **Duração:** de 15:04 a 16:02 (58 min).
  - **Convergência:** 100% das 15 000 réplicas com R̂ < 1,01 (máx. 1,009), ESS mín. 3212.
- [x] 6. `R/04_aplicacao.R`. Corrige os problemas dos scripts de `materiais/trabalho/`:
  - `pasta` indefinida;
  - `dev.off()` sem `pdf()`;
  - `acfplot` sem `print()` (troquei por `acf` em base R);
  - linha sem `#`;
  - critérios dos concorrentes e tabela final que faltavam.
- [x] 7. `R/05_tabelas_figuras.R` gera `relatorio/tabelas/*.tex` e as figuras de simulação.
  - `relatorio/relatorio.pdf`: 24 páginas.
  - `relatorio/slides.pdf`: 19 quadros, um bloco por frente, com o apresentador em cada slide.
  - Para regenerar: rodar `Rscript R/05_tabelas_figuras.R` e depois `pdflatex` duas vezes em `relatorio/`.

## Principais resultados

- **Simulação clássica:**
  - |viés| ≤ 0,03;
  - n·EQM ≈ β²/2, que é a variância assintótica da Gama de forma 2;
  - PC entre 93,6% e 96,3%.
- **Simulação bayesiana (Priori 1):**
  - viés positivo de +0,03 a +0,10 com n = 25, explicado pela assimetria da posteriori (≈ β/(2n));
  - EQM até 7% maior com n = 25 e igual ao do EMV para n ≥ 200;
  - PC entre 93,8% e 96,1%.
- **Aplicação:**
  - **China:** β̂ = 25,05 (EP 2,18); as posterioris com as duas prioris coincidem.
  - **Europa:** β̂ = 393,2. A Gama(0,1; 0,01) puxa a posteriori para 376,8, contra 399,6 com 1/β. A priori "vaga" não é vaga nessa escala.
  - **Comparação de modelos:** a DMEx perde para a geométrica, a binomial negativa e a Weibull discreta em AIC, BIC, DIC, WAIC e LPML, porque seu CV é fixo em ≈ 0,71, contra 0,89 e 1,18 nos dados. O artigo usou outros concorrentes.

## Pontos a confirmar

- [A CONFIRMAR] O "3" fora de ordem nos dados da China, logo depois do 13, também aparece no artigo. Mantido como está, mas há forte evidência de que seja 13:
  - com 13, o EMV dá β̂ = 25,122 e logL = −330,515, exatamente como no artigo;
  - a Poisson do artigo também implica soma 10 unidades maior.
- [A CONFIRMAR] Os números do artigo para a Europa são inconsistentes:
  - o artigo dá β̂ = 409,36 e logL = −257,12; nós obtemos 393,17 e −255,76;
  - a soma implícita na Poisson do artigo é 24357, e a nossa é 24361;
  - os AIC da Tabela 6 do artigo não conferem.
- [A CONFIRMAR] A referência de Dara e Ahmad (2012), livro da Lap Lambert: conferir os dados exatos.
- [A CONFIRMAR] O push para `origin` (Gaius10/bayesiana) é recusado com 403: a conta PauloVincenzi não tem permissão de escrita. O dono do repositório precisa adicioná-la como colaboradora, ou o push precisa ser feito com outra credencial. Os commits estão só locais.
- Decisões tomadas sem consulta, que o grupo pode revisar:
  - os 4 concorrentes vieram dos scripts do professor, não dos 7 do artigo;
  - na simulação, a média a posteriori é o estimador pontual;
  - na aplicação, a FDA preditiva usa só a Priori 1;
  - os slides têm 19 quadros (título e agradecimento incluídos).
