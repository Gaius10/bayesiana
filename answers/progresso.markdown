# Progresso do trabalho (SME0809)

Atualizado em 2026-11-27.

## Referência do modelo (D6)

- **Artigo base (só inferência clássica):** Afify, A. Z.; Ahsan-ul-Haq, M.; Aljohani, H. M.; Alghamdi, A. S.; Babar, A.; Gómez, H. W. (2022). *A new one-parameter discrete exponential distribution: Properties, inference, and applications to COVID-19 data*. Journal of King Saud University – Science, 34(6), 102199. DOI: 10.1016/j.jksus.2022.102199.
- **Modelo:** discrete moment exponential (DMEx): X = floor(Y), Y ~ Gama(forma 2, escala β), que é a moment exponential de Dara e Ahmad (2012). Discretização pela sobrevivência: P(X = x) = S(x) − S(x+1), com S(x) = (1 + x/β) e^{−x/β}.
- **O que o artigo faz:** 7 estimadores clássicos (EMV, Anderson-Darling, Anderson-Darling da cauda direita, Cramér-von Mises, mínimos quadrados, mínimos quadrados ponderados, percentis), simulação e duas aplicações. **Não tem parte bayesiana.** É ela que o trabalho acrescenta.
- **Dados:** mortes diárias por COVID-19.
  - China: 23/01 a 28/03/2020, worldometers, n = 66.
  - Europa: 01/03 a 31/03/2020, OMS, n = 31.
- **EMV do artigo (Tabela 6):**
  - China: β̂ = 25,1237 (EP 2,187), logL = −330,515.
  - Europa: β̂ = 409,362 (EP 52,09), logL = −257,122.
  - Na tabela da Europa, o AIC dos concorrentes não confere com −2 logL + 2k. Vamos refazer tudo e conferir.
- **Ligação com Lindley (1958), que está na pasta:** na Gama com escala β, que é contínua, vale a condição (A) com τ = log β. A posteriori fiducial coincide com a de Bayes sob π(β) ∝ 1/β. Isso motiva a Priori 2. Ao discretizar, (A) deixa de valer exatamente, e por isso vale a pena comparar as duas prioris.
- **Referência Dara e Ahmad (2012):** [A CONFIRMAR] os dados exatos (livro da Lap Lambert, *Recent Advances in Moment Distribution and their Hazard Rates*).

## Etapas

- [x] 0. Skill `.claude/skills/r-enxuto/SKILL.md` criada.
- [ ] 1. Ambiente: **BLOQUEADO**. Há R 4.6.0 em `C:\Program Files\R\R-4.6.0` (fora do PATH). Faltam JAGS, rjags, R2jags, coda e LaTeX (pdflatex + beamer).
- [ ] 2. R/modelo.R
- [ ] 3. answers/derivacoes.markdown
- [ ] 4. R/01_amostras.R e R/02_sim_classica.R
- [ ] 5. R/03_sim_bayes.R (piloto + completa em segundo plano)
- [ ] 6. R/04_aplicacao.R
- [ ] 7. relatorio/relatorio.tex e relatorio/slides.tex

## Pontos a confirmar

- O "3" fora de ordem nos dados da China, logo depois do 13, **também aparece na lista do artigo original**. Mantido como está. Confirmar com o grupo ou com o professor se é erro de digitação do artigo.
- Os nº USP de todos os integrantes foram informados, então não há pendência aqui.
