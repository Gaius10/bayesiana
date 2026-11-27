# Derivações do modelo DMEx

Modelo: $Y \sim \text{Gama}(2, \beta)$ (forma 2, escala $\beta$), $f_Y(y) = \frac{y}{\beta^2} e^{-y/\beta}$, $y > 0$. É a moment exponential de Dara e Ahmad (2012). A variável discreta é $X = \lfloor Y \rfloor$ (Afify et al., 2022).

Notação: $u = 1/\beta$ e $q = e^{-u}$.

## Sobrevivência, fmp e fda

Integrando por partes, a sobrevivência de $Y$ é

$$S(y) = P(Y > y) = \left(1 + \tfrac{y}{\beta}\right) e^{-y/\beta}.$$

Como $X = x \iff x \le Y < x+1$, a fmp é

$$p(x;\beta) = S(x) - S(x+1) = e^{-x/\beta}\Big[\underbrace{1 - (1+u)q}_{P(\text{Gama}(2,1) \le u)} + x\,u\,(1-q)\Big] = e^{-xu}A(x;\beta), \quad x = 0, 1, \dots$$

Essa é a forma usada no código. Ela evita subtrair dois números próximos e tem $A > 0$.

A soma telescópica dá $\sum_x p(x) = S(0) - \lim_{x\to\infty} S(x) = 1$. A fda é

$$F(x) = P(X \le x) = 1 - S(x+1) = P(Y < x+1).$$

## Esperança e variância

Para $X$ inteiro não negativo, $E(X) = \sum_{k \ge 1} P(X \ge k) = \sum_{k\ge1} S(k)$ e $E(X^2) = \sum_{k\ge1}(2k-1)S(k)$. Com $S(k) = q^k(1 + k u)$ e as séries $\sum k q^k = \frac{q}{(1-q)^2}$ e $\sum k^2 q^k = \frac{q(1+q)}{(1-q)^3}$:

$$E(X) = \frac{q}{1-q} + \frac{u\,q}{(1-q)^2},$$

$$E(X^2) = 2\left[\frac{q}{(1-q)^2} + \frac{u\,q(1+q)}{(1-q)^3}\right] - E(X), \qquad \operatorname{Var}(X) = E(X^2) - E(X)^2.$$

Para $\beta$ grande, $E(X) \approx 2\beta - \tfrac12$ e $\operatorname{Var}(X) \approx 2\beta^2 + \tfrac1{12}$. A conferência numérica está em `R/00_checagens.R`.

## Log-verossimilhança

Para uma amostra $x_1, \dots, x_n$:

$$\ell(\beta) = -\frac{1}{\beta}\sum_{i} x_i + \sum_{i} \log A(x_i;\beta).$$

A otimização é feita em $\eta = \log\beta$, ou seja, $\ell(e^\eta)$, o que dispensa a restrição $\beta > 0$.

## Função escore

Derivando $S$ em $\beta$:

$$\frac{\partial S(y)}{\partial\beta} = \frac{y^2}{\beta^3}e^{-y/\beta}.$$

Daí, $\partial p/\partial\beta = \beta^{-3}\left[x^2 e^{-xu} - (x+1)^2 e^{-(x+1)u}\right]$, e o escore é

$$U(\beta) = \frac{\partial \ell}{\partial \beta} = \frac{1}{\beta^3}\sum_{i=1}^n \frac{x_i^2 - (x_i+1)^2 e^{-1/\beta}}{A(x_i;\beta)}.$$

Na escala $\eta$, $U_\eta = \beta\,U(\beta)$. O EMV resolve $U(\hat\beta) = 0$, o que não tem forma fechada. Por isso usamos BFGS com esse gradiente. O erro padrão sai da informação observada em $\eta$ e do método delta, e o IC de Wald, calculado em $\eta$, é sempre positivo.

## Posterioris

**Priori 1:** $\beta \sim \text{Gama}(a = 0{,}1;\ b = 0{,}01)$ (forma, taxa).

$$\pi_1(\beta \mid \mathbf{x}) \propto \beta^{a-1} e^{-b\beta}\, \exp\Big(-\frac{\sum x_i}{\beta}\Big)\prod_{i=1}^n A(x_i;\beta).$$

**Priori 2:** $\pi(\beta) \propto 1/\beta$, isto é, uniforme em $\eta = \log\beta$.

$$\pi_2(\beta \mid \mathbf{x}) \propto \beta^{-1}\, \exp\Big(-\frac{\sum x_i}{\beta}\Big)\prod_{i=1}^n A(x_i;\beta).$$

Nenhuma das duas tem forma conhecida, por isso usamos MCMC (JAGS, ones trick).

A motivação da Priori 2 vem de Lindley (1958). Na Gama contínua, $\beta$ é parâmetro de escala, vale a condição (A) com $\tau = \log\beta$, e a distribuição fiducial coincide com a posteriori sob $1/\beta$. Depois da discretização isso deixa de valer exatamente.

## Por que a posteriori com $1/\beta$ é própria

A priori $1/\beta$ é imprópria: sua integral diverge tanto em $0$ quanto em $\infty$. Basta então verificar que a verossimilhança $L$ compensa nas duas pontas.

- **$\beta \to \infty$ ($u \to 0$):** $1-(1+u)q \approx u^2/2$ e $x u(1-q) \approx x u^2$. Logo $A(x;\beta) \approx u^2(x + \tfrac12)$ e $L(\beta) \approx \beta^{-2n}\prod(x_i + \tfrac12)$. Assim, $\int^\infty \beta^{-1}\beta^{-2n}\,d\beta < \infty$ para todo $n \ge 1$.
- **$\beta \to 0$ ($u \to \infty$):** $A(x;\beta) \to 1 + x u$, logo $L(\beta) \approx e^{-u\sum x_i}\prod(1 + x_i u)$.
  - Se $\sum x_i \ge 1$, $L$ vai a zero mais rápido que qualquer potência de $\beta$, e $\int_0 \beta^{-1}L\,d\beta < \infty$.
  - Se todos os $x_i = 0$, $L \to 1$ e $\int_0 \beta^{-1}d\beta = \infty$.

**Conclusão:** $\pi_2(\beta \mid \mathbf{x})$ é própria **se e somente se** $\sum x_i \ge 1$, ou seja, se pelo menos uma observação é positiva.

A chance de isso falhar é $p(0)^n$. Com $\beta = 2$, $p(0) = 1 - 1{,}5\,e^{-0{,}5} \approx 0{,}090$, e para $n = 25$ isso dá $\approx 10^{-26}$. Nos dados reais, a condição vale. No JAGS, $\eta \sim U(-10, 10)$ trunca a priori em $[e^{-10}, e^{10}]$, o que a torna própria em qualquer caso, sem efeito prático porque a verossimilhança está concentrada bem dentro desse intervalo.

A Priori 1 é própria, então $\pi_1$ é sempre própria.
