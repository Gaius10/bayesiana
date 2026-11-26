# Lindley (1958): as condições (i) e (ii) do sumário

Referência: [lindley.md](../lindley.md)

## O que é uma distribuição de Bayes

É uma **distribuição a posteriori**: o resultado de aplicar o teorema de Bayes a alguma priori $p(\theta)$:

$$
\pi(\theta \mid x) = \frac{f(x \mid \theta)\, p(\theta)}{\int f(x \mid \theta)\, p(\theta)\, d\theta}
$$

Lindley quer saber se a distribuição fiducial pode ser escrita assim. A distribuição fiducial é a proposta de Fisher para obter uma distribuição para $\theta$ sem usar priori:

$$
\phi_x(\theta) = -\frac{\partial}{\partial\theta} F(x \mid \theta)
$$

A ideia de Fisher é tratar a função de distribuição acumulada $F(x \mid \theta)$, com $x$ fixo, como se fosse uma função de distribuição em $\theta$, e derivar em relação a $\theta$.

Dizer que "a fiducial é uma distribuição de Bayes" significa que **existe alguma priori** $p(\theta)$ cuja posteriori é exatamente $\phi_x(\theta)$. Essa priori pode ser imprópria, isto é, sua integral pode não dar 1. Quando ela existe, o método de Fisher é só a inferência bayesiana com uma priori específica.

**Exemplo:** $x \sim N(\theta, 1)$.

- **Fiducial:** $F(x \mid \theta) = \Phi(x - \theta)$, e daí $\phi_x(\theta) = \varphi(x - \theta)$. Ou seja, $\theta \sim N(x, 1)$.
- **Bayes com $p(\theta) \propto 1$:** a posteriori é $\propto \varphi(x - \theta)$, também $N(x, 1)$.

As duas coincidem, então aqui a fiducial é uma distribuição de Bayes.

## O que é um parâmetro de locação

$\theta$ é um **parâmetro de locação** quando a densidade depende de $x$ e $\theta$ apenas pela diferença $x - \theta$:

$$
f(x \mid \theta) = g(x - \theta) \qquad \text{ou, equivalentemente,} \qquad F(x \mid \theta) = G(x - \theta)
$$

Mudar $\theta$ só **desloca a curva** ao longo do eixo; a forma não muda. Exemplos:

- a média $\mu$ da normal com variância conhecida;
- o centro de uma Cauchy;
- $\theta$ em uma $U(\theta, \theta + 1)$.

O contraexemplo típico é o **parâmetro de escala**, em que $f(x \mid \theta) = \frac{1}{\theta}\, g(x/\theta)$. Nesse caso mudar $\theta$ estica ou comprime a curva, como o $\theta$ da exponencial ou o $\sigma$ da normal.

## Condição (A)

Muitas vezes $\theta$ não é de locação na escala original, mas passa a ser depois de uma mudança de variáveis. A condição (A) diz que existem transformações $u = R(x)$ e $\tau = P(\theta)$ tais que

$$
F(x \mid \theta) = G(u - \tau)
$$

Ou seja, $\tau$ é parâmetro de locação para $u$.

O próprio artigo dá o exemplo do parâmetro de escala. Tome $u = \log x$ e $\tau = \log\theta$. Como $x/\theta = e^{u - \tau}$, $\tau$ vira parâmetro de locação.

## Resultado (i)

> A fiducial é uma distribuição de Bayes **se e somente se** vale (A).

Quando vale (A), a priori que reproduz a fiducial é **uniforme em $\tau$**, como no exemplo da normal acima. Para voltar à escala original, use a mudança de variável $p(\theta) = p_\tau(P(\theta)) \cdot |P'(\theta)|$, com $P' = p$.

- **Caso de escala:** se $\tau = \log\theta$ é uniforme, então $p(\theta) \propto 1/\theta$. É a frase do fim da seção 1.

Isso se liga ao curso. A priori uniforme na locação e $1/\theta$ na escala são exatamente as **prioris de Jeffreys** para esses casos. O artigo de Syversveen em [referencias/](../materiais/referencias/syversveen_1998_noninformative_priors.pdf) trata desse assunto.

## Resultado (ii): consistência

Suponha duas amostras independentes, resumidas pelas estatísticas suficientes $x$ e $y$. Há dois caminhos para chegar a uma distribuição para $\theta$:

1. **$\phi_{x,y}$:** aplicar o argumento fiducial à amostra combinada.
2. **$\psi_{x,y}$:** aplicar o argumento fiducial só a $x$, usar o resultado $\phi_x$ como **priori** e atualizar com $y$ pelo teorema de Bayes.

Fisher afirmava que a probabilidade fiducial é "do mesmo tipo" que uma probabilidade a priori. Se isso fosse verdade, os dois caminhos teriam que dar o mesmo resultado.

Lindley mostra que, havendo uma estatística suficiente única, **os dois caminhos só coincidem sob (A)**. Fora de (A), $\phi_x$ não é a posteriori de nenhuma priori. Ao usá-la como priori, o resultado passa a depender de **qual amostra você processa primeiro**. Na seção 2, o artigo dá um exemplo em que $\psi_{x,y} \neq \psi_{y,x}$.

Com o teorema de Bayes isso nunca acontece. Atualizar com $x$ e depois com $y$ dá o mesmo resultado que atualizar com as duas juntas, em qualquer ordem.

Na seção 3, o artigo mostra que, dentro da família exponencial, só a **normal** e a **gama** satisfazem (A), além das distribuições que se transformam nelas.

**Conclusão:** o argumento fiducial só é consistente nos casos em que ele é, na verdade, inferência bayesiana com priori uniforme no parâmetro de locação. Fora desses casos, ele leva a contradições.
