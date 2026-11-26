# Fiducial Distributions and Bayes' Theorem

**D. V. Lindley**
*Statistical Laboratory, Cambridge*
[Received January, 1957]

*Journal of the Royal Statistical Society, Series B*, Vol. 20, No. 1 (1958), pp. 102–107.

---

## Summary

$x$ is a one-dimensional random variable whose distribution depends on a single parameter $\theta$. It is the purpose of this note to establish two results:

1. The necessary and sufficient condition for the fiducial distribution of $\theta$, given $x$, to be a Bayes' distribution is that there exist transformations of $x$ to $u$, and of $\theta$ to $\tau$, such that $\tau$ is a location parameter for $u$. The condition will be referred to as (A). This extends some results of Grundy's (1956).
2. If, for a random sample of any size from the distribution for $x$, there exists a single sufficient statistic for $\theta$ then the fiducial argument is inconsistent unless condition (A) obtains: and when it does, the fiducial argument is equivalent to a Bayesian argument with uniform prior distribution for $\tau$.

The note concludes with an investigation of (A) in the case of the exponential family.

---

## 1. The Necessary and Sufficient Condition for a Fiducial Distribution to be a Bayes' Distribution

Let $F(x \mid \theta)$ be a distribution function of a one-dimensional random variable $x$ for values of a real parameter $\theta$ belonging to a one-dimensional set. Then according to Fisher (1956, p. 70) the fiducial distribution for $\theta$ is given by

$$
\phi_x(\theta) = -\frac{\partial}{\partial\theta} F(x \mid \theta) \tag{1}
$$

For this to be meaningful certain restrictions will have to be imposed on $F$. Apart from the existence of the derivative these are that $F(x \mid \theta)$ decreases with $\theta$ for each $x$,

$$
\lim_{\theta \to U} F(x \mid \theta) = 0 \quad \text{and} \quad \lim_{\theta \to L} F(x \mid \theta) = 1,
$$

where $L$ and $U$ are the upper and lower limits for $\theta$, which may be infinite. It will be supposed that $F$ satisfies these conditions, which we term (B).

We now consider under what conditions (1) can be derived from a Bayesian argument starting from a prior distribution $p(\theta)$. Suppose the probability density

$$
f(x \mid \theta) = \frac{\partial}{\partial x} F(x \mid \theta)
$$

exists, and let

$$
\rho(x) = \int f(x \mid \theta)\, p(\theta)\, d\theta, \tag{2}
$$

then the posterior distribution for $\theta$ is

$$
f(x \mid \theta)\, p(\theta) / \rho(x). \tag{3}
$$

We enquire under what conditions on $F$ can a prior distribution be found such that (1) and (3) are identical. We require to solve the equation

$$
-\frac{\partial}{\partial\theta} F(x \mid \theta) = \frac{\partial}{\partial x} F(x \mid \theta)\, p(\theta) / \rho(x) \tag{4}
$$

for $p(\theta)$, with $\rho(x)$ given by (2). The equation may be rewritten

$$
\frac{-\dfrac{\partial}{\partial\theta} F(x \mid \theta)}{\dfrac{\partial}{\partial x} F(x \mid \theta)} = \frac{p(\theta)}{\rho(x)}. \tag{5}
$$

Hence a necessary condition for (4) to have a solution is that the left-hand side of (5) be equal to the product of a function of $\theta$ and a function of $x$. If this is not possible then no solution of (4) is possible. Let us therefore consider when $F$ satisfies (5) for given $p$ and $\rho$.

We have

$$
\frac{1}{\rho(x)} \frac{\partial F}{\partial x} + \frac{1}{p(\theta)} \frac{\partial F}{\partial \theta} = 0.
$$

The subsidiary equations are

$$
\frac{dF}{0} = \rho(x)\, dx = p(\theta)\, d\theta.
$$

So either $F$ is a constant, which is impossible, or

$$
F = G\big(R(x) - P(\theta)\big)
$$

where $G$ is an arbitrary function and $R(x)$, $P(\theta)$ are the integrals of $\rho(x)$ and $p(\theta)$ respectively.

Consider a new random variable and a new parameter given by

$$
u = R(x) \qquad \tau = P(\theta)
$$

Then

$$
F = G(u - \tau), \tag{6}
$$

when it is easy to see that (B) obtains. It is easy to verify that when $F$ satisfies (6) there is a solution to (4) with $p(\tau) = \text{constant}$, so the condition, which is recognizable as (A), is both sufficient and necessary.

Note that if $\theta$ is a scale parameter for $x$, (A) obtains with $u = \log x$ and $\tau = \log \theta$. The prior distribution for $\theta$ is $p(\theta) = \theta^{-1}$.

---

## 2. A Consistency Condition for Fiducial Probability

In the above analysis let $x$ be a sufficient statistic for $\theta$. Let $y$ be a sufficient statistic for $\theta$ obtained from a second sample, independent of the first. Let $\psi_{x,y}(\theta)$ be the Bayes' posterior distribution for $\theta$ using Bayes' theorem for $y$ with prior distribution $\phi_x(\theta)$. Let $\phi_{x,y}(\theta)$ be the fiducial distribution for $\theta$ given $x$ and $y$, when initially nothing is known about $\theta$. In this section we investigate whether $\phi_{x,y}(\theta) = \psi_{x,y}(\theta)$. According to Fisher (1956) this is true, though no formal proof is given, for he says (p. 51) "The concept of probability involved [in the fiducial argument] is entirely identical with the classical probability of the early writers, such as Bayes". Again he says (p. 125) "This fiducial distribution supplies information of exactly the same sort as would a distribution given *a priori*". $\phi_{x,y}(\theta)$ is the result of applying the fiducial argument to the combined sample of $x$ and $y$; $\psi_{x,y}(\theta)$ is the result of applying the fiducial argument to $x$ and then using the resulting fiducial distribution $\psi_y(\theta)$ as a prior distribution for $y$. Clearly it is necessary for the two to be the same in order that the fiducial argument be consistent: or else Fisher's statements are false.

We shall only investigate the relationship in the case where $x$ and $y$ have the same distribution and this distribution admits of a single sufficient statistic for any size of sample. It is well known that the distribution must then have probability density of the form

$$
f(x)\, g(\theta)\, e^{u(x) h(\theta)},
$$

and, without loss of generality, we may write it

$$
f(x)\, g(\theta)\, e^{x\theta} = f(x \mid \theta). \tag{7}
$$

The joint distribution for $x$ and $y$ is clearly seen to admit $x + y$ as a sufficient statistic and hence $\phi_{x,y}(\theta)$ is a function of $x + y$ and $\theta$, only. We have first then to consider the conditions under which $\psi_{x,y}(\theta)$ is a function of $x + y$ and $\theta$, only. Again it will turn out that these necessary conditions will be sufficient.

Clearly

$$
\psi_{x,y}(\theta) = \frac{\phi_x(\theta)\, f(y \mid \theta)}{\displaystyle\int \phi_x(\theta)\, f(y \mid \theta)\, d\theta}.
$$

If this is to be a function of $x + y$ and $\theta$ the same must be true of

$$
\phi_x(\theta)\, f(y \mid \theta) \,/\, \phi_x(\theta')\, f(y \mid \theta'),
$$

where $\theta' \neq \theta$, and therefore this expression must be unaltered by an interchange of $x$ and $y$. Inserting the expressions for $f(y \mid \theta)$ and $f(y \mid \theta')$ given by (7) we have, on rearrangement, that

$$
\frac{\phi_x(\theta)\, \phi_y(\theta')}{\phi_x(\theta')\, \phi_y(\theta)} = e^{(x-y)(\theta-\theta')}.
$$

By taking logarithms and regarding $\theta'$ and $y$ as constants so that the only variables are $\theta$ and $x$ it is easy to see that

$$
\phi_x(\theta) = F(x)\, G(\theta)\, e^{x\theta}
$$

where $F$ and $G$ are undetermined functions. It thus follows that $\phi_{x,y}(\theta)$ can only be equal to $\psi_{x,y}(\theta)$ if

$$
\frac{\phi_x(\theta)}{f(x \mid \theta)} = \frac{-\dfrac{\partial}{\partial\theta} F(x \mid \theta)}{\dfrac{\partial}{\partial x} F(x \mid \theta)} = \frac{F(x)\, G(\theta)}{f(x)\, g(\theta)},
$$

that is, a product of a function of $x$ and a function of $\theta$. This is (5) above and therefore the equality is only possible under conditions (A). Under these conditions we saw that the fiducial distribution $\phi_x(\theta)$ is a Bayes' distribution with prior distribution given by a uniform distribution for the location parameter. It is then clear that $\phi_{x,y}(\theta) = \psi_{x,y}(\theta)$ because Bayes' theorem is consistent. Hence (A) is both a necessary and sufficient condition. Thus the fiducial argument is consistent if, and only if, it is equivalent to a Bayesian argument under conditions (A) with uniform prior distribution for the location parameter.

Some readers may like to have an example of a breakdown in the relation, since Fisher's examples deal with location parameters. Perhaps the simplest is to take

$$
f(x \mid \theta) = \frac{\theta^2}{\theta + 1}\,(x + 1)\, e^{-x\theta}, \qquad x > 0,\ \theta \geq 0.
$$

The distribution for $z = x + y$ is

$$
f_2(z \mid \theta) = \left(\frac{\theta^2}{\theta + 1}\right)^{2} \left(z + z^2 + \frac{1}{6} z^3\right) e^{-z\theta}
$$

and

$$
\phi_{x,y}(\theta) = \frac{e^{-z\theta}}{(\theta + 1)^3} \left[\theta^3 \left(2z^2 + \frac{4}{3} z^3 + \frac{1}{6} z^4\right) + \theta^4 \left(z^2 + z^3 + \frac{1}{6} z^4\right)\right]. \tag{8}
$$

On the other hand\*

$$
\phi_x(\theta) = \frac{\theta}{(\theta + 1)^2} \big[1 + (1 + \theta)(1 + x)\big]\, x e^{-x\theta}
$$

so that $\psi_{x,y}(\theta)$ is proportional to

$$
\frac{\theta^3}{(\theta + 1)^3}\, x(1 + y) \big[1 + (1 + \theta)(1 + x)\big]\, e^{-(x+y)\theta}, \tag{9}
$$

the constant of proportionality involving $x$ and $y$ only. (8) and (9) are clearly distinct. Notice that (9) is not even equal to $\psi_{y,x}(\theta)$; that is, the order of the two samples would affect the analysis. In fact the idea that $\psi$ must be invariant under a change of order was basic to our proof. It is obvious that (B) is satisfied both for the $x$ and $z$ distributions.

> \* Because of the way $\theta$ appears in this example the usual minus sign in front of $\partial F / \partial \theta$ is not needed in the definition of the fiducial distribution. Here $F(x \mid \theta)$ *increases* with $\theta$ for each $x$.

---

## 3. The Exponential Family and Condition (A)

The analysis in §2 naturally raises the question: under what circumstances does a member of the exponential family, given by (7), satisfy condition (A)? We notice that the family is as wide as the class of functions $f(x)$ which are such that

$$
\int_{-\infty}^{\infty} f(x)\, e^{x\theta}\, dx \tag{10}
$$

converges for all $\theta$ in its permissible range. (10) is then equal to $g(\theta)^{-1}$. Consequently it is sufficient to determine the conditions on $f(x)$. Notice that we are supposing the distributions written in the standard form (7).

We require to find the circumstances under which there exists $u(x)$ and $\tau(\theta)$ such that

$$
\frac{dx}{du}\, f\big(x(u)\big)\, g\big(\theta(\tau)\big)\, e^{\theta(\tau)\, x(u)} \tag{11}
$$

is a function of $u + \tau$ only. (It avoids minus signs to write $u + \tau$ instead of $u - \tau$.) Taking logarithms of (11) and differentiating first with respect to $u$ and then with respect to $\tau$ we easily obtain that

$$
\theta'(\tau)\, x'(u)
$$

is a function of $u + \tau$ only; from which it follows that

$$
\begin{aligned}
\log\big(\theta'(\tau)\big) &= a\tau + b_1, \\
\log\big(x'(u)\big) &= au + b_2,
\end{aligned}
$$

where $a$, $b_1$ and $b_2$ are constants. There are two cases to consider, either $a \neq 0$ or $a = 0$. In the former

$$
\theta(\tau) = c_1 e^{a\tau} + b_1, \qquad x(u) = c_2 e^{au} + b_2 \tag{12}
$$

and in the latter

$$
\theta(\tau) = c_1 \tau + b_1, \qquad x(u) = c_2 u + b_2. \tag{13}
$$

In either case there is no loss in generality in supposing $b_1 = b_2 = 0$ since the terms they give rise to in (11) can always be incorporated in the functions of $u$ and $\tau$ before the exponential. Consider (12) first. We can write (11) as

$$
f_1(u)\, g_1(\tau)\, \exp\{c\, e^{a(u+\tau)}\} \tag{14}
$$

from which it follows that $f_1(u)\, g_1(\tau)$ is a function of $u + \tau$ only and as before

$$
f_1(u) = \beta_1 e^{\alpha u}, \qquad g_1(\tau) = \beta_2 e^{\alpha\tau}.
$$

Hence (14) becomes

$$
\gamma\, e^{\alpha(u+\tau)} \exp\{c\, e^{a(u+\tau)}\}, \tag{15}
$$

or in terms of $x$ and $\theta$ (from (12)), and allowing for the differential element,

$$
C x^{-1} (\theta x)^{\alpha/a}\, e^{\theta x} \tag{16}
$$

where $C$ is chosen to make the integral one. To do this we note that in the usual derivation of the exponential family we are considering the case where the range of $x$ does not depend on the parameter: hence in (15) $u$ and $\tau$ must range over the whole real line and hence in (16) $x$ and $\theta$ go over the half-lines terminating at the origin. We can suppose $x \geq 0$ when, for convergence we must have $\theta < 0$ and $\alpha/a > 0$. Under these conditions $C$ can be found and writing $\alpha/a = \lambda$ and changing the sign of $\theta$ we have

$$
\theta^{\lambda} x^{\lambda - 1} e^{-x\theta} / \Gamma(\lambda) \tag{17}
$$

the well-known $\Gamma$-distribution.

Analogous considerations starting from (13) instead of (12) give in place of (17)

$$
\frac{1}{\sqrt{2\pi}\,\sigma} \exp\left\{-\frac{1}{2} \frac{(x - \theta)^2}{\sigma^2}\right\}, \tag{18}
$$

the normal distribution.

(17) and (18) are the only distributions of the exponential family which satisfy (A) and therefore the fiducial argument is only consistent in the case of a single sufficient statistic when the distribution is of the gamma or normal forms, or transformable thereto.

### Note added in proof

Since this paper was written it has been discovered that the question discussed in this section has also been considered by E. B. Dynkin, "Necessary and sufficient statistics for a family of probability distributions". *Uspehi Matem. Nauk.* (N.S.) **6**, no. 1(41), 68–90, (1951), in Russian. A translation of this paper is available from the Statistical Laboratory, University of Cambridge.

---

## References

- Fisher, R. A. (1956), *Statistical Methods and Scientific Inference*. Edinburgh: Oliver & Boyd.
- Grundy, P. M. (1956), "Fiducial distributions and prior distributions: An example in which the former cannot be associated with the latter", *J. R. Statist. Soc.*, B, **18**, 217–221.
