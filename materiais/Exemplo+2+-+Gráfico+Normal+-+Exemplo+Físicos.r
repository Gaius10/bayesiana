######################################################
# Exemplo 1.3 Figura 1.4
# Dado normal com variancia conhecida e priori normal.
#
# Ricardo Ehlers
######################################################

norm.norm <- function(x,mu0,tau0,s0){
 precisao = 1/tau0 + length(x)/s0
 tau1 = 1/precisao
 w = (1/tau0)/precisao
 mu1 = w*mu0 + (1-w)*mean(x)
 return(list(m=mu1,tau=tau1))
}

curve(dnorm(x,900,20),700,1000,ylab='',ylim=c(0,.022),lty=1,xlab=expression(theta))
curve(dnorm(x,890,17.9),add=T,lty=2)
curve(dnorm(x,800,40),add=T,lty=1)
curve(dnorm(x,840,35.7),add=T,lty=2)
curve(dnorm(x,850,40),add=T,lty=3)
legend(710,.018,lty=c(1,2,3),c("priori","posteriori","verossimilhanca"),bty='n')
text(950,0.015,'Fisico A')
text(720,0.0055,'Fisico B')

