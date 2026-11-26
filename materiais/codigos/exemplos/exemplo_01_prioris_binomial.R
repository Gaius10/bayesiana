#########################################
# Exemplo 1.2 Figuras 1.1, 1.2 e 1.3
# Modelo binomial com prioris diferentes.
#
# Ricardo Ehlers
#########################################


# Priori normal truncada.

dnorm.t <- function(x,mean=0,sd=1){
 aux=pnorm(1,mean,sd)-pnorm(0,mean,sd)
 dnorm(x,mean,sd)/aux
}

curve(dnorm.t(x,0.5,0.5),0,1,lty=1,ylim=c(0,3),ylab=expression(pi(theta)),xlab=expression(theta))
curve(dnorm.t(x,0.0,0.5),0,1,lty=6,add=T)
curve(dnorm.t(x,1.0,0.5),0,1,lty=3,add=T,lwd=2)
curve(dnorm.t(x,2.0,0.5),0,1,lty=5,add=T,lwd=2)
legend(0.3,2.6,c('N(0.5,0.5)','N(0,0.5)','N(1,0.5)','N(2,0.5)'),lty=c(1,6,3,5),lwd=c(1,1,2,2))


# Priori normal-logistica.

dlogist <- function(x,mean,sd){
 z = log(x/(1 - x))
 dnorm(z,mean,sd)/(x-x**2)
}

curve(dlogist(x,-1,.5),0,1,ylab=expression(pi(theta)),xlab=expression(theta))
curve(dlogist(x,1,1),0,1,lty=4,add=T)
curve(dlogist(x,0,4),0,1,lty=5,add=T,lwd=2)
legend(0.45,4,c('N(-1,0.5)','N(1,1)','N(0,4)'),lty=c(1,3,5),lwd=c(1,1,2))

# Priori Beta.

curve(dbeta(x,1.5,4),0,1,lty=1,ylim=c(0,5),ylab=expression(pi(theta)),xlab=expression(theta))
curve(dbeta(x,2,0.7),0,1,lty=6,add=T)
curve(dbeta(x,7,1.5),0,1,lty=3,add=T,lwd=2)
curve(dbeta(x,3,3)  ,0,1,lty=5,add=T,lwd=2)
legend(.2,4,c('Beta(1.5,4)','Beta(2,0.5)','Beta(7,1.5)','Beta(3,3)'),lty=c(1,6,3,5),lwd=c(1,1,2,2))

#curve(dbeta(x,1,1)  ,0,1,lty=8,add=T,lwd=2)