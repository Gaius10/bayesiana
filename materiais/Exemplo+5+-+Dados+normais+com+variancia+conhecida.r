########################################
# Secao 2.3.1
# Dados normais com variancia conhecida.
#
# Ricardo Ehlers
########################################

norm.norm <- function(x,sigma,mu0,tau0,plot=F){
 n = length(x)
 xbar = mean(x)
 ep = sigma/sqrt(n)
 sigma2 = sigma**2
 tau1 = n*(1/sigma2)+(1/tau0)
 mu1 = (n*(1/sigma2)*xbar+(1/tau0)*mu0)/tau1
 if (plot) {
   curve(dnorm(x,xbar,ep),xbar-3*ep,xbar+3*ep)
   curve(dnorm(x,mu0,sqrt(tau0)),add=T,col=2)
   curve(dnorm(x,mu1,1/sqrt(tau1)),add=T,col=3)
   legend(-0.8,1.2,legend=c('veross.','priori','posteriori'),col=1:3,lty=c(1,1,1))
 }
 return(list(mu1=mu1, tau1=tau1)) 
}

norm.norm(x=rnorm(50),sigma=2,mu0=0,tau0=1,plot=T)
