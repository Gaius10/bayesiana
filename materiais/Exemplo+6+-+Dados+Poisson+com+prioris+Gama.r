########################################
# Secao 2.3.2
# Dados Poisson com priori Gama.
#
# Ricardo Ehlers
########################################

poisson.gama <- function(x,a0,b0){
#x,a,b: vetor de observacoes e parametros da priori
 n=length(x)
 tx = sum(x)
 a1= a0+tx
 b1= b0+n
 curve(dgamma(x,a1,b1 ),0,10,col=3,xlab=expression(theta),ylab='') 
 curve(dgamma(x,tx+1,n),0,10,col=1,add=T)
 curve(dgamma(x,a0,b0 ),0,10,col=2,add=T)
 legend(.5,0.5,c('veross','priori','posteriori'),fill=1:3,bty='n')
}


x=c(5,13,6,3,3,4,6,7,7,6)

par(mfrow=c(2,2))
poisson.gama(x,.1,.1)
poisson.gama(x,3,2)
poisson.gama(x,1,3)
poisson.gama(x,4,1)



poisson.gama <- function(x,a0,b0){
#x,a,b: vetor de observacoes e parametros da priori
 n=length(x)
 tx = sum(x)
 a1= a0+tx
 b1= b0+n
 curve(dgamma(x,a1,b1 ),0,10,col=3,xlab=expression(theta),ylab='')
 curve(dgamma(x,tx+1,n),0,10,col=1,add=T)
 curve(dgamma(x,a0,b0 ),0,10,col=2,add=T)
 legend(4.5,1.5,c('veross','priori','posteriori'),fill=1:3,bty='n')
}


#Dados simulados
n=50; y=rpois(n,4)

par(mfrow=c(2,2))
poisson.gama(y,.1,.1)
poisson.gama(y,3,2)
poisson.gama(y,1,3)
poisson.gama(y,4,1)
