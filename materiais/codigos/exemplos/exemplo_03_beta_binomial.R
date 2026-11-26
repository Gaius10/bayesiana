##################################
# Exemplo 2.2 Figura 2.1
# Modelo Binomial com priori Beta.
#
# Ricardo Ehlers
##################################

 library(xtable)
binom.beta <- function(n,y,a,b){
 a1 = a+y
 b1 = b+n-y
 curve(dbeta(x,y+1,n-y+1),0,1,xlab=expression(theta),ylab='',bty='n')
 curve(dbeta(x,a,b),0,1,add=T,lty=2)
 curve(dbeta(x,a1,b1),0,1,add=T,lty=3)
 legend(.00001,3.5,c('veross','priori','posteriori'),lty=1:3,bty='n')
}

par(mfrow=c(1,3))
binom.beta(12,9,1,1)
binom.beta(12,9,2,2)
binom.beta(12,9,1,3)


# Probabilidades preditivas (marginais)

beta.binomial <- function(n,a,b) {
 m = matrix(0,n+1,2)
 m[,1]= 0:n
 for (x in 0:n) 
 m[(x+1),2]=round(choose(n,x)*beta(a+x,b+n-x)/beta(a,b),4)
 m
}

n=12
k=0:n
pred=cbind(k,beta.binomial(n,1,1)[,2],beta.binomial(n,2,2)[,2],beta.binomial(n,1,3)[,2])
colnames(pred)=c('k','Beta(1,1)','Beta(2,2)','Beta(1,3)')

pred

xpred=xtable(pred,display=c('s','d','f','f','f'),digits=4,align=c('c','l','c','c','c'))
print(xpred,include.rownames=F)

# A influencia da priori diminui conforme aumenta a amostra.

par(mfrow=c(2,2))
binom.beta(12 ,     9 ,1,4);title('n=12')
binom.beta(20 ,.75*20 ,1,4);title('n=20')
binom.beta(50 ,.75*50 ,1,4);title('n=50')
binom.beta(100,.75*100,1,4);title('n=100')
