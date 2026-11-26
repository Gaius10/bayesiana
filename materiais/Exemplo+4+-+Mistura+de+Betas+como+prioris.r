###########################################
# Exemplo 2.6 Figura 2.2
# Misturas de densidades Beta como prioris.
#
# Ricardo Ehlers
###########################################

curve(0.33*dbeta(x,4,10)+0.33*dbeta(x,15,28)+0.33*dbeta(x,50,70),0,1,xlab=expression(theta),ylab='',bty='n')
curve(0.25*dbeta(x,3,8) + 0.75*dbeta(x,8,3),0,1,add=T,lty=2)
legend(0.45,3.5,c('.33B(4,10)+.33B(15,28)+.33B(50,70)','.25 B(3,8)+.75 B(8,3)'),lty=1:2,bty='n')

