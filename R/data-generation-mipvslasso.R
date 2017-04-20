##### Função geradora de dados #####
# T = número de observações
# N = Número de variáveis candidatas
# q = Número de variáveis relevantes (as primeiras q variáveis)
# iid = Se TRUE, dados são IID no tempo e no crossection, se FALSE sigma nao é identidade
# e um componente autorregressivo é adicionado como primeira variável relevante. 

# As médias das variáveis candidatas são -1 para a primeira 1, para a segunda, -1 para a terceira...
# Os betas das variáveis relevantes são -0.3 para a primeira, 0.3 para a segunda, -0.3 para a terceira...
# Se houver coeficiente autorregressivo, ele será a primeira variável e terá valor 0.6
require(clusterGeneration)

gensig=function(N,iid=FALSE){
  # require(clusterGeneration)
  if(iid==TRUE){
    sigma=diag(N)
  }else{
    sigma=genPositiveDefMat(N,rangeVar=c(0.1,1))$Sigma
  }
  return(sigma)
}


dgp <- function(T,N,q,iid=FALSE,sigma = NULL){
  
  if (is.null(sigma)) {
      sigma = gensig(N, iid)
  }
  require(clusterGeneration)
  require(mvtnorm)
  
  mu=rep(0,N)
  for(i in 1:N){
    mu[i]=(-0.5/0.5)^i
  }
  
  data=rmvnorm(T+100,mean=mu,sigma=sigma)
  rho=0.6
  beta=mu[1:q]*0.3
  
  if(iid==TRUE){
    y=data[,1:q]%*%beta+rnorm(T+100,0,1)
    data=data[-1,]
  }else{
    y=rep(0,T+100)
    for(i in 2:(T+100)){
      y[i]=rho*y[i-1]+data[i,2:q]%*%beta[2:q]+rnorm(1,0,1)
    }
    data=cbind(y[-length(y)],data[-1,2:ncol(data)])
  }
  
  y=y[-1]
  
  return(list(Y=y[(100):length(y)],X=data[100:nrow(data),],beta = c(beta, rep(0,N-q))))
  
}

