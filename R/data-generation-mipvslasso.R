##### Função geradora de dados #####
# T = número de observações
# N = Número de variáveis candidatas
# q = Número de variáveis relevantes (as primeiras q variáveis)
# iid = Se TRUE, dados são IID no tempo e no crossection, se FALSE sigma nao é identidade
# e um componente autorregressivo é adicionado como primeira variável relevante. 

# As médias das variáveis candidatas são -1 para a primeira 1, para a segunda, -1 para a terceira...
# Os betas das variáveis relevantes são -0.3 para a primeira, 0.3 para a segunda, -0.3 para a terceira...
# Se houver coeficiente autorregressivo, ele será a primeira variável e terá valor 0.6


gensig=function(N,iid=FALSE){
   require(clusterGeneration)
  if(iid==TRUE){
    sigma=diag(N)
  }else{
    sigma = genPositiveDefMat(N,covMethod=c("unifcorrmat"), rangeVar = c(1,1))$Sigma
  }
  return(sigma)
}


dgp <- function(T, Nin, q, iid=FALSE, autoregressive = TRUE, sigma = NULL, rho=0, sigma_e = 1){


  if (is.null(sigma)) {
      sigma = gensig(N, iid)
  }
  
  mu = rep(0,N)
  for(i in 1:N) {
    mu[i] = (-0.5/0.5)^i
  }
  
  data=rmvnorm(T+100,mean=mu,sigma=sigma)
  beta=mu[1:q]*0.3
  
  if(autoregressive==FALSE){
    y=data[,1:q]%*%beta + rnorm(T+100,0,sqrt(sigma_e))
    data=data[-1,]
    rho = 0
  }else{
    
    y=rep(0,T+100)
    for(i in 2:(T+100)){
      y[i]=rho*y[i-1]+data[i,1:q]%*%beta[1:q]+rnorm(1,0,sqrt(sigma_e))
    }
    data=cbind(y[-length(y)],data[-1,])
    # beta = c(beta[1:q], rep(0,N-q))
  }
  
  y=y[-1]
  
  return(list( Y = y[(100):length(y)],
               X=data[100:nrow(data),],
               rho = rho,     beta = c(beta, rep(0,N-q)),
               sigma = sigma))
  
}



require(clusterGeneration)
require(mvtnorm)
require(normtest)


# T = 10000; N = 3; q = 3; iid = FALSE; rho = 0.0; autoregressive = FALSE; sigma_e = 0.3
# sigma = genPositiveDefMat(N,covMethod=c("unifcorrmat"), rangeVar = c(1,1))$Sigma
# sigma = gensig(N, iid)
# data = rmvnorm(T+100,mean=mu,sigma=sigma)
# dados = dgp(T,N,q, iid, rho = rho, autoregressive =  autoregressive,sigma_e = sigma_e)
# X = dados$X
# y = dados$Y
# sigma = dados$sigma
# apply(X,2,sd)
# apply(X,2,mean)
# beta_true = dados$beta
# if (autoregressive) {
#   print(mean(y))
#   # jb.norm.test(y)
#   print(E_yt <- (q)*0.3/(1-rho))
#   print(var(y))
#   print(V_yt <- ( (t(beta_true) %*% sigma %*% beta_true)[1,1] + sigma_e) / (1-rho^2))
# } else {
#   print(mean(y))
#   # jb.norm.test(y)
#   print(E_yt <- (q)*0.3)
#   print(var(y))
#   print(V_yt <- ( (t(beta_true) %*% sigma %*% beta_true)[1,1] + sigma_e))
# }
# 
# 
# plot(density(y))
# lines(valores_x <- seq(E_yt - 6*sqrt(V_yt),  E_yt + 6*sqrt(V_yt), 0.001), dnorm(valores_x, E_yt, sqrt(V_yt)), col = 4  )
# abline(v = E_yt, col = 2)
# 
# 
