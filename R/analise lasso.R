#### Este código lê as tabelas fornecidas pelas otimizações do LASSO e MIP (para
#### seleção das variáveis da regressão quantílica) e exibe um gráfico de tudo


# Functions ---------------------------------------------------------------

loss.function <- function(y,beta, alpha, P = 12) {
  n <- length(y)
  loss <- 0
  for (i in 1:(n-P)) {
    y.hat = sum(rev(y[i:(i+P-1)]) * beta[2:(P+1)]) + beta[1]
    epsilon = y[i+P] - y.hat
    loss = loss + (0.5 - (0.5 - alpha)*(sign(epsilon))) * abs(epsilon)
  }
  return(loss)
}


SIC <- function(y,beta, alpha, P = 12) {
  n <- length(y)
  nvar <- sum(abs(beta) > 0.00000001)
  return(n * log(loss.function(y, beta, alpha, P)/n) + 1/2 * nvar * log(n))
}

# SIC(y,beta,alpha, P=12)


lagmatrix <- function(x, lags) {
  n = length(x)
  I = (lags+1):(n)
  x_estim = matrix(0,n-lags, lags)
  for (i in 1:lags) {
    x_estim[ ,i] = x[I - i]
  }
  return(x_estim)
}



library(stringr); library(ggfortify); library(rjulia)
tabela <- read.csv("RegressãoQuantílica_STREET/table-betas-sellassonorm-alpha-095.csv", header = FALSE)
destino <- "Documento Regressao Quantilica/Figuras/"
icaraizinho <- read.csv("RegressãoQuantílica_STREET/icaraizinho.csv", header = FALSE)[,1]
# View(abs(tabela) > 0.0000001)
# y <- icaraizinho - mean(icaraizinho)
# 
# qplot(-log(lambdas), nvar, type = "l")
# qplot(-log(lambdas), nvar)
# 
# 
# tabela.coef <- t(tabela[-1, ])
# colnames(tabela.coef) <- str_c("beta_", 1:12, "")
# rownames(tabela.coef) <- lambdas
# autoplot(ts(tabela.coef), facets = FALSE)
# 


# Gráficos e BIC ----------------------------------------------------------

library(stringr); library(xtable); library(quantreg); library(dplyr)


# Finds 

pasta.origem <- "RegressãoQuantílica_STREET/"
# pasta.destino <- "Documento Regressao Quantilica/Figuras/selecao-lasso/"
nomes <- list.files(path = pasta.origem, pattern = "table-betas-postlasso-")
# nomes <- list.files(path = pasta.origem, pattern = "table-betas-sellassonorm-")
P = 12
epsilon = 0.000000000000001
maxK = P
arquivos.tabelas.lasso <- str_c(pasta.origem,nomes)

nomes <- list.files(path = pasta.origem, pattern = "table-betas-selecaointeira-")
arquivos.tabelas.mip <- str_c(pasta.origem,nomes)

# lambdas = exp(seq(-5,5,0.1))  # Lasso normalizando
# lambdas = exp(seq(-2,10,0.3))  # Lasso sem normalizar
lambdas = exp(seq(-2,9,0.1))


matriz_com_lags =  lagmatrix(icaraizinho, 12)
matriz_correlacao = cor(matriz_com_lags)



i=1
for (i in 1:length(arquivos.tabelas.lasso)) {
  tabela <- read.table(arquivos.tabelas.lasso[i], header = FALSE, sep = ",")
  nvar <- colSums(abs(tabela) > 0.000000000001)
  rownames(tabela) <- str_c("beta_", 0:12, "")
  colnames(tabela) <- lambdas
  
  # extracts alpha from name of table
  tmp <- str_extract(nomes[i], pattern = "\\d+")
  tmp.strlen <- str_length(tmp)
  tmp.parts <- str_sub(tmp, c(1,2), c(1,tmp.strlen))
  alpha <- as.numeric(str_c(tmp.parts[1], ".", tmp.parts[2]))
  
  # fills sic for each value of lambda
  keep.sic.lasso <- rep(0, ncol(tabela))
  keep.loss.fn <- rep(0, ncol(tabela))
  for (coluna in 1:ncol(tabela)) {
    beta <- tabela[, coluna]
    keep.sic.lasso[coluna] <- SIC(icaraizinho, beta, alpha, P)
    keep.loss.fn[coluna] <- loss.function(icaraizinho, beta, alpha, P)
  }
  
  # plot(keep.sic.lasso, main = str_c("sic ", alpha))
  # plot(log(keep.loss.fn), main = str_c("loss fn ", alpha))
  # plot(nvar)
  # 
  # min.sic <- which.min(keep.sic.lasso)
  # nvar[min.sic]
  
  # For each number of variable, finds the one with the best SIC for the LASSO
  # keep.best.K.i <- rep(0, maxK)
  # for (k in 1:maxK) {
  #   i.k <- which(nvar == k)
  #   if (length(i.k)==0) { # if there is no estimation with k values, puts NA
  #     keep.best.K.i[k] <- NA
  #   } else {
  #     min.sic <-  which.min(keep.sic.lasso[i.k])
  #     keep.best.K.i[k] <- i.k[min.sic]  
  #   }
  # }
  

  # For each number of variable, finds the one with the best SIC for the post-LASSO
  keep.best.K.i <- rep(0, maxK)
  for (k in 1:maxK) {
    i.k <- which(nvar == (k+1))
    if (length(i.k)==0) { # if there is no estimation with k values, puts NA
      keep.best.K.i[k] <- NA
    } else {
      min.sic <-  which.min(keep.sic.lasso[i.k])
      keep.best.K.i[k] <- i.k[min.sic]  
    }
  }
  
  
  tabela.mip <- read.csv( str_c(pasta.origem, "table-betas-selecaointeira-alpha-",
                                str_replace(alpha, "\\.", ""), ".csv"), header = FALSE)
  
  # which mip table contains the desired alpha 
  i.mip <-  str_detect(arquivos.tabelas.mip, str_c("-",str_replace(alpha, "\\.", ""), "\\."))
  tabela.mip <- read.table(arquivos.tabelas.mip[i.mip], header = FALSE, sep = ",")
  keep.sic.mip <- rep(0,maxK)
  # For each number of variables, finds the SIC
  for (k in 1:maxK) {
    beta <- tabela.mip[, k]
    keep.sic.mip[k] <- SIC(icaraizinho, beta, alpha, P)
  }
  
  distance <- rep(0, maxK)
  distance2 <- rep(0,maxK)
  k = 3
  for(k in 1:maxK) {
    index.sol <- keep.best.K.i[k]
    if (is.na(index.sol)) {
      distance[k] <- NA    
    } else {
      beta.lasso <- tabela[keep.best.K.i[k]]
      beta.mip <- tabela.mip[k]
      sel_lasso <- (abs(beta.lasso) > epsilon)[-1] # variables that keep which covariates are selecret
      sel_mip <- (abs(beta.mip) > epsilon)[-1]
      distance[k] <- 1/(2*k) * sum(abs(sel_lasso - sel_mip))
      
      L_mip = which(sel_mip)
      L_lasso = which(sel_lasso)
      #### Chama código em Júlia para fazer a otimização
      julia_init()
      r2j(matriz_correlacao, "rho")
      r2j(L_mip, "L_mip")
      r2j(L_lasso, "L_lasso")
      # r2j(arquivo_procedimento, "arquivo_procedimento")
      # Código júlia
      julia_void_eval("print(convert(Array{Int64,1},L_lasso))")
      julia_void_eval("@show size(L_lasso)[1]")
      
      julia_void_eval("using JuMP, Gurobi;\\
                      @show  pwd();\\
                      include(\"R/procedure_distancia_R.jl\");\\
                      @show rho;\\
                      L_lasso = convert(Array{Int64,1},L_lasso); L_mip = convert(Array{Int64,1},L_mip);\\
                      distancia, matriz_delta = procedure_distancia(L_lasso, L_mip, rho);") #                
      saidas <- j2r("(distancia, matriz_delta)")              
      distance2[k] <- saidas[[1]]
      
    }
  }
  
  
########### Begin plot ################
  alpha_out <- str_replace(alpha, "\\.", "")
  pdf(str_c(destino,"SIC",alpha_out, ".pdf"), width = 6, height = 4.5)
  ## add extra space to right margin of plot within frame
  par(mar=c(5, 4, 4, 6) + 0.1)
  limits=c(min(keep.sic.mip),max(keep.sic.lasso[keep.best.K.i], na.rm = TRUE))


  plot(distance2, xlab="", ylab="", ylim=c(0,1),
       axes=FALSE, type="h", col=rgb(0, 1, 0,0.5), lwd=10)
  axis(4, ylim=c(0,1), las=1)

  par(new=TRUE)

  ## Plot first set of data and draw its axis
  ## Mark the best points of both MIP and LASSO
  i.best.lasso <- which.min(keep.sic.lasso[keep.best.K.i])
  i.best.mip <- which.min(keep.sic.mip)
  
  plot(keep.sic.mip , axes=FALSE, ylim=limits, ylab = "",
       xlab = "", main = str_c("alpha = ", alpha), col = 4)
  # points(keep.sic.mip)
  points(i.best.mip, keep.sic.mip[i.best.mip], col = 2)
  points(keep.sic.lasso[keep.best.K.i], pch=3, col=4)
  points(i.best.lasso, keep.sic.lasso[keep.best.K.i][i.best.lasso], pch = 3, col = 2)
  
  
  legend(x="topright", c("MIP","Post-Lasso","Distance"), pch=c(1,3,15), col=c(4,4,rgb(0, 1, 0,0.5)))
  axis(2, ylim=limits,col="black",las=1)  ## las=1 makes horizontal labels
  mtext("Schwarz Information Criteria",side=2,line=2.5)
  box()

  
  ## Allow a second plot on the same graph

  ## Plot the second plot and put axis scale on right
  ## a little farther out (line=4) to make room for labels
  mtext("Distance",side=4,col=1,line=4)
  # axis(4, ylim=c(0,1), col="red",col.axis="red",las=1)

  yLabels <- seq(0.2, 0.8, 0.2)

  ## Draw the time axis
  axis(1,1:maxK)
  mtext("K (Number of covariates included)",side=1,col="black",line=2.5)
  dev.off()
  ########### end plot ################

  
}








# LIXO --------------------------------------------------------------------

# 
# 
# 
# 
# f=function(x) (x-2)^3
# f = function(x) 10*x
# f1 = function(x) {
#   if (x > 0) return(x^2)
#   return(-3*x^2)
# }; f = function(x) sapply(x,f1)
# tmp <- rnorm(1000,mean = 1)
# plot(density(tmp))
# plot(density(f(tmp)))
# 
# # x <- seq(-4,6,0.01)
# plot(tmp,f(tmp))
# 
# 
# 
# x.bar <- rep(0,100000)
# for(i in 1:100000) {
#   x.bar[i] = mean(rgamma(100000,2))
# }
# plot(density(x.bar))
# 
# 
# 
# 
# 
# 
# 
# 
# tabela
# 
# library(rjulia)
# r2j(tabela, "table")
# julia_eval("tmp = pwd(); x = 1")
# tmp = j2r("tmp")
# tmp
# 

