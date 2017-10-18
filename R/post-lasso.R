## This script takes all tables from Lasso estimations and use them as input for
## the post-lasso estimation

lagmatrix <- function(x, lags) {
  n = length(x)
  I = (lags+1):(n)
  x_estim = matrix(0,n-lags, lags)
  for (i in 1:lags) {
    x_estim[ ,i] = x[I - i]
  }
  return(x_estim)
}


library(stringr); library(quantreg)

icaraizinho <- read.csv("RegressãoQuantílica_STREET/icaraizinho.csv", header = FALSE)[,1]
X <- lagmatrix(icaraizinho,12)
pasta.origem <- "RegressãoQuantílica_STREET/"
# pasta.destino <- "Documento Regressao Quantilica/Figuras/selecao-lasso/"
nomes <- list.files(path = pasta.origem, pattern = "table-betas-selecaolasso-")
# nomes <- list.files(path = pasta.origem, pattern = "table-betas-sellassonorm-")
P = 12
maxK = P
arquivos.tabelas.lasso <- str_c(pasta.origem,nomes)

# lambdas = exp(seq(-2,10,0.3))  # Lasso sem normalizar
# lambdas = exp(seq(-2,9,0.1))


i=1
for (i in 1:length(arquivos.tabelas.lasso)) {
  tabela <- read.table(arquivos.tabelas.lasso[i], header = FALSE, sep = ",")
  vars <- abs(tabela[-1,]) > 0.0000001
  rownames(vars) <- str_c("beta_", 1:12, "")
  # colnames(vars) <- lambdas
  
  table.results <- matrix(0,nrow(tabela), ncol(tabela))
  # extracts alpha from name of table
  tmp <- str_extract(nomes[i], pattern = "\\d+")
  tmp.strlen <- str_length(tmp)
  tmp.parts <- str_sub(tmp, c(1,2), c(1,tmp.strlen))
  alpha <- as.numeric(str_c(tmp.parts[1], ".", tmp.parts[2]))
  for (j in 1:ncol(tabela)) {
    include.model <- vars[, j]
    # fit.rq <- rq(formula = icaraizinho[13:372] ~ X[, include.model], tau = 0.05)
    if (any(include.model)) {
      fit.rq <- rq(formula = icaraizinho[13:372] ~ X[, include.model], tau = alpha)
    } else {
      fit.rq <- rq(formula = icaraizinho[13:372] ~ 1, tau = alpha)
    }  
    table.results[c(TRUE,include.model), j] <- as.vector(coef(fit.rq))
  }
  alpha_out <- str_replace(alpha, "\\.", "")
  write.table(table.results, file = str_c(pasta.origem,"table-betas-postlasso-alpha-", alpha_out, ".csv"), sep= ",", row.names = FALSE, col.names = FALSE)
}



#  include.model <- (beta.lasso > 0)[-1]
#  fit.rq <- rq(formula = icaraizinho[13:372] ~ X[, include.model], tau = 0.9)
