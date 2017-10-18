library(stringr, dplyr); library(xtable)
pasta.origem <- "RegressãoQuantílica_STREET/"
pasta.destino <- "Documento Regressao Quantilica/Tabelas/"
nomes <- list.files(path = pasta.origem, pattern = "table-")
arquivos.tabelas <- str_c(pasta.origem,nomes)
for (i in 1:length(arquivos.tabelas)) {
  tmp <- read.table(arquivos.tabelas[i], header = FALSE, sep = ",")
  rownames(tmp) <- str_c("$\\beta_{", 0:12, "}$")
  colnames(tmp) <- str_c("K=",1:12)
  tabela.latex <- xtable(tmp)
  print(tabela.latex, type = "latex", file = str_c(pasta.destino, 
          str_replace(nomes[i], ".csv", ".tex")), sanitize.text.function=function(x){x})
}


# \caption{\textbf{$\alpha = 0.05.} Coefficients $\beta_i$ for each value of $K$, where $K$ is the number of nonzero coefficients, excluding the intercept $\beta_0 which is always included.$}



