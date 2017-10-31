library(ggplot2)
library(ggfortify)
x1 = seq(0,1,0.01)
x2 = seq(1,2,0.01)
x = c(x1,x2)
y1 = 1 + rnorm(101,sd = 0.05)
y2 = x2*1 + rnorm(101,sd = 0.05)
y = c(y1,y2)
plot(x, y)
library(quantreg)
fit = rq(y1 ~ x1,tau=0.2)
lines(x1, coef(fit)[1] + coef(fit)[2] * x1, type = "l", col = 2)
fit = rq(y2 ~ x2,tau=0.2)
lines(x2, coef(fit)[1] + coef(fit)[2] * x2, type = "l", col = 2)

fit = rq(y ~ x,tau=0.2)
qplot(x, y) + geom_abline(intercept = coef(fit)[1], slope  =  coef(fit)[2], colour = "red")

write.table(x, file = "../RegressãoQuantílica_STREET/x_nonlin.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(y, file = "../RegressãoQuantílica_STREET/y_nonlin.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# Arima
ts.sim <- arima.sim(list(order = c(1,0,3), ar = 0.7, ma = c(0.2, 0.6, 0.3)), n = 200)
ts.plot(ts.sim)

library(ggfortify)
autoplot(ts.sim)
qplot(ts.sim[-length(ts.sim)],ts.sim[-1])


yt.ini = arima.sim(list(order = c(4,0,3), ar = c(0.3,0.02,0.02, 0.3), ma = c(0.2, 0.6, 0.3)), n = 200)
# plot(yt.ini)
xt = yt.ini[-length(yt.ini)]; yt = yt.ini[-1]
ordem = order(xt)
x = xt[ordem]
y = yt[ordem]
# plot(x,y)


xt = rnorm(500); yt=rnorm(500)
ordem = order(xt)
x = xt[ordem]
y = yt[ordem]
plot(x,y)



# Gráfico coeficientes a partir das tabelas

library(stringr, dplyr); library(xtable)
pasta.origem <- "RegressãoQuantílica_STREET/"
pasta.destino <- "Documento Regressao Quantilica/Figuras/selecao-lasso/"
nomes <- list.files(path = pasta.origem, pattern = "table-betas-selecaointeira-")
arquivos.tabelas <- str_c(pasta.origem,nomes)
i=1
for (i in 1:length(arquivos.tabelas)) {
  tmp <- read.table(arquivos.tabelas[i], header = FALSE, sep = ",")
  rownames(tmp) <- str_c("$\\beta_{", 0:12, "}$")
  colnames(tmp) <- str_c("K=",1:12)
  x11()
  plot(tmp)
  plot(ts(tmp[, 1:10]))
  print(tabela.latex, type = "latex", file = str_c(pasta.destino, 
                                                   str_replace(nomes[i], ".csv", ".tex")), sanitize.text.function=function(x){x})
}


