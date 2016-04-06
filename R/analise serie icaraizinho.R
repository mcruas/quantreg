library(ggfortify, ggplot2)
x = read.csv("RegressãoQuantílica_STREET/icaraizinho.csv", header = FALSE)[,1]
plot.ts(x)
pacf(x, main = "Série Icaraizinho")
acf(x)
qplot(x[-length(x)],x[-1])

     
     lag(1:10)
     