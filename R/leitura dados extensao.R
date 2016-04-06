extJRAI_FC <- read.csv("R/Dados extensão/JRAI/extJRAI_FC.csv")
names(extJRAI_FC) <- c("Ano", "Mes", "Dia", "Hora", "value")

yt.ini <- extJRAI_FC["value"]
xt = yt.ini[-length(yt.ini)]; yt = yt.ini[-1]
ordem = order(xt)
x = xt[ordem]
y = yt[ordem]
plot(x,y)
