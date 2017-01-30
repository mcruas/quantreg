
# setwd("/Users/henrique/Desktop/Tese GAS/Resultados univariados/ENANCEL_COL_16")
# setwd("/Users/henrique/Desktop/Tese GAS/Resultados univariados/ICARAIZINHO_COL_4")
# setwd("/Users/henrique/Desktop/Tese GAS/Resultados univariados/PRAIA DO MORGADO_COL_13")

########################################################################################################
bd=read.csv("RegressãoQuantílica_STREET/icaraizinho.csv", header = FALSE,sep=";",dec=".")

#bd=read.csv("novosdados2.csv", header = FALSE,sep=";",dec=",")

#rm(list = ls())
##########################################################################################################
# Criando o banco de dados m?s a m?s
##########################################################################################################
prod.mensal.certif=data.matrix(bd)
bd.mes.a.mes = matrix(NA,31,12)
alpha=0
for(m in 1:12){
  cont=1
  for (a in 1:31){
    bd.mes.a.mes[a,m]=prod.mensal.certif[cont+alpha,1]
    cont=cont+12 
  }
  alpha=alpha+1
}

bd.mes.a.mes
colnames(bd.mes.a.mes) = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
rownames(bd.mes.a.mes) = 1981:2011
pdf("Documento Regressao Quantilica/Figuras/Icaraizinho/icaraizinho-boxplot.pdf")
boxplot(bd.mes.a.mes,col =2, ylab = "Mean MW", main = "Icaraizinho monthly data boxplot")
dev.off()


meltdf = melt(bd.mes.a.mes); colnames(meltdf) = c("Year", "Month", "Mean MW")

require(reshape2); require(ggplot2)
ggplot(data = meltdf, aes(x=Month, y=`Mean MW`)) +
  geom_boxplot(aes(fill=Month)) + theme(legend.position="none") + xlab("")

ggplot(meltdf,aes(x=Month,y=`Mean MW`,colour=Year,group=Year)) + geom_line() + xlab("")


# Faz tabela de icaraizinho mes a mes
library(xtable)
print(xtable(bd.mes.a.mes), type = "latex", file = "Documento Regressao Quantilica/Tabelas/tabela-icaraizinho.tex")

## scatterplot mensal em icaraizinho (vento)
library(latex2exp)
pdf("Documento Regressao Quantilica/Figuras/Icaraizinho/scatterplot.pdf", width = 4, height = 4)
# x11()
qplot(bd[-1,],head(bd,-1)) + xlab(TeX("y_{t-1}")) + ylab(TeX("y_t")) + ggtitle("Icaraizinho monthly data")
dev.off()




################# Mesmo gráficos com série solar
library(readxl); library(dplyr)
bd2 <- read_excel("Dados Climaticos/Solar-tubarao/tubarao solar.xlsx")[, 1:8] %>% select(yt0, yt1,mes)

pdf("Documento Regressao Quantilica/Figuras/Solar-exemplos/tubarao-boxplot.pdf")
boxplot(yt0 ~ mes, data = bd2, col = 2, names = c("Jan","Feb","Mar","Apr",
                          "May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), ylab = "Mean MW",
                          main = "Tubarão hourly data boxplot")
dev.off()
serie = bd2[[1]]

pdf("Documento Regressao Quantilica/Figuras/Solar-exemplos/scatterplot.pdf", width = 3, height = 3)
# x11()
qplot(bd2$yt0,bd2$yt1) + xlab(TeX("y_{t-1}")) + ylab(TeX("y_t")) + ggtitle("Tubarão hourly data")
dev.off()

# bd.mes.a.mes <- bd %>% group_by(mes) %>% summarise(blau = mean(yt0, na.rm = TRUE))


## scatterplot diário em tubarão (solar)










##########################################################################################################
# Gráfico Time serie
##########################################################################################################
#y.EN = prod.mensal.certif
#y.PM = prod.mensal.certif
#y.IC = prod.mensal.certif

length(y.EN)
y.plot.EN = ts(y.EN,start=c(1981,1),end=c(2011,12),frequency=12)

plot.ts(y.plot.EN, type="l", col=4, ylab="FC", main="Enacel",ylim=c(0,100),xlab="Anos")
abline(v = 348, col = 2,lty=2,lwd = 4)

setEPS()
postscript("serie_EN.eps")
plot.ts(y.plot.EN, type="l", col=4, ylab="FC", main="Enacel",ylim=c(0,100),xlab="Anos")
abline(v = 348, col = 2,lty=2,lwd = 4)
dev.off()


y.plot.IC = ts(y.IC,start=c(1981,1),end=c(2011,12),frequency=12)

plot.ts(y.plot.IC, type="l", col=4, ylab="FC", main="Icaraizinho",ylim=c(0,100),xlab="Anos")
abline(v = 348, col = 2,lty=2,lwd = 4)

setEPS()
postscript("serie_IC.eps")
plot.ts(y.plot.EN, type="l", col=4, ylab="FC", main="Icaraizinho",ylim=c(0,100),xlab="Anos")
abline(v = 348, col = 2,lty=2,lwd = 4)
dev.off()


y.plot.PM = ts(y.IC,start=c(1981,1),end=c(2011,12),frequency=12)

plot.ts(y.plot.PM, type="l", col=4, ylab="FC", main="Praia do Morgado",ylim=c(0,100),xlab="Anos")
abline(v = 348, col = 2,lty=2,lwd = 4)

setEPS()
postscript("serie_PM.eps")
plot.ts(y.plot.PM, type="l", col=4, ylab="FC", main="Praia do Morgado",ylim=c(0,100),xlab="Anos")
abline(v = 348, col = 2,lty=2,lwd = 4)
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#1 = EN
#2 = PM
#3 = IC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#bd.out.1 = ts(c(prod.mensal.certif[361:length(prod.mensal.certif)])) #EN
#bd.out.2 = ts(c(prod.mensal.certif[361:length(prod.mensal.certif)])) #PM
#bd.out.3 = ts(c(prod.mensal.certif[361:length(prod.mensal.certif)])) #IC



##########################################################################################################
# Gráfico BOX Plot mensal
##########################################################################################################

library(ggplot2)
library(reshape2)
matriz.GF.2=bd.mes.a.mes
colnames(bd.mes.a.mes) = c("Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez")

#---------------------------------------------------------------------------------------------------------------#
# Box-plots das series mensais
#---------------------------------------------------------------------------------------------------------------#
setEPS()
postscript("BOXPLOT_PM.eps")
teste<-factor(colnames(bd.mes.a.mes),levels=c("Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez"))
a.plotar<-melt(as.data.frame(matriz.GF.2),id.vars=colnames(colnames(matriz.GF.2)))
ggplot(a.plotar,aes(variable,value))+ geom_boxplot() + coord_cartesian(ylim = c(0, 100))+ scale_y_continuous(name="FC",breaks=seq(0, 100, 10))+ scale_x_discrete(name="Meses")+ theme(axis.title.x = element_text(face="bold", size=15),
                                                                                                                                                                                            axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
                                                                                                                                                                                              axis.title.y = element_text(face="bold", size=15),
                                                                                                                                                                                              axis.text.y  = element_text(size=16)) + labs(title = "Box-Plot Mensal - Praia do Morgado",size=40)
dev.off()


#---------------------------------------------------------------------------------------------------------------#
# Time series plots das series mensais
#---------------------------------------------------------------------------------------------------------------#
{
# setEPS()
# postscript("Grafico_SERIES.eps")
# 
# tdat <- data.frame(Site = rep(paste0("Usina ", c("Enacel","Icaraizinho","Praia do Morgado")),
#                               each = 372),
#                    Date = rep(seq(Sys.Date(), by = "1 day", length = 372), 3),
#                    Fitted = c(y.EN,y.IC,y.PM),
#                    Signif = rep(NA, 372*3))
# 
# ## select 1 region per Site as signif
# tdat$Signif[360:372] <- tdat$Fitted[360:372]
# tdat$Signif[(360+372):(372+372)] <- tdat$Fitted[(360+372):(372+372)]
# tdat$Signif[(360+372*2):(372+372*2)] <- tdat$Fitted[(360+372*2):(372+372*2)]
# 
# 
# ggplot(tdat, aes(x = Date, y = Fitted, group = Site)) +
# geom_line() + 
# labs(size= "Nitrogen",
#        x = "Data",
#        y = "FC eólico",
#        title = "Gráfico das séries das três usinas")+
# #geom_line(mapping = aes(y = Upper), lty = "dashed") +
# #geom_line(mapping = aes(y = Lower), lty = "dashed") +
# geom_line(mapping = aes(y = Signif), lwd = 1.3, colour = "red") +
# facet_wrap( ~ Site)
# 
# dev.off()
}

tdat <- data.frame(Site = rep(paste0("Usina ", c("Enacel")),
                              each = 372),
                   Date = rep(seq(as.Date(c("1981-01-01")), by = "1 month", length = 372), 1),
                   Fitted = c(y.EN),
                   Signif = rep(NA, 372))

## select 1 region per Site as signif
tdat$Signif[360:372] <- tdat$Fitted[360:372]
#tdat$Signif[(360+372):(372+372)] <- tdat$Fitted[(360+372):(372+372)]
#tdat$Signif[(360+372*2):(372+372*2)] <- tdat$Fitted[(360+372*2):(372+372*2)]

setEPS()
postscript("Grafico_SERIES_EN.eps")

ggplot(tdat, aes(x = Date, y = Fitted, group = Site)) +
  geom_line() + scale_y_continuous(name="FC",breaks=seq(0, 100, 10))+coord_cartesian(ylim = c(0, 100))+
  labs(size= "Nitrogen",
       x = "Ano",
       y = "FC eólico",
       title = "")+
  theme(axis.title.x = element_text(face="bold", size=15),
                                          axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
                                            axis.title.y = element_text(face="bold", size=15),
                                            axis.text.y  = element_text(size=16))+
  
  #geom_line(mapping = aes(y = Upper), lty = "dashed") +
  #geom_line(mapping = aes(y = Lower), lty = "dashed") +
  geom_line(mapping = aes(y = Signif), lwd = 1.3, colour = "red") +
  facet_wrap( ~ Site)

dev.off()

#---------------------------------------------------------------------------------------------------

tdat <- data.frame(Site = rep(paste0("Usina ", c("Icaraizinho")),
                              each = 372),
                   Date = rep(seq(as.Date(c("1981-01-01")), by = "1 month", length = 372), 1),
                   Fitted = c(y.IC),
                   Signif = rep(NA, 372))

## select 1 region per Site as signif
tdat$Signif[360:372] <- tdat$Fitted[360:372]
#tdat$Signif[(360+372):(372+372)] <- tdat$Fitted[(360+372):(372+372)]
#tdat$Signif[(360+372*2):(372+372*2)] <- tdat$Fitted[(360+372*2):(372+372*2)]
setEPS()
postscript("Grafico_SERIES_IC.eps")

ggplot(tdat, aes(x = Date, y = Fitted, group = Site)) +
  geom_line() + scale_y_continuous(name="FC",breaks=seq(0, 100, 10))+coord_cartesian(ylim = c(0, 100))+
  labs(size= "Nitrogen",
       x = "Ano",
       y = "FC eólico",
       title = "")+
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
          axis.title.y = element_text(face="bold", size=15),
          axis.text.y  = element_text(size=16))+
  
  #geom_line(mapping = aes(y = Upper), lty = "dashed") +
  #geom_line(mapping = aes(y = Lower), lty = "dashed") +
  geom_line(mapping = aes(y = Signif), lwd = 1.3, colour = "red") +
  facet_wrap( ~ Site)

dev.off()
#---------------------------------------------------------------------------------------------------

tdat <- data.frame(Site = rep(paste0("Usina ", c("Praia do Morgado")),
                              each = 372),
                   Date = rep(seq(as.Date(c("1981-01-01")), by = "1 month", length = 372), 1),
                   Fitted = c(y.PM),
                   Signif = rep(NA, 372))

## select 1 region per Site as signif
tdat$Signif[360:372] <- tdat$Fitted[360:372]
#tdat$Signif[(360+372):(372+372)] <- tdat$Fitted[(360+372):(372+372)]
#tdat$Signif[(360+372*2):(372+372*2)] <- tdat$Fitted[(360+372*2):(372+372*2)]
setEPS()
postscript("Grafico_SERIES_PM.eps")

ggplot(tdat, aes(x = Date, y = Fitted, group = Site)) +
  geom_line() + scale_y_continuous(name="FC",breaks=seq(0, 100, 10))+coord_cartesian(ylim = c(0, 100))+
  labs(size= "Nitrogen",
       x = "Ano",
       y = "FC eólico",
       title = "")+
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
          axis.title.y = element_text(face="bold", size=15),
          axis.text.y  = element_text(size=16))+
  
  #geom_line(mapping = aes(y = Upper), lty = "dashed") +
  #geom_line(mapping = aes(y = Lower), lty = "dashed") +
  geom_line(mapping = aes(y = Signif), lwd = 1.3, colour = "red") +
  facet_wrap( ~ Site)

dev.off()



#---------------------------------------------------------------------------------------------------------------#
# Scatterplot entre as séries
#---------------------------------------------------------------------------------------------------------------#

tdat_PM_IC <- data.frame(Site = rep(paste0("Usina ", c("Praia do Morgado", "Icaraizinho")),each = 372),
                   xvar = c(y.PM),
                   yvar = c(y.IC))
                   
setEPS()
postscript("Scater_IC_PM.eps")

ggplot(tdat_PM_IC, aes(x=xvar, y=yvar)) +
  geom_point(shape=1)  +
  scale_y_continuous(name="FC Icaraizinho",breaks=seq(0, 85, 10))+
  coord_cartesian(ylim = c(0, 85))+
  
  scale_x_continuous(name="FC Praia do Morgado",breaks=seq(0, 70, 10))+
  coord_cartesian(xlim = c(0, 70))+
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
          axis.title.y = element_text(face="bold", size=15),
          axis.text.y  = element_text(size=16))+
  geom_smooth(method=lm)   

dev.off()

#---------------------------------------------------------------------------------------------------------------#
tdat_PM_EN <- data.frame(Site = rep(paste0("Usina ", c("Praia do Morgado", "Enacel")),each = 372),
                         xvar = c(y.PM),
                         yvar = c(y.EN))

setEPS()
postscript("Scater_PM_EN.eps")

ggplot(tdat_PM_EN, aes(x=xvar, y=yvar)) +
  geom_point(shape=1)  +
  scale_y_continuous(name="FC Enacel",breaks=seq(0, 70, 10))+
  coord_cartesian(ylim = c(0, 70))+
  
  scale_x_continuous(name="FC Praia do Morgado",breaks=seq(0, 70, 10))+
  coord_cartesian(xlim = c(0, 70))+
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
          axis.title.y = element_text(face="bold", size=15),
          axis.text.y  = element_text(size=16))+
  geom_smooth(method=lm)   

dev.off()


#---------------------------------------------------------------------------------------------------------------#

tdat_IC_EN <- data.frame(Site = rep(paste0("Usina ", c("Icaraizinho", "Enacel")),each = 372),
                         xvar = c(y.IC),
                         yvar = c(y.EN))

setEPS()
postscript("Scater_IC_EN.eps")

ggplot(tdat_IC_EN, aes(x=xvar, y=yvar)) +
  geom_point(shape=1)  +
  scale_y_continuous(name="FC Enacel",breaks=seq(0, 70, 10))+
  coord_cartesian(ylim = c(0, 70))+
  
  scale_x_continuous(name="FC Icaraizinho",breaks=seq(0, 85, 10))+
  coord_cartesian(xlim = c(0, 85))+
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))+theme(
          axis.title.y = element_text(face="bold", size=15),
          axis.text.y  = element_text(size=16))+
  geom_smooth(method=lm)   

dev.off()


#---------------------------------------------------------------------------------------------------------------#
# Scatterplot entre as PITS
#---------------------------------------------------------------------------------------------------------------#

dados_matriz_Pits <- data.frame(c(pit.RF), c(pit.IC), c(pit.EN))
colnames(dados_matriz_Pits) =c("Rio do Fogo", "Icaraizinho", "Enacel")
pairs(cbind(pit.RF, pit.IC, pit.EN))

setEPS()
postscript("Matriz_PITS.eps")
plotmatrix(dados_matriz_Pits, colour="gray20")
dev.off()







