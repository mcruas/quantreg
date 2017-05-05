source("R/biblioteca-funcoes-npquantile.R")

ipak(c("readxl", "dplyr", "lattice"))
dados <- read_excel(path = "Dados Climaticos/Solar-tubarao/tubarao solar-28.467_-49.005_uncorrected.xlsx")[,1:6]
dados.filtrados <- dados %>% filter(y_t0 > 0 , y_t1 > 0, hour > 9, hour < 18) %>% arrange(y_t0) %>% select(y_t0, y_t1)
boxplot(y_t0 ~ month, data = dados.filtrados)

xyplot(y_t0 ~ hour, data = dados, type = c("p","g"))
boxplot(y_t0 ~ hour, data = dados, col = "gray", main = "Hourly predicted solar power", xlab = "Hour")


boxplot(y_t0 ~ month, data = dados, col = "gray", main = "Monthly predicted solar power", xlab = "Hour")



xyplot(y_t0 ~ hour | factor(month), data = dados, type = c("p","g"))
densityplot(y_t0 ~ hour | factor(month), data = dados, type = "l")
densityplot( ~  y_t0 | factor(hour), data = dados.filtrados, type = "l")

for (i in 10:16) {
  dado.mensal <- dados %>% filter(hour == i) %>% select(y_t0) 
  plot(density(dado.mensal[[1]]))
}


plot(dados.filtrados$y_t0, dados.filtrados$y_t1, col = rainbow(12)[dados.filtrados$month], pch = 16)
legend("topleft",legend = 1:12, cex = 0.5)


densityplot(y_t0 ~ `time UTC`, group = factor(month), data = dados, type = "l", auto.key = list(column = 5))



xyplot(y_t0 ~ y_t1, groups = factor(month), data = dados.filtrados, auto.key = list(column = 5))

gray(1:100/100)



# Gerar dados otimização --------------------------------------------------
vetor_horas <- 0:23; i = 12
hora <- vetor_horas[i]
dados_por_hora <- dados %>% filter(y_t0 > 0 , y_t1 > 0, hour == hora) %>% 
  arrange(y_t0) %>% select(y_t0, y_t1)
x = dados_por_hora$y_t0
y = dados_por_hora$y_t1
plot(x,y)




  	library('dplyr')
		library('readxl')
		library('lattice')
		# dados <- read_excel(path = '../Dados Climaticos/Solar-tubarao/tubarao solar-28.467_-49.005_uncorrected.xlsx')[,1:6]
 		dados_por_hora <- dados %>% filter(y_t0 > 0 , y_t1 > 0, hour == hora)
		n = nrow(dados_por_hora)
		dados_por_hora %>% mutate(y_t0modif = y_t0 + rnorm(n)*10) %>% arrange(y_t0modif) %>% select(y_t0modif, y_t1)
  	x = dados_por_hora$y_t0modif
 		y = dados_por_hora$y_t1


