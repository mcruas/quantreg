source("R/biblioteca-funcoes-npquantile.R")
ipak(c("dplyr", "ggfortify", "lattice", "readr", "lubridate"))

dados <- read_csv("Dados Climaticos/Dados-kaggle/windforecasts_wf3.csv")
dados <- dados %>% mutate(date = ymd(floor(date/100))) # Changes date to lubridate format
dados <- dados %>% mutate(mes = month(date), ano = year(date), dia = day(date), hors = as.double(hors)) # Includes year, month and day columns
table(dados$mes)
table(dados$ano)


boxplot(ws ~ mes, data = dados, col = "gray", main = "Monthly predicted solar power", xlab = "Hour")

boxplot(ws ~ hors, data = dados, col = "gray", main = "Hourly predicted solar power", xlab = "Hour")

dados.plot <- dados %>% filter(hors > 40) %>% select(mes,hors,ws)
dados.plot %>% summarise(group_by(mes,hors), ws= mean(ws))
bwplot(hors ~ ws | factor(mes), data = dados.plot, col = "gray", 
       main = "Hourly predicted solar power", xlab = "Hour", horizontal = FALSE)


xyplot(ws ~ hour | factor(month), data = dados, type = c("p","g"))
densityplot(ws ~ hour | factor(month), data = dados, type = "l")
densityplot( ~  ws | factor(hour), data = dados.filtrados, type = "l")
