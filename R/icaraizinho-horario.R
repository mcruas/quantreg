library(lubridate)
library(dplyr)

# Reads data and creates new columns for month, day, year and hour
transform.data <- function(df) {
  df$time <- ymd_hms(df$FIELDS)
  df$day <- day(df$time)
  df$year <- year(df$time)
  df$month <- month(df$time)
  df$hour <- hour(df$time)
  df$FIELDS <- NULL
  return(df)
}

ERAI <- read.csv("~/Dropbox/Pesquisa Doutorado/Paper NPQuantile/R/ICARAIZINHO/0-UEE0_ICARAIZINHO-ERAI.csv")
MERRA <- read.csv("~/Dropbox/Pesquisa Doutorado/Paper NPQuantile/R/ICARAIZINHO/0-UEE0_ICARAIZINHO-MERRA.csv")
NNRP <- read.csv("~/Dropbox/Pesquisa Doutorado/Paper NPQuantile/R/ICARAIZINHO/0-UEE0_ICARAIZINHO-NNRP.csv")

ERAI <- transform.data(ERAI)
MERRA <- transform.data(MERRA)
NNRP <- transform.data(NNRP)

plot(ERAI$windspeed_100m, lag(ERAI$windspeed_100m, n = 1))
plot(ERAI$windspeed_100m, lag(ERAI$windspeed_100m, n = 12))
plot(ERAI$windspeed_100m, lag(ERAI$windspeed_100m, n = 24))

plot(ERAI$windspeed_100m, lag(ERAI$windspeed_100m, n = 4320))

hist(ERAI$windspeed_100m)
hist(MERRA$windspeed_100m)
hist(NNRP$windspeed_100m)

plot(ERAI$windspeed_100m, ERAI$winddirection_100m)
plot(MERRA$windspeed_100m, MERRA$winddirection_100m)

pairs(ERAI[, c(1,2,3,4,8)])

library(car)
scatterplotMatrix(ERAI[, c(1,2,3,4,8)],
                  diagonal="histogram",
                  smooth=TRUE)

# Carrega pacotes e leitura inicial dos dados -----------------------------

library(ggfortify)
library(lubridate)
library(dplyr)
library(stringr)

RealDataES <- read.csv("Dados complatt/RealDataES.csv", na.strings="n.a.", sep = ";")
RealDataPT <- read.csv("Dados complatt/RealDataPT.csv", na.strings="n.a.")
RealMarketPriceDataPT <- read.csv("Dados complatt/RealMarketPriceDataPT.csv", sep = ";")
historical_weather <- read.csv("Dados complatt/historical_weather.csv")


# Transformação dos tipos das datas em POIXlt, por meio das funções do pacote 'lubridate'
RealDataPT$date <- ymd_hms(RealDataPT$date)
RealDataES$data..UTC. <- dmy_hm(as.character(RealDataES$data..UTC.))
RealMarketPriceDataPT$date..UTC. <- dmy_hm(as.character(RealMarketPriceDataPT$date..UTC.))
historical_weather$prediction_date <- dmy_hm(as.character(historical_weather$prediction_date))

names(RealDataPT) <- str_c(names(RealDataPT),"_PT")
names(RealDataES) <- str_c(names(RealDataES),"_ES")


dados <- RealDataPT %>% select(-id_PT) %>% inner_join(RealDataES, by = c("date_PT" = "data..UTC._ES"))
dados <- dados %>% inner_join(RealMarketPriceDataPT, by = c("date_PT" = "date..UTC."))

# para cada ponto de observações dos dados de clima, cria uma coluna distinta
for (i in unique(historical_weather$point)) { 
  tmp <- historical_weather %>% filter(point == i) %>% select(-point)
  names(tmp)[-1] <- str_c( names(tmp)[-1], "_p" , i) # adiciona o valor do ponto no nome da coluna
  dados <- dados %>% inner_join(tmp, by = c("date_PT" = "prediction_date"))
}

write.csv(dados, file = "dados_compilados.csv", row.names = FALSE)

