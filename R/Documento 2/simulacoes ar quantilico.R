source("R/biblioteca-funcoes-npquantile.R")
hora = 8
hora_ini = 0:23

ipak(c("quantreg", "dplyr", "stringr", "readxl", "lattice", "lubridate"))
localidade = "saojosemipibu" # Escolher entre saojosemipibu ; tubarao ; tabocas

caminho_arquivo = sprintf("Dados Climaticos/Solar-%s/%s.xlsx",localidade, localidade) # caminho do arquivo contendo os dados
serie <- read_excel(path = caminho_arquivo, sheet = "original data") %>% mutate(hora = hour(UTC))
dados_por_hora_ini <- dados %>% filter(hour %in% hora_ini) %>% select(y_t0)
dados_por_hora_prev <- dados %>% filter(hour %in% hora_prev) %>% select(y_t0)
# dados_por_hora <- dados %>% filter(y_t0 > 0 , y_t1 > 0, hour %in% hora)

n = nrow(dados_por_hora)
dados_por_hora =	dados_por_hora %>% mutate(y_t0modif = y_t0 + runif(n)*0.00001) %>% arrange(y_t0modif) %>% select(y_t0modif, y_t1)
x = dados_por_hora$y_t0modif
y = dados_por_hora$y_t1


qr1 <- rq(y ~ x, tau = 0.5)
plot(x,y)
abline(qr1,col=2)


qr1 <- rq(y ~ x, tau = seq(0.01,0.99,0.01))
