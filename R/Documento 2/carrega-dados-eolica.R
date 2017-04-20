setwd('/home/mcruas/Dropbox/Pesquisa Doutorado/Paper-NPQuantile')
source('R/biblioteca-funcoes-npquantile.R')

ipak(c('readxl', 'dplyr', 'lattice'))
dados <- read_excel(path = 'Dados Climaticos/Solar-tubarao/tubarao solar.xlsx')[,1:6]
dados_filtrados <- dados %>% select(yt0, yt1, hora, mes) %>% as.data.frame
# colnames(dados_filtrados) = NULL
class(dados_filtrados)
attach(dados_filtrados)
tmp = matrix(dados_filtrados %>% as.vector, nrow(dados_filtrados))
yt0 = dados_filtrados$yt0
yt1 = dados_filtrados$yt1
hora = dados_filtrados$hora
mes = dados_filtrados$mes

hist(yt0[hora == ])
# boxplot(yt0 ~ hora, data = dados_filtrados, col = 2)



# bwplot(hora ~ yt0 | factor(mes), data = dados_filtrados, horizontal = TRUE)


#
# ipak(c('dplyr', 'stringr', 'readxl', 'lattice', 'lubridate'))
# localidade = 'saojosemipibu' # Escolher entre saojosemipibu ; tubarao ; tabocas
#
# caminho_arquivo = sprintf('Dados Climaticos/Solar-%s/%s.xlsx',localidade, localidade) # caminho do arquivo contendo os dados
# serie <- read_excel(path = caminho_arquivo, sheet = 'original data') %>% mutate(hora = hour(UTC))
# (serie %>% select(kW))[[1]] %>% plot
