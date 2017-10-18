source("R/biblioteca-funcoes-npquantile.R")
hora = 8
hora_ini = 0:23

ipak(c("quantreg", "dplyr", "stringr", "readxl", "lattice", "lubridate"))
localidade = "saojosemipibu" # Escolher entre saojosemipibu ; tubarao ; tabocas

caminho_arquivo = sprintf("Dados Climaticos/Solar-%s/%s.xlsx",localidade, localidade) # caminho do arquivo contendo os dados
serie <- read_excel(path = caminho_arquivo, sheet = "original data") %>% mutate(hora = hour(UTC))
