# include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/experimentos_grupo.jl")

#######################################################
# Este script gera gráficos para o relatório, mostrando como uma variação
# em x_t afeta o quantil do período seguinte.

################# Gráficos linear #####################
vetor_x = [-1 0 1 2 3 10] ; n = 100
# X_tau = 500 # Qual o valor de X utilizado para a previsão
# limx = 1500; limy = 2000; # para ENA sudeste
limx = 55; limy = 55; # para icaraizinho
# limx = NaN; limy = NaN # quando não se sabe

using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx #, Distributions



usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

#####################################################
############# Carregar dados em R ########################
pwd();
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
# include(pwd()*"/RegressãoQuantílica_STREET/par-multi.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
n = 100;
# @rput n pasta_trabalho;
# R"
# serie = arima.sim(n = n, list(ar = c(0.9,-0.2), ma = NULL),
#           sd = sqrt(0.1796))
# setwd(pasta_trabalho)
# source('R/biblioteca-funcoes-npquantile.R') # Carrega dados por script do R
# # source('R/biblioteca-funcoes-npquantile.R')
# # ipak(c('readxl', 'dplyr', 'lattice'))
# # dados <- read_excel(path = 'Dados Climaticos/Solar-tubarao/tubarao solar.xlsx')[,1:6]
# # dados_filtrados <- dados %>% select(yt0, yt1, hora, mes) %>% as.matrix
# # boxplot(yt0 ~ hora, data = dados.filtrados, col = 2)
# "
#
# @rget serie dados_filtrados
# serie = dados_filtrados[:,1].values;
################################################################

############# Carregar dados ###################################

# para icaraizinho
serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"; x_new = 30
vetor_X_tau = collect(10:8:42)
# serie = readcsv("Dados Climaticos/Enas/enas sudeste 1931-2014.csv"); nomeserie = "enasudeste"; serie = serie ./ 100; x_new = 500 # ENA Sudeste
# vetor_X_tau = collect(200:50:700) # para ENA sudeste


# serie = readcsv("Dados Climaticos/Solar-tabocas/tabocas.csv"); nomeserie = "solartabocas"; x_new = 0.4 # Para Tabocas


n_tmp = length(serie); # keeps size of series before cutting them
### Inicialização

# Does adaptation for a AR(2) model

# X = [serie[2:n_tmp-1] serie[1:n_tmp-2]]; # AR(2)

# X = Array{Float64,2} # AR(1)
tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X = zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]

y = serie[3:n_tmp];
n = length(y);
T = 1:n;

# X_lags = lagmatrix(serie,0:2)
# Alphas = collect(0.005:0.005:0.995);


# Y = [ y * ones(1,n_cenarios);
#       zeros(max_sim,n_cenarios)];

limx = isnan(limx) ? maximum(serie) : limx
limy = isnan(limy) ? maximum(serie) : limy

plot(serie)


################################################################################################
################################################################################################
# Experimentos com o seleção inteira
X_lags = lagmatrix(serie,0:12)
Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
TimeLimit = 600

betas0, betas = rq_par_mip(X_lags[:,1], X_lags[:, 2:end], Alphas; non_cross = true, max_K = 3, TimeLimit = TimeLimit, MIPGap = 0.00)

betas0, betas = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
          non_cross = true, max_K = 3, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = 1)

heatmap(Alphas, betas .== 0)
savefig("Coeficientes MIP apos 6000s.pdf")
