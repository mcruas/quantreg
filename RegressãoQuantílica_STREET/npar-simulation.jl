using JuMP, DataFrames, Plots, RCall, Interpolations, Dierckx #, Distributions

n = 1000
gr()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
cd(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/")
# cd("C:/Users/mcruas/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET")
# pwd()
pwd();
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/par-multi.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");

# @rput n
# R"
# # serie = arima.sim(n = n, list(ar = c(0.9,-0.2), ma = NULL),
# #           sd = sqrt(0.1796))
# source('R/biblioteca-funcoes-npquantile.R')
#
# ipak(c('readxl', 'dplyr', 'lattice'))
# dados <- read_excel(path = 'Dados Climaticos/Solar-tubarao/tubarao solar.xlsx')[,1:6]
# dados_filtrados <- dados %>% select(yt0, yt1, hora, mes) %>% as.matrix
# # boxplot(yt0 ~ hora, data = dados.filtrados, col = 2)
# "
# @rget dados

serie_tmp = readcsv("Dados Climaticos/Enas/enas sudeste 1931-2014.csv"); nomeserie = "enasudeste"
serie = serie_tmp ./ 100;



max_sim= 10
n_cenarios = 100
n_tmp = length(serie) # keeps size of series before cutting them
### Inicialização

# Does adaptation for a AR(2) model

X = [rand(n_tmp-2) serie[2:n_tmp-1] serie[1:n_tmp-2]] # X = [y_t-1 y_t-2]
X = [serie[2:n_tmp-1] serie[1:n_tmp-2]]

y = serie[3:n_tmp]

n = length(y)
T = 1:n

Alphas = collect(0.005:0.005:0.995)

Y = [ y * ones(1,n_cenarios);
      zeros(max_sim,n_cenarios)]


### Step 1
tau = n +1


### Step 2
X_tau = [y[tau-1] ; y[tau-2]] #não é generico. nao aceitar.

betas0, betas = rq_par(y,X, Alphas)
# plot(betas')
Q_hat = (betas0 + X_tau' * betas)[1,:]

### Step 3
unif = rand(n_cenarios)
y_simulated = Q(unif, Q_hat, Alphas)
Y[tau, :] = y_simulated

### Step 4
for s in 1:n_cenarios, tau in n+2:n+max_sim

  # (step 2)
  X_tau = [Y[tau-1,s] ; Y[tau-2,s]] #não é generico. nao aceitar.
  # betas0, betas = rq(y,X, Alphas)
  Q_hat = (betas0 + X_tau' * betas)[1,:]

  # (step 4)
  unif = rand(1)
  y_simulated = Q(unif, Q_hat, Alphas)[1,1]
  Y[tau, s] = y_simulated


end
