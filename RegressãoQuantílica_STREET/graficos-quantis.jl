
#######################################################
# Este script gera gráficos para o relatório, mostrando como uma variação
# em x_t afeta o quantil do período seguinte.

################# Gráficos linear #####################
vetor_x = [-1 0 1 2 3 10]

using JuMP, DataFrames, Plots, RCall, Interpolations, Dierckx #, Distributions
gr()
usesolver = "mosek"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)
# cd("C:/Users/mcruas/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET")
# pwd()
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl")
include(pwd()*"/RegressãoQuantílica_STREET/par-multi.jl")
@rput n
R"
# serie = arima.sim(n = n, list(ar = c(0.9,-0.2), ma = NULL),
#           sd = sqrt(0.1796))
source('R/biblioteca-funcoes-npquantile.R')

ipak(c('readxl', 'dplyr', 'lattice'))
dados <- read_excel(path = 'Dados Climaticos/Solar-tubarao/tubarao solar.xlsx')[,1:6]
dados_filtrados <- dados %>% select(yt0, yt1, hora, mes) %>% as.matrix
# boxplot(yt0 ~ hora, data = dados.filtrados, col = 2)
"

@rget serie dados_filtrados
serie
plot(serie)

dados_filtrados[1,1] + 4.

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

betas0, betas = rq(y,X, Alphas)
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

@rput Y n n_cenarios max_sim
R"
# n = 100
# n_cenarios = 100
# max_sim = 10
 frequencia = 1
# Y = matrix(rnorm((n + max_sim)*n_cenarios),(n + max_sim))

ano_ini = 2006; mes_ini = 1; level = c(95,99.9999); limites.y = NULL; title = NULL
dem = ts(Y[1:n,1], start = c(1,1), frequency = frequencia)
sim = ts(tail(Y, max_sim), start = c(n+1,1), frequency = frequencia)
inicio.serie <- start(dem)
inicio.previsao <- start(sim)

lower <- apply(sim, 1, quantile, 0.5 - level/200, type = 8)
upper <- apply(sim, 1, quantile, 0.5 + level/200, type = 8)
ymax <- max(sim); ymin <- min(sim)
cenario.medio <- rowMeans(sim)
limites.y = c(ymin, ymax)
plot(ts(0, start = inicio.serie), xlim=c(start(dem)[1],end(sim)[1]+1), ylim = limites.y, main = title, ylab = 'MW')
lines(ts(c(dem,cenario.medio[1]), start = inicio.serie, frequency = 1), col = 3)

# dem.ts <- ts(dem, start = inicio.serie, frequency = 12)

lines(ts(cenario.medio, start = inicio.previsao, frequency = frequencia), col = 1);
# lines(ts(lower[1,], start = inicio.previsao, frequency = 12), col = 5, lty = 2);
# lines(ts(upper[1,], start = inicio.previsao, frequency = 12), col = 5, lty = 2);
lines(ts(lower[1,], start = inicio.previsao, frequency = frequencia), col = 5, lty = 2);
lines(ts(upper[1,], start = inicio.previsao, frequency = frequencia), col = 5, lty = 2);
lines(ts(lower[2,], start = inicio.previsao, frequency = frequencia), col = 2);
lines(ts(upper[2,], start = inicio.previsao, frequency = frequencia), col = 2);
# lines(ts(, start = inicio.previsao, frequency = 12), col = 5);
grid()
"
