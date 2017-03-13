# include("par-selecaolasso.jl")

# Verificar:
# - numero de lags
# - tipos de funcoes incluidas
# - os dados gerados
# - o nome dos arquivos
mode_distribution = "par"   # choose between "par" for parametric
                            # and "np" for nonparametric


tic()
n = 1000
using JuMP, DataFrames, Distributions, Plots, RCall, Interpolations, Dierckx #, Distributions
gr()
usesolver = "mosek"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
cd(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/")
# cd("C:/Users/mcruas/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET")
# pwd()
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl")
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl")

################################################################
############# Carregar dados ###################################
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
# @rget serie dados_filtrados
# serie
# plot(serie)

# serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"
serie = readcsv("Dados Climaticos/Enas/enas sudeste 1931-2014.csv"); nomeserie = "enasudeste"
serie = serie ./ 100
max_sim= 100;
n_cenarios = 100;
n_tmp = length(serie); # keeps size of series before cutting them





############################################################
### Inicialização

# Does adaptation for a AR(2) model

max_lag = 1
X = [rand(n_tmp-2) serie[2:n_tmp-1] serie[1:n_tmp-2]] # X = [y_t-1 y_t-2]
X = [serie[2:n_tmp-1] serie[1:n_tmp-2]]

X = zeros(length(1:n_tmp-1),1)
X[:,1] = serie[1:n_tmp-1]

y = serie[max_lag + 1:n_tmp]

n = length(y)
T = 1:n

Alphas = collect(0.05:0.05:0.95)

Y = [ y * ones(1,n_cenarios);
      zeros(max_sim,n_cenarios)]


### Step 1
tau = n +1


### Step 2
X_tau = [y[tau-1]] #não é generico. nao aceitar.

if mode_distribution == "par"
  Q_hat = Estimar_Q_hat_par(y,X,Alphas, X_tau)
else
  Q_hat = Estimar_Q_hat_np(y,X[:,1],Alphas, 0, 300, x_new)
end




###################################################
######################## TMP

# rq_np(y,X[:,1],Alphas, 0, 100)
#
# Q_hat = Estimar_Q_hat_np(y,X[:,1],Alphas, 0, 100, X_tau)
#
# thetas, x_ord, y_ord = rq_np(y,X[:,1],Alphas, 0, 100)
#
#
# Q_hat = Estimar_Q_hat_par(y,X,Alphas, x_new)
# plot(Q_hat)
# Q_hat = Estimar_Q_hat_np(y,X[:,1],Alphas, lambda1, lambda2, x_new,
#                           range_data = range_data,  non_cross = non_cross)
# plot!(Q_hat)
# #
# x = X[:,1]
# non_cross= true
# lambda1 = 0; lambda2 = 300
# range_data = NaN
# degree_splines = 2
# x_new = 500
#
# thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2,
#   range_data = range_data,  non_cross = non_cross)
# Q_hat = zeros(length(Alphas))
# for col_alpha in 1:length(Alphas)
#   sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = degree_splines)
#   Q_hat[col_alpha] = evaluate(sp1, x_new)
# end
#
# tic()
# thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2,
#   range_data = range_data,  non_cross = non_cross)
# scatter(x_ord, y_ord, leg = false)
# plot!(x_ord ,thetas)
# toc()
#
# Q_hat = zeros(length(Alphas))
# for col_alpha in 1:length(Alphas)
#   sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = degree_splines)
#   Q_hat[col_alpha] = evaluate(sp1, x_new)
# end
# return Q_hat
###########################

### Step 3
unif = rand(n_cenarios)
y_simulated = Q(unif, Q_hat, Alphas)
Y[tau, :] = y_simulated

### Step 4
for s in 1:n_cenarios, tau in n+2:n+max_sim

  # (step 2)
  X_tau = [Y[tau-1,s]] #não é generico. nao aceitar.
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
