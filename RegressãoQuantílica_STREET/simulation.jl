# include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/simulation.jl")

# Verificar:
# - numero de lags
# - tipos de funcoes incluidas
# - os dados gerados
# - o nome dos arquivos
mode_distribution = "np"   # choose between "par" for parametric
                            # and "np" for nonparametric

lambda1 = 0; lambda2 = 300
non_cross = true

# serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"

max_sim= 10;
n_cenarios = 20;



tic()
using JuMP, DataFrames, Distributions, Plots, RCall, Interpolations, Dierckx #, Distributions
gr()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
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
serie = readcsv("Dados Climaticos/Enas/enas sudeste 1931-2014.csv"); nomeserie = "enasudeste"
serie = serie ./ 100

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


#############



#############




### Step 1
tau = n +1


### Step 2
X_tau = [y[tau-1]] #não é generico. nao aceitar.

if mode_distribution == "par"
  betas0, betas = rq_par(y,X,Alphas, non_cross = non_cross)
  Q_hat = Estimar_Q_hat_par2(betas0, betas, X_tau)
else
  thetas, x_ord, y_ord = rq_np(y,X[:,1],Alphas, lambda1, lambda2, non_cross = non_cross)
  Q_hat  = Estimar_Q_hat_np2(x_ord,thetas,Alphas, X_tau[1])
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

  if mode_distribution == "par"
    Q_hat = Estimar_Q_hat_par2(betas0, betas, X_tau)
  else
    Q_hat  = Estimar_Q_hat_np2(x_ord,thetas,Alphas, X_tau[1])
  end


  # (step 4)
  unif = rand(1)
  y_simulated = Q(unif, Q_hat, Alphas)[1,1]
  Y[tau, s] = y_simulated


end



plot(Y[1:n,1], legend = false, title = "Scenarios $mode_distribution")
media = mapslices(mean, Y[n+1:end, :], 2)
plot!(n+1:n+max_sim, media)
p_quantis = [0.75,0.90, 0.95];
p_quantis = vcat(reverse(1 .- p_quantis), p_quantis)'
quantis = zeros(max_sim, length(p_quantis))
for i in 1:max_sim
  quantis[i,:] = quantile(Y[n+i,:], p_quantis)
end
# plot!(n+1:n+max_sim, quantis[:, [1;4]])
for i in 1:convert(Int64,length(p_quantis)/2)
  plot!(n+1:n+max_sim, quantis[:,[i,length(p_quantis)-i+1]], color = i+2)
end
plot!()

savefig("Documento Regressao Quantilica/Figuras/simulation/$nomeserie-$mode_distribution.pdf")

#
#
# @rput Y n n_cenarios max_sim
# R"
# # n = 100
# # n_cenarios = 100
# # max_sim = 10
#  frequencia = 1
# # Y = matrix(rnorm((n + max_sim)*n_cenarios),(n + max_sim))
#
# ano_ini = 2006; mes_ini = 1; level = c(95,99.9999); limites.y = NULL; title = NULL
# dem = ts(Y[1:n,1], start = c(1,1), frequency = frequencia)
# sim = ts(tail(Y, max_sim), start = c(n+1,1), frequency = frequencia)
# inicio.serie <- start(dem)
# inicio.previsao <- start(sim)
#
# lower <- apply(sim, 1, quantile, 0.5 - level/200, type = 8)
# upper <- apply(sim, 1, quantile, 0.5 + level/200, type = 8)
# ymax <- max(sim); ymin <- min(sim)
# cenario.medio <- rowMeans(sim)
# limites.y = c(ymin, ymax)
# plot(ts(0, start = inicio.serie), xlim=c(start(dem)[1],end(sim)[1]+1), ylim = limites.y, main = title, ylab = 'MW')
# lines(ts(c(dem,cenario.medio[1]), start = inicio.serie, frequency = 1), col = 3)
#
# # dem.ts <- ts(dem, start = inicio.serie, frequency = 12)
#
# lines(ts(cenario.medio, start = inicio.previsao, frequency = frequencia), col = 1);
# # lines(ts(lower[1,], start = inicio.previsao, frequency = 12), col = 5, lty = 2);
# # lines(ts(upper[1,], start = inicio.previsao, frequency = 12), col = 5, lty = 2);
# lines(ts(lower[1,], start = inicio.previsao, frequency = frequencia), col = 5, lty = 2);
# lines(ts(upper[1,], start = inicio.previsao, frequency = frequencia), col = 5, lty = 2);
# lines(ts(lower[2,], start = inicio.previsao, frequency = frequencia), col = 2);
# lines(ts(upper[2,], start = inicio.previsao, frequency = frequencia), col = 2);
# # lines(ts(, start = inicio.previsao, frequency = 12), col = 5);
# grid()
# "
# using RCall
# R"
# library(stringr)
# pasta.origem <- 'RegressãoQuantílica_STREET/'
# # pasta.destino <- 'Documento Regressao Quantilica/Figuras/selecao-lasso/'
# nomes <- list.files(path = pasta.origem, pattern = 'table-betas-postlasso-')
# P = 12
# epsilon = 0.000000000000001
# maxK = P
# arquivos_tabelas_lasso <- str_c(pasta.origem,nomes)
#
# nomes <- list.files(path = pasta.origem, pattern = 'table-betas-selecaointeira-')
# arquivos_tabelas_mip <- str_c(pasta.origem,nomes)
# "
# @rget arquivos_tabelas_lasso arquivos_tabelas_mip
#
#
#
# macro rget(args...)
#     blk = Expr(:block)
#     for a in args
#         if isa(a,Symbol)
#             v = a
#             push!(blk.args,:($(esc(v)) = rcopy(Const.GlobalEnv[$(QuoteNode(v))])))
#         elseif isa(a,Expr) && a.head == :(::)
#             v = a.args[1]
#             T = a.args[2]
#             push!(blk.args,:($(esc(v)) = rcopy($(esc(T)),Const.GlobalEnv[$(QuoteNode(v))])))
#         else
#             error("Incorrect usage of @rget")
#         end
#     end
#     blk
# end
#
# :($(esc(x))
#
#
# macro blau(args...)
#     blk = Expr(:block)
#     for a in args
#         if isa(a,Symbol)
#           print("$a is a symbol!")
#         else
#           print("$a is NOT a symbol!")
#
#         end
#     end
#     blk
# end




tmp = Symbol("x")
@blau tmp
#
#
# typeof(Symbol("x"))
# @rget Symbol("x")
#
#
# reval("t")
# rcopy(parse(t[1]))
# t[1]
# @show symbol(t[1])
