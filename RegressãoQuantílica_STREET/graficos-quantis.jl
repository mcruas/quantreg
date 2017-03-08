# include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/graficos-quantis.jl")

#######################################################
# Este script gera gráficos para o relatório, mostrando como uma variação
# em x_t afeta o quantil do período seguinte.

################# Gráficos linear #####################
vetor_x = [-1 0 1 2 3 10] ; n = 100
limx = 1500; limy = 2000;
using JuMP, DataFrames, Plots, RCall, Interpolations, LaTeXStrings, Dierckx #, Distributions



usesolver = "mosek"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

#####################################################
############# Carregar dados em R ########################
pwd();
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/par-multi.jl");
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


serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"
# serie = readcsv("Dados Climaticos/Enas/enas sudeste 1931-2014.csv"); nomeserie = "enasudeste"
serie = serie ./ 100
max_sim= 10;
n_cenarios = 100;
n_tmp = length(serie); # keeps size of series before cutting them
### Inicialização

# Does adaptation for a AR(2) model

X = [serie[2:n_tmp-1] serie[1:n_tmp-2]]; # AR(2)

# X = Array{Float64,2} # AR(1)
tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X = zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]

y = serie[3:n_tmp];
n = length(y);
T = 1:n;

# Alphas = collect(0.005:0.005:0.995);


Y = [ y * ones(1,n_cenarios);
      zeros(max_sim,n_cenarios)];



#######################################################
############## PARAMETRICO ############################
#######################################################
#######################################################
#######################################################

######### Teste diferentes Alphas
X_tau = 10000
gr()
alpha_plot = collect(0:0.001:1);
Pouco_Alphas = collect(0.05:0.05:0.95);
betas0, betas = rq_linear(y,X, Pouco_Alphas); # coeficientes do modelo linear
Q_hat_pouco_alpha = (betas0 + X_tau' * betas)[1,:];
q_plot = Q(alpha_plot, Q_hat_pouco_alpha, Pouco_Alphas)
plot(alpha_plot, q_plot,ylabel = "Q_{y_{t}|y_{t-1}}", xlabel = "\\alpha",  leg = false, grid = false) # ylim = (0,limy)


Mto_Alphas = collect(0.005:0.005:0.995);
betas0, betas = rq_linear(y,X, Mto_Alphas); # coeficientes do modelo linear
Q_hat_mto_alpha = (betas0 + X_tau' * betas)[1,:];
q_plot = Q(alpha_plot, Q_hat_mto_alpha, Mto_Alphas)
plot!(alpha_plot, q_plot)

# Pouco_Alphas = collect(0.07:0.1:0.97);
# betas0, betas = rq_linear(y,X, Pouco_Alphas); # coeficientes do modelo linear
# Q_hat_pouco_alpha = (betas0 + X_tau' * betas)[1,:];
# q_plot = Q(alpha_plot, Q_hat_pouco_alpha, Pouco_Alphas)
# plot!(alpha_plot, q_plot)
#
# Pouco_Alphas = collect(0.05:0.1:0.95);
# betas0, betas = rq_linear(y,X, Pouco_Alphas); # coeficientes do modelo linear
# Q_hat_pouco_alpha = (betas0 + X_tau' * betas)[1,:];
# q_plot = Q(alpha_plot, Q_hat_pouco_alpha, Pouco_Alphas)
# plot!(alpha_plot, q_plot)
plot!()
savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/quantile-vs-alphas-linear.pdf")


############ Gráficos dos quantis

# x = X[:,1]
# Alphas = collect(0.05:0.05:0.95);
# thetas, x_ord, y_ord = rq_np(y,x,Alphas, 100000000, true)
#
# scatter(x,y, leg = false)
# plot!(sort(x) ,thetas)
#
#








########## Gráficos de Q com a variação do X para modelo linear
gr()

vetor_X_tau = collect(10:5:40)
alpha_plot = collect(0:0.001:1);
Alphas = collect(0.05:0.05:0.95);
betas0, betas = rq_linear(y,X, Alphas, true); # coeficientes do modelo linear
plot(leg = false, ylabel = "Q_{y_{t}|y_{t-1}}", xlabel = "\\alpha", ylim = (0,limy), grid = false)
for i in 1:length(vetor_X_tau)
  X_i = [vetor_X_tau[i]]
  Q_hat = (betas0 + X_i' * betas)[1,:];
  q_plot = Q(alpha_plot, Q_hat, Alphas)
  plot!(alpha_plot, q_plot, label = "X = $X_i")
end
plot!()

savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/quantile-linear.pdf")


############ Scatterplot com os quantis com modelo paramétrico
scatter(X, y, leg = false, xlabel = "y_{t-1}", ylabel = "y_t", xlim = (0,limx), ylim= (0,limy), grid = false)
range_x = collect(minimum(X):0.01:maximum(X))
plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/quantile-linear-scatter.pdf")

# X_i = [vetor_X_tau[1]]




#######################################################
############## NÃO PARAMETRICO ########################
#######################################################
#######################################################
#######################################################

########## Gráficos de Q com a variação de X para modelo não-paramétrico

Alphas = collect(0.05:0.05:0.95);
thetas, x_ord, y_ord = rq_np(y,x,Alphas,1, 10, true)


# Finds Q_hat given a new value of x
x_new = 30
Q_hat = zeros(length(Alphas))
plot(leg = false, ylabel = "Q_{y_{t}|y_{t-1}}", xlabel = "\\alpha", title = "", ylim = (0,limy), grid=false)
for x_new in vetor_X_tau
  for col_alpha in 1:length(Alphas)
    sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = 2)
    Q_hat[col_alpha] = evaluate(sp1, x_new)
  end
  # scatter!(repmat([x_new], length(Alphas)), Q_hat)

  # Finds the whole distribution for a sequence of values
  q_plot = Q(alpha_plot, Q_hat, Alphas)
  plot!(alpha_plot, q_plot)
end
plot!()

savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/quantile-nonpar.pdf")



scatter(x,y, leg = false, xlabel = "y_{t-1}", ylabel = "y_t", xlims = (0,limx), ylim= (0,limy), grid = false)
plot!(sort(x) ,thetas)
savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/quantile-nonpar-scatter.pdf")


######################## Formato de Q variando granularidade de A
x_new = 10000
Q_hat = zeros(length(Alphas))

plot(leg = false, ylabel = "Q_{y_{t}|y_{t-1}}", xlabel = "\\alpha", title = "", ylim = (0,limy), grid=false)

Alphas = collect(0.05:0.05:0.95);
thetas, x_ord, y_ord = rq_np(y,x,Alphas, 0, 10, true)
for col_alpha in 1:length(Alphas)
  sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = 2)
  Q_hat[col_alpha] = evaluate(sp1, x_new)
end
# Finds the whole distribution for a sequence of values
q_plot = Q(alpha_plot, Q_hat, Alphas)
plot!(alpha_plot, q_plot)

Alphas = collect(0.005:0.005:0.995);
# thetas, x_ord, y_ord = rq_np(y,x,Alphas, 0, 10, non_cross=true)
# Q_hat = zeros(length(Alphas))
# for col_alpha in 1:length(Alphas)
#   sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = 2)
#   Q_hat[col_alpha] = evaluate(sp1, x_new)
# end
Q_hat  = Estimar_Q_hat_np(y,x,Alphas, 0,10,x_new)
q_plot = Q(alpha_plot, Q_hat, Alphas)
plot!(alpha_plot, q_plot)


# Finds the whole distribution for a sequence of values
q_plot = Q(alpha_plot, Q_hat, Alphas)
plot!(alpha_plot, q_plot)

savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/quantile-vs-alphas-nonpar.pdf")



plot!()

#############################################################################
######################## Gráficos do quantis (scatterplot) não-paramétrico variando os lambdas

Alphas = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99];
# lambdas1 = [0.1, 0.3, 1, 3, 10, 200]
lambdas1 = [0]

# Lambdas2 para utilizar com icaraizinho
# lambdas2 = [0.1, 0.3, 1, 3, 10, 200]
# nomes = ["01", "03", "1", "3", "10", "200"]

# lambdas2 para utilizar com a Ena sudeste
lambdas2 = [0.1, 0.3, 1, 3, 10, 200]
nomes = ["01", "03", "1", "3", "10", "200"]


for i in 1:length(lambdas1), j in 1:length(lambdas2)
  lambda1 = lambdas1[i]
  lambda2 = lambdas2[j]
  thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2, true)
  scatter(x,y, leg = false, xlabel = "y_{t-1}", ylabel = "y_t", xlim = (0,limx), ylim= (0,limy), grid = false)
  plot!(sort(x) ,thetas)
  # savefig("Documento Regressao Quantilica/Figuras/regressao-quantilica/$nomeserie-crossing-" * nomes[j] * ".pdf")
end



pyplot()
lambda1= 0
lambda2 = 10000
thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2, range_data = [0, 15000])
scatter(x_ord,y_ord, leg = false, xlabel = "y_{t-1}", ylabel = "y_t", grid = false,  xlim = (0,limx), ylim= (0,limy))
plot!(sort(x_ord) ,thetas)


if any(!isnan(range_data))
  print("hey")
end
