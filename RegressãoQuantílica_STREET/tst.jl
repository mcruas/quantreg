
#######################################################
# Este script gera gráficos para o relatório, mostrando como uma variação
# em x_t afeta o quantil do período seguinte.

################# Gráficos linear #####################
vetor_x = [-1 0 1 2 3 10] ; n = 100
# X_tau = 500 # Qual o valor de X utilizado para a previsão
# limx = 1500; limy = 2000; # para ENA sudeste
limx = 55; limy = 55; # para icaraizinho
# limx = NaN; limy = NaN # quando não se sabe

using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx, JLD, RCall #, Distributions



usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
# max_K = 4
cd(pasta_trabalho)

#####################################################
############# Carregar dados em R ########################
pwd();
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
# include(pwd()*"/RegressãoQuantílica_STREET/par-multi.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
n = 100;

serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"; x_new = 30

n_tmp = length(serie); # keeps size of series before cutting them
### Inicialização

tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X =zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]

y = serie[3:n_tmp];
n = length(y);55
T = 1:n;

# X_lags = lagmatrix(serie,0:2)
# Alphas = collect(0.005:0.005:0.995);


# Y = [ y * ones(1,n_cenarios);
#       zeros(max_sim,n_cenarios)];

limx = isnan(limx) ? maximum(serie) : limx
limy = isnan(limy) ? maximum(serie) : limy

# plot(serie)


################################################################################################
################################################################################################
# Experimentos com o seleção inteira
X_lags = lagmatrix(serie,0:12)
Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
TimeLimit = 200
vetor_TimeLimit = [200, 600, 1800]
vetor_Grupos = [1,2,3,10]

vetor_max_K = [1,5,6] 
# vetor_TimeLimit = [200, 600, 1800, 5400]


# vetor_Grupos = [1,2,3,10]
max_K = 5
Grupos = 2
# gr()
# Guarda_Registros = Registros[];
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");

 y=X_lags[:,1]; X=X_lags[:, 2:end]; non_cross = true; max_K = max_K; TimeLimit = 200; MIPGap = 0.00; Grupos = Grupos


vetor_lambdas = 0:100:10000
keep_betas = Array{Float64,2}[]
# initialize Arrays
for alf in 1:length(Alphas) 
    push!(keep_betas, zeros(12,length(vetor_lambdas)))
end

for i_lambda in 1:length(vetor_lambdas)
        lambda = vetor_lambdas[i_lambda]
        betas0opt, betasopt, objectiveValue, status, solvetime = rq_par_lasso(X_lags[:,1], X_lags[:, 2:end], Alphas;
              lambda = lambda,  non_cross = true)       
        for alf in 1:length(Alphas) 
            keep_betas[alf][:,i_lambda] = betasopt[:,alf]
        end

end

unicodeplots()

# betas0opt + X * betasopt  
for alf in 1:length(Alphas) 
            println("Probabilidade $(Alphas[alf])")
            display(plot(keep_betas[alf]'))
end
unicodeplots()
plot(vetor_lambdas,  keep_betas[18]')

R"x = rnorm(1999)
  hist(x)"