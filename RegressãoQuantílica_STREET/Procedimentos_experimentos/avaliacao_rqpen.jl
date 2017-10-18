using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx, JLD, RCall #, Distributions


usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
n = 100;

# para icaraizinho
serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"; x_new = 30

n_tmp = length(serie); # keeps size of series before cutting them

### Inicialização dos dados

tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X =zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]
y = serie[3:n_tmp];
n = length(y);55
T = 1:n;


################################################################################################
################################################################################################
# Experimentos com o seleção inteira
X_lags = lagmatrix(serie,0:12)
Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
TimeLimit = 157
vetor_TimeLimit = [200, 600, 1800]


# Recolocar os grupos
vetor_Grupos = [1,2,3,10, length(Alphas)]
vetor_max_K = [2] 
 vetor_TimeLimit = [18000] # Máximo de dois dias de simulação


# TEMP
vetor_Grupos = [10]
vetor_max_K = [2] 
 vetor_TimeLimit = [157] # Máximo de dois dias de simulação



max_K = 2
Grupos = 36


y = X_lags[:,1]; X = X_lags[:, 2:end];   non_cross = true; MIPGap = 0.00


lambda = 0.1
R"

lambda = $lambda
library(rqPen)
# results_r = rq.lasso.fit.mult($X,$y,0.5, 1)
set.seed(123)
X <- matrix(rnorm(8000),ncol=10)
y <- 1 + X[,1] - 3*X[,5] + rnorm(800)
results_r <- rq.lasso.fit.mult(X,y,lambda=lambda)
rq.lasso.fit(X,y,lambda=lambda, tau = 0.7)
"
@rget X y results_r
# gr()
# Guarda_Registros = Registros[];
results = rq_par_lasso(y, X, [0.1; 0.3; 0.5; 0.7; 0.9]; lambda = 0.2, non_cross = false, Save = NaN)

results_r[:coefficients]

results[:1]
results[:2]
results[:3]