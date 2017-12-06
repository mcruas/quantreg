

using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall #, Distributions


usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");

# para icaraizinho

serie_tmp = readtable("Dados Climaticos/Dados-kaggle/WindData_melt.csv"); nomeserie = "W1"
serie_tmp = serie_tmp[1:50000,:]

serie = serie_tmp[serie_tmp[:location] .== "W1", :wind]
n_tmp = length(serie); # keeps size of series before cutting them

### Inicialização dos dados

tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X =zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]
y = serie[3:n_tmp];
n = length(y);
T = 1:n;

################################################################################################
################################################################################################
# Experimentos com o seleção inteira
X_lags = lagmatrix(serie,0:12)
Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
Alphas= collect(0.05:0.05:0.95)
TimeLimit = 157
vetor_TimeLimit = [200, 600, 1800]


vetor_TimeLimit = [18000] # Máximo de dois dias de simulação


# TEMP
vetor_Grupos = [10]
vetor_max_K = [2] 
vetor_TimeLimit = [157] # Máximo de dois dias de simulação




y = X_lags[:,1]; X = X_lags[:, 2:end];   non_cross = true; MIPGap = 0.00

pyplot()
lambda = 0.00000 ; gamma = 0.0; unicodeplots()
# nome_pasta = replace("n $n", ".0", "")
nome_pasta = "50000"
        # lambda = 5; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        
        results_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = false) # estimates the normal lasso
        p = plot(Alphas, results_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(I)\$\\lambda=$lambda \\quad \\gamma_1=0\$")
        results_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = gamma, non_cross = false) # estimates lasso with penalization of derivatives
        q = plot(Alphas, results_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(II)\$\\gamma_1 = $gamma\$")
        
        # w1 = calculate_w_as_norm(results_penlasso[:2])
        w_penlasso = calculate_w_as_weight(results_penlasso[:2])
        w_lasso = calculate_w_as_weight(results_lasso[:2])
        w_lasso_wn = calculate_w_as_weighted_norm(results_lasso[:2])
        w_penlasso_wn = calculate_w_as_weighted_norm(results_penlasso[:2])
        # results_adap1 = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w1, gamma = gamma, non_cross = false)
        results_adap_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso,  non_cross = false)
        r = plot(Alphas, results_adap_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(I) \$\\gamma_2 = 0\$")
        
        results_adap_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = gamma, non_cross = false)
        s = plot(Alphas, results_adap_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = $gamma\$")
        
        results_adap_penlasso_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = 0, non_cross = false)
        t = plot(Alphas, results_adap_penlasso_wn_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = 0\$")
        
        results_adap_lasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso_wn,  non_cross = false)
        u = plot(Alphas, results_adap_lasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(I) \$\\gamma_2 = 0\$ ")
        
        results_adap_penlasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = gamma, non_cross = false)
        v = plot(Alphas, results_adap_penlasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = $gamma\$")
        
        results_adap_penlasso_wn_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = 0, non_cross = false)
        v1 = plot(Alphas, results_adap_penlasso_wn_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = 0\$")
        

        # r = plot(Alphas, results_adap1[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$w\$ as norm")
        
        plot(p,q,r,s,t,u,v,v1,size = (1200,700))

        savefig("Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis/Experiments/.pdf")






################################### TESTA NORMALIZAÇÂO #########################
Alphas = [0.1,0.25,0.5,0.75,0.9]; lambda = 0.1; n = 100000
@rput n
R"
lambda = $lambda
library(rqPen)
# results_r = rq.lasso.fit.mult($X,$y,0.5, 1)
set.seed(123)
X <- cbind(rnorm(n,2,1),runif(n,1,2),rnorm(n,0,1.5))
y <- 2 + 3*X[,1] - 2*X[,2] + rnorm(n)
lm(y~X)
"
@rget X y 






lambda = 0 
tmp_lasso1 = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
tmp_lasso2 = rq_par_lasso_


tmp_lasso2 = rq_par_lasso(y, X_norm, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
betas_til = tmp_lasso2[:2]
beta0_til = tmp_lasso2[:1]

betas_new = betas * 0;
for p = 1:size(betas)[1], j = 1:size(betas)[2]
    betas_new[p,j] = betas_til[p,j] / X_sd[p]
end
beta0_new = beta0 * 0;
for j = 1:length(beta0)

betas_new - betas



