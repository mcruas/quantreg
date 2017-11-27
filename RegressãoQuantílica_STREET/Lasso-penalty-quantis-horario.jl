###########################################################################
#  Este arquivo produz os gráficos disponíveis na pasta 'Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis-horario'
#  Ele testa estimação com o LASSO e ADALASSO comparando o resultado "cru" com a metodologia em que fazemos regularização 
#  nos quantis. 
#  DADOS HORÁRIOS PROVENIENTES DO KAGGLE 

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
try mkdir("Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis-horario/$nome_pasta/") end
for lambda = [1000.0, 10000.0]   
    for gamma = [0.0]
    
        print(replace("Lambda$lambda-gamma$gamma",".",""))
########## USAR CTRL + SHIFT + Q #############


        # lambda = 0.2; gamma = 1.0 # TIRAR DEPOIS DOS EXPERIMENTOS        
        results_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
        unicodeplots()
        p = plot(Alphas, results_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(I) lambda=$lambda   gamma=$gamma     n=$n")
        pyplot()
        p = plot(Alphas, results_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(I)\$\\lambda=$lambda \\quad \\gamma_1=0 \\quad n=$n\$")
        # results_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = gamma, non_cross = true) # estimates lasso with penalization of derivatives
        # q = plot(Alphas, results_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(II)\$\\lambda=$lambda \\quad \\gamma_1 = $gamma\$")
        
        # # w1 = calculate_w_as_norm(results_penlasso[:2])
        # w_penlasso = calculate_w_as_weight(results_penlasso[:2])
        # w_lasso = calculate_w_as_weight(results_lasso[:2])
        # w_lasso_wn = calculate_w_as_weighted_norm(results_lasso[:2])
        # w_penlasso_wn = calculate_w_as_weighted_norm(results_penlasso[:2])
        # # results_adap1 = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w1, gamma = gamma, non_cross = true)
        # results_adap_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso,  non_cross = true)
        # r = plot(Alphas, results_adap_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(I) \$\\gamma_2 = 0\$")
        
        # results_adap_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = gamma, non_cross = true)
        # s = plot(Alphas, results_adap_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = $gamma\$")
        
        # results_adap_penlasso_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = 0, non_cross = true)
        # t = plot(Alphas, results_adap_penlasso_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = 0\$")
        
        # results_adap_lasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso_wn,  non_cross = true)
        # u = plot(Alphas, results_adap_lasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(I) \$\\gamma_2 = 0\$ ")
        
        # results_adap_penlasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = gamma, non_cross = true)
        # v = plot(Alphas, results_adap_penlasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = $gamma\$")
        
        # results_adap_penlasso_wn_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = 0, non_cross = true)
        # v1 = plot(Alphas, results_adap_penlasso_wn_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = 0\$")
        

        # # r = plot(Alphas, results_adap1[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$w\$ as norm")
        
        # plot(p,q,r,s,t,u,v,v1,size = (1200,700), ylim = (-3.5,1.5))

        name_file = replace("Lambda$lambda-gamma$gamma",".","")
        savefig("Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis-horario/$nome_pasta/$name_file.pdf")
        savefig("Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis-horario/$nome_pasta/$name_file.png")

    end

end

    # lambda = 20; gamma = 0.5
    # # Teste adap e regquant
    # include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
    # results = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
    # results_regquantile = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = gamma, non_cross = true) # estimates lasso with penalization of derivatives
    # results_regquantile[:2]
    # w1 = calculate_w_as_norm(results_regquantile[:2])
    # w2 = calculate_w_as_weight(results_regquantile[:2])
    # results_adap1 = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w1, gamma = gamma, non_cross = true)
    # results_adap2 = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w2, gamma = gamma, non_cross = true)


    # p = plot(Alphas, results[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$\\lambda=$lambda \\quad \\gamma=0\$")
    # s = plot(Alphas, results_regquantile[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$\\gamma = $gamma\$")
    # r = plot(Alphas, results_adap1[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$w\$ as norm")
    # q = plot(Alphas, results_adap2[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$w\$ as weight")

    # plot(p,s,r,q,size = (900,700))


    # # results_r[:coefficients]

    # # results[:1]
    # # results2[:1]
    # # results[:2]
    # # results2[:2]
    # # results[:3]

    # p_tmp = 2; n_tmp = 100000
    # beta_tmp = 1:p_tmp
    # X_tmp = rand(n_tmp,p_tmp)
    # X_scale, X_mean, X_sd = normalization(X_tmp)
    # @rput X_tmp X_scale X_mean X_sd beta_tmp n_tmp
    # R" 

    # X_tmp = cbind(rnorm(n_tmp,1,2), rnorm(n_tmp, 3,4))
    # y_tmp = X_tmp %*% beta_tmp + rnorm(nrow(X_tmp))
    # X_scale = scale(X_tmp)
    # res1 = lm(y_tmp ~ X_tmp)
    # X_mean
    # X_sd
    # res2 = lm(y_tmp ~ X_scale)
    # "







#######################################################################################
# Conclusões parciais do estudo

# - PAra n = 500 (testar outras depois)
# Sem a penalização, os lambdas que fazem os metodos com gamma = 0 selecionar apenas as variaveis
# corretas tambem desconfigura a forma da função beta(alpha), que deveria ser uma reta em -3 e 1

