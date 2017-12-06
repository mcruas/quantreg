###########################################################################
#  Este arquivo produz os gráficos disponíveis na pasta 'Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis-arima'
#  Ele testa estimação com o LASSO e ADALASSO comparando o resultado "cru" com a metodologia em que fazemos regularização 
#  nos quantis. 
#  DADOS SIMULADOS DE UM PROCESSO ARIMA


using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall, Dierckx #, Distributions


# Function to generate values that follow a QAR process
function simqar(Alphas, n, seed) # n = 100 ; seed = 123
    n_ini = 1000 # number of initial observations for 
    n_burnin = 3000
    
    # At first, puts some initial value for the time series y
    R"
    set.seed($seed)
    y = arima.sim(n = $n_ini, list(ar = c(0.75, -0.5, rep(0,9), 0.15)), sd = 1)    
    #y = rep(0,n)
    "
    @rget y
    y = [y;zeros(n + n_burnin)]
    

    # create quantile function
    Alphas = collect(0.05:0.05:0.95)
    Q_hat = zeros(length(Alphas))
    beta1 =  0.3 * ones(length(Alphas))
    beta12 = (0.7 - 0.3 * Alphas)
    beta2 = +0.6 * (Alphas .- (Alphas .^ 2)) .- 0.3
    for t = n_ini+1:length(y)
        Q_hat = 0.9 .+ beta1 * y[t-1] .+ beta12 * y[t-12] .+ beta2 * y[t-2]
        unif = rand(1)
        y[t] = (Q(unif, Q_hat, Alphas) + 0.5 * randn(1))[1]
    end 
    return y[end-n:end]
end   



nome_arquivo = "Lasso-penalty-quantis-qar"


pyplot()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");


################################################################################################
################################################################################################
# Experimentos com o seleção inteira

Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
Alphas= collect(0.05:0.05:0.95)



for n = [100, 250, 500,1000]  # n = 5000; lambda = 3; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        
    serie = simqar(Alphas, n, n)
    X_lags = lagmatrix(serie,0:12)
    y = X_lags[:,1]; X = X_lags[:, 2:end]; 
    
    # gr()
    # Guarda_Registros = Registros[];
    
    nome_pasta = replace("n $n", ".0", "")
    try mkdir("Documento Regressao Quantilica/Figuras/$nome_arquivo/$nome_pasta/") end
    for lambda = [0.01, 0.03, 0.1, 0.3, 1.0,3.0,10.0,20.0,30.0,50.0,100.0]   # lambda = 1.0
        for gamma = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0]

            # lambda = 5; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        
            # lambda = 3; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        
            results_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
            p = plot(Alphas, results_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(I)\$\\lambda=$lambda \\quad \\gamma_1=0 \\quad n=$n\$")
            results_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = gamma, non_cross = true) # estimates lasso with penalization of derivatives
            q = plot(Alphas, results_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(II)\$\\lambda=$lambda \\quad \\gamma_1 = $gamma\$")
            
            # w1 = calculate_w_as_norm(results_penlasso[:2])
            w_penlasso = calculate_w_as_weight(results_penlasso[:2])
            w_lasso = calculate_w_as_weight(results_lasso[:2])
            w_lasso_wn = calculate_w_as_weighted_norm(results_lasso[:2])
            w_penlasso_wn = calculate_w_as_weighted_norm(results_penlasso[:2])
            # results_adap1 = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w1, gamma = gamma, non_cross = true)
            results_adap_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso,  non_cross = true)
            r = plot(Alphas, results_adap_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(I) \$\\gamma_2 = 0\$")
            
            results_adap_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = gamma, non_cross = true)
            s = plot(Alphas, results_adap_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = $gamma\$")
            
            results_adap_penlasso_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = 0, non_cross = true)
            t = plot(Alphas, results_adap_penlasso_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = 0\$")
            
            results_adap_lasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso_wn,  non_cross = true)
            u = plot(Alphas, results_adap_lasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(I) \$\\gamma_2 = 0\$ ")
            
            results_adap_penlasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = gamma, non_cross = true)
            v = plot(Alphas, results_adap_penlasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = $gamma\$")
            
            results_adap_penlasso_wn_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = 0, non_cross = true)
            v1 = plot(Alphas, results_adap_penlasso_wn_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = 0\$")
            
            results_qr = rq_par(y, X, Alphas; non_cross = true)
            v2 = plot(Alphas, results_qr[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Quantile Regression")

            # r = plot(Alphas, results_adap1[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$w\$ as norm")
            
            plot(p,r,u, q ,s,t,v,v1,v2, size = (1200,700), ylim = (-3.5,1.5))
    
            name_file = replace("Lambda$lambda-gamma$gamma",".","")
            savefig("Documento Regressao Quantilica/Figuras/$nome_arquivo/$nome_pasta/$name_file.pdf")
            savefig("Documento Regressao Quantilica/Figuras/$nome_arquivo/$nome_pasta/$name_file.png")

        end

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

