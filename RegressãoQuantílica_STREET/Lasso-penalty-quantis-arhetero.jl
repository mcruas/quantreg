    ###########################################################################
    #  Este arquivo produz os gráficos disponíveis na pasta 'Documento Regressao Quantilica/Figuras/Lasso-penalty-quantis-arima'
    #  Ele testa estimação com o LASSO e ADALASSO comparando o resultado "cru" com a metodologia em que fazemos regularização 
    #  nos quantis. 
    #  DADOS SIMULADOS DE UM PROCESSO ARIMA


    nome_arquivo = "Lasso-penalty-quantis-arhetero"
    try mkdir("Documento Regressao Quantilica/Figuras/$nome_arquivo/") end
    using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall, Dierckx #, Distributions

    pyplot()
    usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
    pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
    cd(pasta_trabalho)

    include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
    include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");


    Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
    Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
    Alphas = collect(0.05:0.05:0.95)



for n = [100, 250, 500,1000]  # n = 100; lambda = 3; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        


    X = rand(n)
    y = 0.3 * X .+ (X .+ randn(n))
    X = X[:,:] #makes X a matrix

    # scatter(X,y, xlab = "\$x_t\$", ylab =  "\$y_{t}\$")

        
    nome_pasta = replace("n $n", ".0", "")
    try mkdir("Documento Regressao Quantilica/Figuras/$nome_arquivo/$nome_pasta/") end
    for lambda = [0.1, 0.3, 1.0,3.0,10.0,20.0,30.0,50.0,100.0]   # lambda = 1.0
        for gamma = [0.1, 0.3, 1.0, 3.0, 10.0,20.0]

            # lambda = 5; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        
            # lambda = 3; gamma = 0.5 # TIRAR DEPOIS DOS EXPERIMENTOS        
            ylim = (-1.0, 2.5)
            range_x = collect(minimum(X):0.01:maximum(X))
            

            results_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
            p_lasso = plot(Alphas, results_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(I)\$\\lambda=$lambda \\quad \\gamma_1=0 \\quad n=$n\$", ylim = ylim)
            betas0 = results_lasso[:1] ; betas = results_lasso[:2]
            scatter(X,y, xlab = "\$y_{t-1}\$", ylab =  "\$y_{t}\$", leg = false, grid = false)
            s_lasso = plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
            
            
            
            results_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = gamma, non_cross = true) # estimates lasso with penalization of derivatives
            # q = plot(Alphas, results_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(II)\$\\lambda=$lambda \\quad \\gamma_1 = $gamma\$")
            
            # # w1 = calculate_w_as_norm(results_penlasso[:2])
            w_penlasso = calculate_w_as_weight(results_penlasso[:2])
            w_lasso = calculate_w_as_weight(results_lasso[:2])
            w_lasso_wn = calculate_w_as_weighted_norm(results_lasso[:2])
            w_penlasso_wn = calculate_w_as_weighted_norm(results_penlasso[:2])
            # # results_adap1 = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w1, gamma = gamma, non_cross = true)
            
            
            results_adap_lasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_lasso,  non_cross = true)
            p_adap_lasso = plot(Alphas, results_adap_lasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(I) \$\\gamma_2 = 0\$", ylim = ylim)
            betas0 = results_adap_lasso[:1] ; betas = results_adap_lasso[:2]
            scatter(X,y, xlab = "\$y_{t-1}\$", ylab =  "\$y_{t}\$", leg = false, grid = false)
            s_adap_lasso = plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
            


            results_adap_penlasso = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = gamma, non_cross = true)
            p_adap_penlasso = plot(Alphas, results_adap_penlasso[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = $gamma\$", ylim = ylim)
            betas0 = results_adap_penlasso[:1] ; betas = results_adap_penlasso[:2]
            scatter(X,y, xlab = "\$y_{t-1}\$", ylab =  "\$y_{t}\$", leg = false, grid = false)
            s_adap_penlasso = plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
              



            # results_adap_penlasso_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso, gamma = 0, non_cross = true)
            # t = plot(Alphas, results_adap_penlasso_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive w=(II) \$\\gamma_2 = 0\$")
            
            
            # results_adap_penlasso_wn = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = gamma, non_cross = true)
            # v = plot(Alphas, results_adap_penlasso_wn[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = $gamma\$")
            
            # results_adap_penlasso_wn_notpen = rq_par_lasso(y, X, Alphas; lambda = lambda, w = w_penlasso_wn, gamma = 0, non_cross = true)
            # v1 = plot(Alphas, results_adap_penlasso_wn_notpen[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Adaptive wn=(II) \$\\gamma_2 = 0\$")
            

            results_qr = rq_par(y, X, Alphas; non_cross = true)
            betas0 = results_qr[:1] ; betas = results_qr[:2]
            p_qr = plot(Alphas, betas', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Quantile Regression", ylim = ylim)
            scatter(X,y, xlab = "\$y_{t-1}\$", ylab =  "\$y_{t}\$", leg = false, grid = false)
            s_qr = plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
            
             
    
            # r = plot(Alphas, results_adap1[:2]', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "\$w\$ as norm")
            
            plot(plot(p_qr,p_lasso,s_qr,s_lasso) ,plot(p_adap_lasso,p_adap_penlasso,s_adap_lasso,s_adap_penlasso), size = (1200,700))
            
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

