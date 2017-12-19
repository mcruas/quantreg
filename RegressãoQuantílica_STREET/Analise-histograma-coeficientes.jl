### Este arquivo possui um 


using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall, Dierckx, JLD

pyplot()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");


Alphas = collect(0.05:0.05:0.95)




vector_gamma = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0,20.0] 

K_folds = 10
n_iter = 1003
n = 400

keep_betas = zeros(n_iter, length(vector_gamma), length(Alphas))
keep_APD = zeros(n_iter, length(vector_gamma), length(Alphas))
keep_ar = zeros(n_iter)
name_file = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/variables_$(n)_$(n_iter).jld"

tic()

for iter = 1:n_iter  # iter = 1
    R"
    library(forecast)
    set.seed($iter)
    serie <- arima.sim(n = $n, list(ar = c(0.3), sd = 1)) # generates time serie
    coef_ar <- coef(arima(serie, order = c(1,0,0)))[1] # gets ar(1) coefficient
    "
    @rget serie coef_ar
    
    keep_ar[iter] = coef_ar

    X_lags = lagmatrix(serie,0:1)
    y = X_lags[:,1]; X = X_lags[:, 2:end]; 
    
    range_x = collect(minimum(X):0.01:maximum(X))
    
    folds = rcopy(R"sample(1:$K_folds, $(length(y)), replace = TRUE)")


    
    for i_gamma =  1:length(vector_gamma)   # For each value of the parameter, performs CV
          # k = 1 ; i_gamma = 1
        gamma = vector_gamma[i_gamma]

        #### Inicia aqui o CV
        for k = 1:K_folds # for each fold, estimates and gets the values 
          

            # separates fold
            y_fold = y[find(folds .!= k)]
            X_fold = X[find(folds .!= k), :]

            y_test_fold = y[find(folds .== k)]
            X_test_fold = X[find(folds .== k), :]

            # Estimate coefficients for fold
            results_qr = rq_par_lasso(y_fold, X_fold, Alphas; lambda = 0, gamma = gamma, non_cross = true) # estimates QR with quantile regularization
            beta0 = results_qr[:1]; betas = results_qr[:2]
            keep_betas[iter,i_gamma, :] = betas[1,:]
               


            # Predicts on the fold
            q_predicted_fold = (beta0' * ones(length(X_test_fold))'  + betas' * X_test_fold')'


            # Evaluates fit by testing, for each alpha, how close are the number of obs on the expected range
            for i_alpha = 1:length(Alphas)
                alpha = Alphas[i_alpha]
                keep_APD[iter, i_gamma, i_alpha] = abs(sum(q_predicted_fold[:,i_alpha] .> y_test_fold)/length(y_test_fold) - alpha)
            end
            
            

        end # end folds



    end # end gamma

    if (iter % 50) == 0
        print(iter)
        name_file = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/variables_$(n)_$(n_iter).jld"
        save(name_file, "keep_betas", keep_betas, "keep_APD", keep_APD, "keep_ar", keep_ar)
    end
        

end # end iter


toc()

save(name_file, "keep_betas", keep_betas, "keep_APD", keep_APD, "keep_ar", keep_ar)
