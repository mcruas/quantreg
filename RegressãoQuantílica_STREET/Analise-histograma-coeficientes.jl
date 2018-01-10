### Este arquivo possui um experimento que compara a estimação por parte da nossa metodologia, em que consideramos uma
# penalização no 

# Loads packages and source files
using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall, Dierckx, JLD
pyplot()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");


# Parameters of simulations
Alphas = collect(0.05:0.05:0.95)
vector_gamma = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0,20.0] 
K_folds = 10
n_iter = 1000
n = 400
rho = 0.3

## Initialize vectors and matrix
keep_betas = zeros(n_iter, length(vector_gamma), length(Alphas)) # values of beta used in each fold on the CV
betas_CV = zeros((n_iter, length(Alphas)))  # values of beta estimated after selecting gamma on the CV
betas_koenker = zeros((n_iter, length(Alphas)))  # values of beta for gamma = 0, as in Koenker (1978)
keep_APD = zeros(n_iter, length(vector_gamma), length(Alphas)) # values of APD for selecting best fold
keep_ar = zeros(n_iter)   # estimated coefficient of autoregressive model
betas_CV = zeros(n_iter, length(Alphas)) # betas that were selected by CV
keep_best_gamma_CV = zeros(n_iter)




name_file = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/variables_$(n)_$(n_iter).jld"
for iter = 1:n_iter  # iter = 1
    tic()
    R"
    #library(forecast)
    set.seed($iter)
    serie <- arima.sim(n = $n, list(ar = c($rho), sd = 1)) # generates time serie
    coef_ar <- coef(arima(serie, order = c(1,0,0)))[1] # gets ar(1) coefficient
    "
    @rget serie coef_ar
    keep_ar[iter] = coef_ar  # keeps value of AR(1) estimation
    X_lags = lagmatrix(serie,0:1)  # Makes a matrix matching y_t and its regressor y_{t-1}
    y = X_lags[:,1]; X = X_lags[:, 2:end]; # Assigns y and X values
    folds = rcopy(R"sample(1:$K_folds, $(length(y)), replace = TRUE)")   # assigns randomly a fold to each 
    for i_gamma =  1:length(vector_gamma)   # For each value of the parameter, performs CV
          # k = 1 ; i_gamma = 3
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
                keep_APD[iter, i_gamma, i_alpha] += abs(sum(q_predicted_fold[:,i_alpha] .> y_test_fold)/length(y_test_fold) - alpha)
                # println("iter $iter    gamma $i_gamma   alpha $i_alpha     APD $(keep_APD[iter, i_gamma, i_alpha])")
            end
            
            

        end # end folds
                
    end # end gamma

    CV_gamma = sum(keep_APD[iter,:,:], 2)  # Finds the CV value for each gamma, by summing every alpha
    index_best_CV = findmin(CV_gamma)[:2]  # Finds the value of index of gamma which minimizes the CV function
    keep_best_gamma_CV[iter] = vector_gamma[index_best_CV]   # Keeps value of best gamma 

    # Reestimates QR with the whole dataset, using the best value of gamma
    results_qr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = vector_gamma[index_best_CV], non_cross = true) 
    beta0 = results_qr[:1]; betas = results_qr[:2]
    betas_CV[iter,:] = betas[1,:]

    # Estimates and registers estimation by koenker (1978)
    results_qr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = 0, non_cross = true) # estimates QR with quantile regularization
    beta0 = results_qr[:1]; betas = results_qr[:2]
    betas_koenker[iter, :] = betas[1,:]

    toc()
    if (iter % 50) == 0
        print(iter)
        save(name_file, "keep_betas", keep_betas ,"betas_CV",  betas_CV  ,"betas_koenker", betas_koenker ,"keep_APD", keep_APD ,"keep_ar",  keep_ar  ,"keep_best_gamma_CV", keep_best_gamma_CV)
    end
        

end # end iter



save(name_file, "keep_betas", keep_betas ,"betas_CV",  betas_CV  ,"betas_koenker", betas_koenker ,"keep_APD", keep_APD ,"keep_ar",  keep_ar  ,"keep_best_gamma_CV", keep_best_gamma_CV)