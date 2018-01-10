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


# Define os dados como serão
seq_X = 0:0.01:5
X = seq_X[:,:]


# Parameters of simulations
n_methods = 3 # quantity of methods we use; (qr , ar and ours)
Alphas = collect(0.05:0.05:0.95)
vector_gamma = [0.0,0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0,30.0] 
# vector_gamma = [0.0,0.001, 0.01, 0.1, 1] 
K_folds = 10
n_iter = 10
n = 400
n_test = 400 

# β = 0.3 - 0.6 * (Alphas .- (Alphas .^ 2)) 
β = 0.3 + 0.6 * (Alphas .- (Alphas .^ 2))  # plot(β')
β = reshape(β, 1, length(β))
beta0_sim,betas_sim = SimB0Qar(X, Alphas, β)
plot(seq_X, (beta0_sim' * ones(length(seq_X))'  + betas_sim' * seq_X')', leg = false) # plot qar process

## Initialize vectors and matrix
keep_betas = zeros(n_iter, length(vector_gamma), length(Alphas)) # values of beta used in each fold on the CV
betas_CV = zeros((n_iter, length(Alphas)))  # values of beta estimated after selecting gamma on the CV
betas_koenker = zeros((n_iter, length(Alphas)))  # values of beta for gamma = 0, as in Koenker (1978)
keep_APD = zeros(length(vector_gamma), length(Alphas)) # values of APD for selecting best fold
keep_score_APD = zeros(n_iter, n_methods) # values of APD for each method and iteration
keep_ar = zeros(n_iter)   # estimated coefficient of autoregressive model
betas_CV = zeros(n_iter, length(Alphas)) # betas that were selected by CV
keep_best_gamma_CV = zeros(n_iter)




name_file = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients-qar/variables_$(n)_$(n_iter).jld"
for iter = 1:n_iter  # iter = 1
    tic()   
    ##### generates series
    serie_total = simqar(Alphas, n + n_test + 1, iter, beta0_sim, betas_sim; y0 = mean(seq_X)) # simulates series for the training and test data
    X_lags_total = lagmatrix(serie_total,0:1)  # Makes a matrix matching y_t and its regressor y_{t-1}
    serie = serie_total[1:n]
    y = X_lags_total[1:n,1]; X = X_lags_total[1:n, 2:end]; # Assigns y and X values
    y_test = X_lags_total[n+1:end,1]; X_test = X_lags_total[n+1:end, 2:end]; # Assigns y and X values to the test data



    ##### The following code does cross validation selection to find the best value of
    # gamma for prediction with our method
    folds = rcopy(R"sample(1:$K_folds, $(length(y)), replace = TRUE)")   # assigns randomly a fold to each 
    keep_APD = zeros(length(vector_gamma), length(Alphas)) # resets value of APD for selecting best fold
    
    ## Begining of Cross-Validation
    for i_gamma =  1:length(vector_gamma)   # For each value of the parameter, performs CV
          # k = 2 ; i_gamma = 1
        gamma = vector_gamma[i_gamma]

        for k = 1:K_folds # for each fold, estimates and gets the values 
          
            # separates fold
            y_fold = y[find(folds .!= k)]
            X_fold = X[find(folds .!= k), :]
            y_test_fold = y[find(folds .== k)]
            X_test_fold = X[find(folds .== k), :]

            # Estimate coefficients for fold
            results_qr = rq_par_lasso(y_fold, X_fold, Alphas; lambda = 0, gamma = gamma, non_cross = true) # estimates QR with quantile regularization
            beta0 = results_qr[:1]; betas = results_qr[:2]
            keep_betas[iter,i_gamma, :] = betas[:]
               
            # Predicts on the fold
            q_predicted_fold = (beta0' * ones(length(X_test_fold))'  + betas' * X_test_fold')'

            # Evaluates fit by testing, for each alpha, how close are the number of obs on the expected range
            keep_APD[i_gamma, :] += Calculate_APD(y_test_fold, q_predicted_fold, Alphas)


        end # end folds
                
    end # end gamma

    CV_gamma = sum(keep_APD, 2)  # Finds the CV value for each gamma, by summing every alpha
    index_best_CV = findmin(CV_gamma)[:2]  # Finds the value of index of gamma which minimizes the CV function
    keep_best_gamma_CV[iter] = vector_gamma[index_best_CV]   # Keeps value of best gamma 

    # Reestimates QR with the whole dataset, using the best value of gamma
    results_qr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = vector_gamma[index_best_CV], non_cross = true) 
    beta0 = results_qr[:1]; betas = results_qr[:2]
    betas_CV[iter,:] = betas[:]

    # With the best coefficients, makes prediction for the test set
    q_predicted_test = (beta0' * ones(length(X_test))'  + betas' * X_test')'
    keep_score_APD[iter, 1] = mean(Calculate_APD(y_test, q_predicted_test, Alphas))  # keeps the APD for this method  

    ###### Estimates and registers estimation by koenker (1978)
    results_qr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = 0, non_cross = true) # estimates QR with quantile regularization
    beta0 = results_qr[:1]; betas = results_qr[:2]
    q_predicted_test = (beta0' * ones(length(X_test))'  + betas' * X_test')'
    betas_koenker[iter, :] = betas[:]
    keep_score_APD[iter, 2] = mean(Calculate_APD(y_test, q_predicted_test, Alphas))  # keeps the APD for this method  


    ###### Estimates AR coefficient
    R"coef_ar <- coef(lm($y ~ $X))"
    @rget coef_ar
    keep_ar[iter] = coef_ar[2] 
    y_predicted_test = coef_ar[2,:] .* ones(length(X_test)) .+ coef_ar[1] * X_test # Calculates E(y_t|y_t-1)
    quantiles_normal = rcopy(R"qnorm($Alphas)")
    q_predicted_test = y_predicted_test[:,:] * ones(1,length(Alphas))  .+ ones(length(X_test)) * quantiles_normal[:,:]'
    keep_score_APD[iter, 3] = mean(Calculate_APD(y_test, q_predicted_test, Alphas))  # keeps the APD for this method  
        

    toc()
    if (iter % 50) == 0
        print(iter)
        save(name_file, "keep_betas", keep_betas ,"betas_CV",  betas_CV  ,
            "betas_koenker", betas_koenker ,"keep_APD", keep_APD ,"keep_ar",  
            keep_ar  ,"keep_best_gamma_CV", keep_best_gamma_CV, "keep_score_APD", 
            keep_score_APD, "beta0_sim", beta0_sim, "betas_sim" , betas_sim)
    end

end # end iter

save(name_file, "keep_betas", keep_betas ,"betas_CV",  betas_CV  ,
    "betas_koenker", betas_koenker ,"keep_APD", keep_APD ,"keep_ar",  
    keep_ar  ,"keep_best_gamma_CV", keep_best_gamma_CV, "keep_score_APD", 
    keep_score_APD, "beta0_sim", beta0_sim, "betas_sim" , betas_sim)