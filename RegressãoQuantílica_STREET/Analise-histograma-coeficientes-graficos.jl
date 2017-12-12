
########################## Gráficos ########################
# keep_betas[:,1,:]
n =  400; n_iter = 100
nome_arquivo = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/variables_$(n)_$(n_iter).jld"
keep_betas, keep_APD = load(nome_arquivo, "keep_betas", "keep_apd")

Alphas = collect(0.05:0.05:0.95)
vector_gamma = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0,20.0] 
i_gamma = 9 ;            gamma = vector_gamma[i_gamma]
# # 


R"""
keep_betas_melt = melt($keep_betas, varnames = c("Iteration", "Gamma", "Alpha"), value.name = "beta")
bwplot(beta ~ factor(Alpha) | Gamma, data = keep_betas_melt)
"""

betas_CV = zeros(n_iter, length(Alphas))# betas that were selected by CV
keep_gammas = zeros(n_iter)
for iter = 1:n_iter
    betas_iter = keep_betas[iter,:,:]     # gets the values of
    CV_gamma = sum(keep_APD[iter,:,:], 2)  # Finds the CV for each gamma, by summing every alpha
    index_best_CV = findmin(CV_gamma)[:2]  # Finds the value of index of gamma which minimizes the CV function
    betas_CV = keep_betas[iter,index_best_CV, :]
end

# using StatPlots
# violin(map(x-> string(x), Alphas)',matrix_betas)
# violin(["1" "2"], matrix_betas[:,1:5])

# rcopy(["1" "2"], matrix_betas)
boxplot(matrix_betas)

# blau = rand(100,4) # Four series of 100 points each
# violin(["Series 1" "Series 2" "Series 3" "Series 4"],matrix_betas,leg=false)


######## Investigate iteration / gamma
iter = 10 ;  i_gamma = 1;

R"
set.seed($iter)
serie <- arima.sim(n = $n, list(ar = c(0.6), sd = 1))
"
@rget serie
gamma = vector_gamma[i_gamma]

folds = rcopy(R"sample(1:$K_folds, $(length(y)), replace = TRUE)")

# separates fold
y_fold = y[find(folds .!= k)]
X_fold = X[find(folds .!= k), :]

y_test_fold = y[find(folds .== k)]
X_test_fold = X[find(folds .== k), :]



# Estimate coefficients for fold
results_qr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = gamma, non_cross = true) # estimates QR with quantile regularization
beta0 = results_qr[:1]; betas = results_qr[:2]
plot(Alphas, keep_betas[iter,i_gamma, :], legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "(I)\$\\lambda=$lambda \\quad \\gamma_1=0 \\quad n=$n\$", ylim = ylim)    
# Predicts on the fold
q_predicted_fold = (betas0' * ones(length(X_test_fold))'  + betas' * X_test_fold')'


# Evaluates fit by testing, for each alpha, how close are the number of obs on the expected range
for i_alpha = 1:length(Alphas)
    alpha = Alphas[i_alpha]
    keep_APD[iter, i_gamma, i_alpha] = sum(q_predicted_fold[:,i_alpha] .> y_test_fold)/length(y_test_fold) - alpha
end
            
scatter(X_fold,y_fold, xlab = "\$y_{t-1}\$", ylab =  "\$y_{t}\$", leg = false, grid = false)
plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
s_adap_penlasso = scatter!(X_test_fold , y_test_fold, color = :red)

