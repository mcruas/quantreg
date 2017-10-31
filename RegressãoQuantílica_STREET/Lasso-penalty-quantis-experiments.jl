

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
Alphas = [0.1,0.25,0.5,0.75,0.9]; lambda = 0.1; n = 100
@rput n
R"
lambda = $lambda
library(rqPen)
# results_r = rq.lasso.fit.mult($X,$y,0.5, 1)
set.seed(123)
X <- cbind(rnorm(n,2,1),runif(n,1,2),rnorm(n,0,1.5))
y <- 1 + X[,1] - 3*X[,2] + rnorm(n)
"
@rget X y 


function denormalize(betas_til, beta0_til, X_mean, X_sd)
    # programar depois
    beta[p] = beta_til[p]/s[p]
    beta0_ = beta0_til - sum(x_mean[p],beta[p] for p = P)
end



function normalization(X)
    R"
    X_scale = scale($X)
    X_mean = colMeans($X)
    X_sd = apply($X,2,sd)
    "
    @rget X_scale X_mean X_sd
    return X_scale, X_mean, X_sd

  end

lambda = 0
tmp_lasso1 = rq_par_lasso(y, X, Alphas; lambda = lambda, gamma = 0, non_cross = true) # estimates the normal lasso
betas = tmp_lasso1[:2]
beta0 = tmp_lasso1[:1]
X_norm , X_mean, X_sd = normalization(X)

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