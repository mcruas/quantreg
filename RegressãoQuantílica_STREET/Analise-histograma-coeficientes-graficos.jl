using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall, Dierckx, JLD

########################## Gráficos ########################
# keep_betas[:,1,:]
n =  400; n_iter = 1000
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)
nome_arquivo = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/variables_$(n)_$(n_iter).jld"
pasta_raiz = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/"
keep_betas, keep_APD = load(nome_arquivo, "keep_betas", "keep_APD")
keep_betas, keep_APD, keep_ar = load(nome_arquivo, "keep_betas", "keep_APD", "keep_ar")
betas_koenker = keep_betas[:,1,:] # Takes coefficients with gamma = 0 to compare


Alphas = collect(0.05:0.05:0.95)
vector_gamma = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0,20.0] 


### Executes this if using a partial version of JLD files
n_iter = find(keep_betas[:,1,1] .== 0)[1] - 1; 
keep_betas = keep_betas[1:n_iter, : , :]
keep_APD = keep_APD[1:n_iter, :_APD, keep_ar = , :]
betas_koenker = keep_betas[:,1,:]


## Seleciona apenas os betas por CV para cada iteração
betas_CV = zeros(n_iter, length(Alphas)) # betas that were selected by CV
keep_gammas = zeros(n_iter)
indexes_best_CV = zeros(n_iter)
for iter = 1:n_iter
    betas_iter = keep_betas[iter,:,:]     # gets the values of
    CV_gamma = sum(keep_APD[iter,:,:], 2)  # Finds the CV for each gamma, by summing every alpha
    index_best_CV = findmin(CV_gamma)[:2]  # Finds the value of index of gamma which minimizes the CV function
    indexes_best_CV[iter] = index_best_CV
    betas_CV[iter,:] = keep_betas[iter,index_best_CV, :]
end

## Shows the frequency of the gamma selected as the best ones on the CV
best_gammas = map(x->(vector_gamma[convert(Int64,x)]), indexes_best_CV)
R"barplot(table($best_gammas), main = 'Frequency of selected gammas')"

## Gammas Boxplots 
R"""
library(reshape2); library(lattice)
keep_betas_melt = melt($keep_betas, varnames = c("Iteration", "Gamma", "Alpha"), value.name = "beta")
bwplot(beta ~ factor(Alpha, labels = $Alphas) | factor(Gamma, labels = $vector_gamma), data = keep_betas_melt, xlab = $Alphas)
"""
# Análise:
# De modo geral, o erro dos coeficientes é aproximadamente o mesmo para todos 
# os valores de gamma

## Teste de hipóteses de diferença de medias




## Teste de hipóteses de variância
pvalues = zeros(length(Alphas))
for i_alpha = 1:length(Alphas)
    var_CV = var(betas_CV[:,i_alpha])
    var_0 = var(keep_betas[:,1,i_alpha])
    pvalues[i_alpha] = rcopy(R"pf($var_0/$var_CV,$n_iter-1,$n_iter-1)")
end

nome_tex = "$(pasta_raiz)p-values-ar1.tex"
## Exports table of pvalues to .tex
R"""
library(xtable)
table_pvalues = xtable(data.frame(`Probability` = $Alphas, `P-values` = $pvalues))
print(table_pvalues, include.rownames=FALSE, file = $nome_tex)
"""

## Checks for difference on the mean
sum(abs(colMeans(betas_CV ) - 0.3) - abs(colMeans(betas_koenker) - 0.3))


R"""
limites = c(min($betas_CV),max($betas_CV))
dados = $betas_CV; colnames(dados) = $Alphas
boxplot($(keep_betas[:,1,:]), col = "gray", main = "Gamma = 0")
boxplot(dados, ylim = limites, col = "gray", main = "Gamma selection with CV")
"""



output_plot = "$(pasta_raiz)boxplot-ar1.pdf"


@rput Alphas keep_ar output_plot
# Boxplot final comparando cross-validation com gamma = 0 e Ar(1)
R"""
quantiles_ar = quantile(keep_ar,c(0.001,0.25,0.5,0.75,0.999))
library(ggplot2)
t1= cbind.data.frame(melt($betas_CV, varnames = c("Iteration", "Alpha"), value.name = "beta") , Type = "Gamma CV", stringsAsFactors = FALSE)
t2 = cbind.data.frame(melt($(keep_betas[:,1,:]), varnames = c("Iteration", "Alpha"), value.name = "beta"), Type = "Gamma = 0", stringsAsFactors = FALSE)
tnew = rbind(t1,t2)
ggplot(data = tnew, aes(x = factor(Alpha, labels = Alphas), y = beta)) + 
    geom_boxplot(aes(fill = factor(Type)), width = 0.8) + theme_bw() + xlab(expression(alpha)) + ylab(expression(beta)) +
    scale_fill_discrete(name = "Type") + geom_abline(intercept = quantiles_ar[1], slope = 0, linetype = 2, colour = "grey") +
    geom_abline(intercept = quantiles_ar[2], slope = 0, linetype = 2) + geom_abline(intercept = quantiles_ar[3], slope = 0, linetype = 2, colour = "blue") +
    geom_abline(intercept = quantiles_ar[4], slope = 0, linetype = 2) + geom_abline(intercept = quantiles_ar[5], slope = 0, linetype = 2, colour = "grey")

ggplot(data = tnew, aes(x = factor(Alpha, labels = Alphas), y = beta)) + 
    geom_boxplot(aes(fill = factor(Type)), width = 0.8) + theme_bw() + xlab(expression(alpha)) + ylab(expression(beta)) +
    scale_fill_discrete(name = "Type") 
    
"""




R"""

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
limites = c(min($betas_CV),max($betas_CV))

p1 = ggplot(data = tnew, aes(x = factor(Alpha, labels = Alphas), y = beta)) + 
geom_boxplot(aes(fill = factor(Type)), width = 0.8) + theme_bw() + xlab(expression(alpha)) + ylab(expression(beta)) +
scale_fill_discrete(name = "Type") + ylim(limites) + geom_abline(intercept = 0.7, slope = 0)
p2 = qplot(y=$keep_ar, x= 1, geom = "boxplot")
multiplot(p1, p2, cols = 2)


"""

########################## Gráficos ########################
# keep_betas[:,1,:]
n =  400; n_iter = 50
nome_arquivo = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Analise-histograma-coeficients/variables_$(n)_$(n_iter).jld"
keep_betas, keep_APD, keep_ar = load(nome_arquivo, "keep_betas", "keep_APD", "keep_ar")

Alphas = collect(0.05:0.05:0.95)
vector_gamma = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0,20.0] 
# # 


## Boxplots de gammas 
R"""
library(reshape2); library(lattice)
keep_betas_melt = melt($keep_betas, varnames = c("Iteration", "Gamma", "Alpha"), value.name = "beta")
bwplot(beta ~ factor(Alpha, labels = $Alphas) | factor(Gamma, labels = $vector_gamma), data = keep_betas_melt, xlab = $Alphas)
"""
# Análise:
# De modo geral, o erro dos coeficientes é aproximadamente o mesmo para todos 
# os valores de gamma




### Simulação ARIMA




# using StatPlots
# violin(map(x-> string(x), Alphas)',matrix_betas)
# violin(["1" "2"], matrix_betas[:,1:5])

# rcopy(["1" "2"], matrix_betas)

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



## Boxplots final
R"""
library(reshape2); library(lattice)
betas_CV_melt = melt($betas_CV, varnames = c("Iteration", "Alpha"), value.name = "beta")
bwplot(beta ~ factor(Alpha), data = betas_CV_melt, xlab = $Alphas)
"""



### Simulação ARIMA






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

