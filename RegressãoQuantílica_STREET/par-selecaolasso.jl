# include("par-selecaolasso.jl")

# Verificar:
# - numero de lags
# - tipos de funcoes incluidas
# - os dados gerados
# - o nome dos arquivos
tic()

usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
cd("/home/mcruas/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
# cd("C:/Users/mcruas/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")


include("funcoes_npqar.jl")

using JuMP, DataFrames, Winston, RCall #, Distributions

x = convert(Matrix,readtable("icaraizinho.csv", header = false))

# x = rcopy("arima.sim(list(order = c(6,0,3), ar = c(0.0,0.2,0.0, 0.05, 0.0, 0.2), ma = c(0.2, 0.6, 0.3)), n = 200)")

P = 12 # number of lags to include

# takes a time series and transforms its lags

n = length(x)
#y = rand(Normal(), n)
I = (P+1):(n)
I1 = 1:(n-P)
y = x[I]
# lambdas = 10:10:500
lambdas = e .^(-2:0.005:6.9)
# lambda = 500

# Create matrix with all p lags.
x_estim = lagmatrix(x,12)

# plot(x_estim[:,1], y, "r.")


n_var = (P)

keep_betas = zeros(P+1,length(lambdas))  # matrix to keep values of beta estimated for


alphas = [0.95] ;   Alf = 1:length(alphas)
alpha = alphas[1]

ilamb=1; lambda = 0.000000001
for ilamb in 1:length(lambdas) # Itera sob os valores de lambda
	lambda = lambdas[ilamb]

	##################### Inicio laço principal #######################
	m = Model(solver = solvfunc)

	@defVar(m, beta[1:n_var, Alf])
	@defVar(m, beta0[Alf])
	@defVar(m, delta_mais[I1, Alf] >= 0)
	@defVar(m, delta_menos[I1, Alf] >= 0)
	@defVar(m, abs_beta[1:n_var, Alf] >= 0)

	# Objective Function
	@setObjective(m, Min, sum{alphas[j] *  delta_mais[i, j] + (1-alphas[j]) *
						delta_menos[i, j], i = I1, j = Alf } + lambda *
						sum{abs_beta[i,a] , i = 1:n_var, a = Alf}
						)

	# Absolute value of beta contraints
	@addConstraint(m, absolute_beta_pos[i = 1:n_var, a = Alf],
		abs_beta[i,a] >= beta[i,a])
	@addConstraint(m, absolute_beta_neg[i = 1:n_var, a = Alf],
		abs_beta[i,a] >= - beta[i,a])

	########## Evitar cruzamento de quantis
	#@addConstraint(m, evita_cross[i = I1, j = 2:length(alphas)], theta[i,j] >= theta[i,j-1])

	# Dar valores ao delta_mais e delta_menos
    @addConstraint(m, deltas[i = I1, a = Alf], delta_mais[i,a] - delta_menos[i,a] == y[i] - beta0[a] - sum{beta[l,a] * x_estim[i,l], l = 1:n_var})

	status = solve(m)

	betasopt = getValueMatrixCoco(beta, 12, 1)
	beta0opt = getValueVectorCoco(beta0, 1)
	delta_mais_opt = getValueMatrixCoco(delta_mais, 360,1)
	delta_menos_opt = getValueMatrixCoco(delta_menos, 360,1)



	# plot(delta_mais_opt[1:40], "r.")
	# oplot(delta_menos_opt[1:40], "b.")

	# (1-alpha) * sum(delta_menos_opt) - alpha * sum(delta_mais_opt)

	# loss = delta

	keep_betas[2:(P+1), ilamb] = betasopt
	keep_betas[1, ilamb] = beta0opt[1]

	# betasopt
	# R"  # INICIO R
	# $x_estim
	# variaveis.incluidas <- $betasopt > 0.0001
	# fit <- lm($y ~ $x_estim[, variaveis.incluidas])
 #  	-BIC(fit)
 #  	"  # FIM R


	# plot(x_estim[:,1],y, "r." , xlabel =  "y_{t-1}", ylabel = "y_t")
	# for alf in 1:length(alphas)
	# #	print(alf)
	# 	oplot(x[I],thetasoptMat[I,alf], "b-")
	# end
	# oplot(x,y, "r.")
	# savefig("arima - crossing -" * string(lambda) * ".pdf")

#savefig("sem ccruzar quantile.pdf")
end # Fim do laço dos lambdas

tablebetas = keep_betas[2:13,:]'
plot(log10(lambdas), tablebetas, "-", xlabel = "log_{10}(\\lambda)", ylabel = "size of coefficients", title = "Alpha = $alpha")
legend(["y_{t-1}","y_{t-2}","y_{t-3}","y_{t-4}","y_{t-5}","y_{t-6}","y_{t-7}","y_{t-8}","y_{t-9}","y_{t-10}","y_{t-11}","y_{t-12}"])
alpha_out = replace("$alpha", ".", "")
savefig("par-sellasso-$alpha_out.pdf")

writetable("table-betas-selecaolasso-alpha-$alpha_out.csv", convert(DataFrame, keep_betas), header = false)

toc()




