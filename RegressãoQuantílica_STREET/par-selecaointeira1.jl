# include("par-selecaointeira1.jl")

# Verificar:
# - numero de lags
# - tipos de funcoes incluidas
# - os dados gerados
# - o nome dos arquivos


usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
include("funcoes_npqar.jl")

using JuMP, DataFrames, Winston, RCall #, Distributions

x = convert(Matrix,readtable("icaraizinho.csv", header = false))


P = 12 # number of lags to include

# takes a time series and transforms its lags

n = length(x)
#y = rand(Normal(), n)
I = (P+1):(n)
I1 = 1:(n-P)
y = x[I]
# lambdas = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
# lambdas = [200]


# Create matrix with all p lags. 
x_estim = zeros(n-P, P)
for (i in 1:P)
	x_estim[:,i] = x[I - i]
end

# plot(x_estim[:,1], y, ".")



max_K = 11 # number of variables to include on the model
n_var = (P)
Mu = 30 # Maximum value a coeficient may assume



#### These alpha values are those that will be estimated together
# alphas = [0.05, 0.1, 0.5, 0.9, 0.95]
alphas = [0.90] # alphas junto
Alf = 1:length(alphas)
	alpha = "09"

#### These alpha values are those that will enter a loop, each of them
#at once

keep_betas = zeros(P+1,12)


	for max_K in 1:12 # Itera os maiores de K utilizados
		
		##################### Inicio laço principal #######################
		m = Model(solver = solvfunc)

		@defVar(m, z[1:n_var, Alf], Bin)
		@defVar(m, beta[1:n_var, Alf])
		@defVar(m, beta0[Alf])
		@defVar(m, delta_mais[I1, Alf] >= 0)
		@defVar(m, delta_menos[I1, Alf] >= 0)
		
		# Objective Function
		@setObjective(m, Min, sum{alphas[j] *  delta_mais[i, j] + (1-alphas[j]) *
							delta_menos[i, j], i = I1, j = Alf })

		# inclusion or not inclusion of beta_i on the model
		@addConstraint(m, range_beta_inf[i = 1:n_var, a = Alf], - Mu * z[i,a] <= beta[i,a])
		@addConstraint(m, range_beta_sup[i = 1:n_var, a = Alf], beta[i,a] <= Mu * z[i,a])

		# Only includes max_K variables on the model
		@addConstraint(m, max_var[a = Alf], sum{z[i,a], i = 1:n_var} <= max_K)


		########## Evitar cruzamento de quantis
		#@addConstraint(m, evita_cross[i = I1, j = 2:length(alphas)], theta[i,j] >= theta[i,j-1])

		# Dar valores ao delta_mais e delta_menos
	    @addConstraint(m, deltas[i = I1, a = Alf], delta_mais[i,a] - delta_menos[i,a] == y[i] - beta0[a] - sum{beta[l,a] * x_estim[i,l], l = 1:n_var})

		tic(); status = solve(m); toc()

		tmp = getValue(beta)

		betasopt = getValueMatrixCoco(beta, P, 1)
		beta0opt = getValueVectorCoco(beta0, 1)


		keep_betas[2:(P+1), max_K] = betasopt
		keep_betas[1, max_K] = beta0opt[1]
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
	plot(tablebetas, xlabel = "Number or variables included", ylabel = "size of coefficients", title = "Alpha = $alpha")
	# legend(["y_{t-1}","y_{t-2}","y_{t-3}","y_{t-4}","y_{t-5}","y_{t-6}","y_{t-7}","y_{t-8}","y_{t-9}","y_{t-10}","y_{t-11}","y_{t-12}"])
	savefig("par-selinteira-$alpha.pdf")

	writetable("betas-selecaointeira-alpha-$alpha.csv", convert(DataFrame, keep_betas), header = false)





