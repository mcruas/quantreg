# include("npqar-multi.jl")

# Verificar:
# - crossing ou nao crossing
# - os valores de alpha e lambda desejados
# - os dados gerados
# - o nome dos arquivos
tic()
usesolver = "glpk"    # Escolher entre os valores 'mosek' ou 'gurobi'
cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
include("funcoes_npqar.jl")

using JuMP, DataFrames, Winston, RCall #, Distributions


# Nao linear
# rcopy(" set.seed(1)
# 		n = 1001
# 		n1 = (n-1)/2
# 		x1 = seq(0,1,1/n1)
# 		x2 = seq(1+1/n1,2,1/n1)
# 		x = c(x1,x2)
# 		y1 = 1 + rnorm(n/2,sd = 0.2)
# 		y2 = x2*1 + rnorm(n/2,sd = 0.2)
# 		y = c(y1,y2)")

# Arima
# rcopy("		set.seed(1)
# 			yt.ini = arima.sim(list(order = c(4,0,3), ar = c(0.3,0.02,0.02, 0.3), ma = c(0.2, 0.6, 0.3)), n = 1000)
# 			# plot(yt.ini)
# 			xt = yt.ini[-length(yt.ini)]; yt = yt.ini[-1]
# 			ordem = order(xt)
# 			x = xt[ordem]
# 			y = yt[ordem]
# 			# plot(x,y)
# 			")

# Normal bivariada
# rcopy("
# set.seed(1)
# xt = rnorm(500); yt=rnorm(500)
# ordem = order(xt)
# x = xt[ordem]
# y = yt[ordem]
# ")


# x = rcopy("x")
# y = rcopy("y")
# plot(x,y, ".")
x = convert(Matrix,readtable("x.csv", header = false))
y = convert(Matrix,readtable("y.csv", header = false))
#x = convert(Matrix,readtable("x_nonlin.csv", header = false))
#y = convert(Matrix,readtable("y_nonlin.csv", header = false))
plot(x,y, ".")

p = 1 # lag

# takes a time series and transforms its lags

n = length(x)
#y = rand(Normal(), n)
I = (p+2):(n-1)
I2 = 1:n
# lambdas = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
lambdas = [50]

lambda = 0.1
alphas = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
#alphas = 0.05:0.05:0.95
#alphas = 0.01:0.01:0.99

#alphas = [0.05]
Alf = 1:length(alphas)

for ilamb in 1:length(lambdas) # Itera so bos valores de lambda
	lambda = lambdas[ilamb] #  ilamb = 1
	##################### Inicio laço principal #######################
	m = Model(solver = solvfunc)

	@defVar(m, delta_mais[I, Alf] >= 0)
	@defVar(m, delta_menos[I, Alf] >= 0)
	@defVar(m, theta[1:n, Alf])
	@defVar(m, D2_theta[I2, Alf])

	@setObjective(m, Min, sum{alphas[j] *  delta_mais[i, j] + (1-alphas[j]) *
						delta_menos[i, j], i = I, j = Alf } + lambda * sum{ D2_theta[i, j]  , i=I, j = Alf})
		# D2_theta ser igual ao valor absoluta da segunda diferenca de theta

	#@addConstraint(m, absolut_posit[i = 3:n, j = Alf], D2_theta[i, j] >=
	#		( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) ) / (x[i]-x[i-1]) ))
	@addConstraint(m, absolut_posit[i = 3:n, j = Alf], D2_theta[i, j] >=
			( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) ) ))


	@addConstraint(m, absolut_negat[i = 3:n, j = Alf], D2_theta[i,j] >= 
		-( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) )  ) )

	########## Evitar cruzamento de quantis
	#@addConstraint(m, evita_cross[i = 1:n, j = 2:length(alphas)], theta[i,j] >= theta[i,j-1])

		# Dar valores ao delta_mais e delta_menos
	@addConstraint(m, deltas[i = I, j = Alf], delta_mais[i,j] - delta_menos[i,j] == y[i] - theta[i,j])

	status = solve(m)

	thetasopt = getValue(theta)

	thetasoptMat = zeros(n, length(alphas))
	for i in 1:n , j in 1:length(alphas)
	  thetasoptMat[i,j] = thetasopt[i,j]
	end
	plot(x,y, "r." , xlabel = "y_t", ylabel = "y_{t-1}")
	for alf in 1:length(alphas)
	#	print(alf)
		oplot(x[I],thetasoptMat[I,alf], "b-")
	end
	oplot(x,y, "r.")
	savefig("icaraizinho-crossing-" * string(lambda) * ".pdf")

#savefig("sem ccruzar quantile.pdf")
end # Fim do laço dos lambdas


end

# tst()

toc()