# for lambda in [1,10,30,100]
#      include("npqar.jl")
#       end



usesolver = "mosek"    # Escolher entre os valores 'mosek' ou 'gurobi'
cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
include("funcoes_npqar.jl")

using JuMP, DataFrames, Winston
x = convert(Matrix,readtable("x.csv", header = false))
y = convert(Matrix,readtable("y.csv", header = false))


p = 1 # lag
n = length(x)
I = (p+2):n
I2 = 1:n
#lambda = 1

#alphas = [0.05, 0.25, 0.5, 0.75, 0.95]
#alphas = 0.05:0.05:0.95
alphas = 0.01:0.01:0.99
thetasopt = zeros(n, length(alphas))

##################### Inicio laço principal #######################
for alf in 1:length(alphas)

	m = Model(solver = solvfunc)
	alpha = alphas[alf]


	@defVar(m, delta_mais[I] >= 0)
	@defVar(m, delta_menos[I] >= 0)
	@defVar(m, theta[1:n])
	@defVar(m, D2_theta[I2])

	@setObjective(m, Min, sum{alpha *  delta_mais[i] + (1-alpha) *
					delta_menos[i], i = I } + lambda * sum{ D2_theta[i]  , i=I})
	# D2_theta ser igual ao valor absoluta da segunda diferenca de theta
	@addConstraint(m, absolut_posit[i = 3:n], D2_theta[i] >=
		( ( (theta[i]-theta[i-1])/(x[i]-x[i-1]) - (theta[i-1]-theta[i-2])/(x[i-1]-x[i-2]) ) / (x[i]-x[i-1]) ))

	@addConstraint(m, absolut_negat[i = 3:n], D2_theta[i] >= -( ( (theta[i]-theta[i-1])/(x[i]-x[i-1]) - (theta[i-1]-theta[i-2])/(x[i-1]-x[i-2]) ) / (x[i]-x[i-1]) ) )

	# Dar valores ao delta_mais e delta_menos
	@addConstraint(m, deltas[i = I], delta_mais[i] - delta_menos[i] == y[i] - theta[i])

	status = solve(m)

	thetasopt[:, alf] = getValue(theta)
end

plot(x,y, "m." , xlabel = "y_t", ylabel = "y_{t-1}")
for alf in 1:length(alphas)
#	print(alf)
	oplot(x,thetasopt[:,alf], "b-")
end
#oplot()
savefig(string(lambda) * ".pdf")

savefig("exemplo cruzar quantile.pdf")
