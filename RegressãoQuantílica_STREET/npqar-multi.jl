# include("npqar-multi.jl")
# for hora in 0:23
# 	include("npqar-multi.jl")
# end
# hora = 0:23
# hora = 9

# Verificar:
# - crossing ou nao crossing
# - os valores de alpha e lambda desejados
# - os dados gerados
# - o nome dos arquivos
tic()
usesolver = "mosek"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/
# 																	RegressãoQuantílica_STREET")
cd("/home/mcruas/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
# cd("/home/mcruas/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET/")
# cd("C:/Users/mcruas/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET/")

include("funcoes_npqar.jl")

using JuMP, DataFrames, Plots, RCall, Interpolations #, Distributions

localidade = "saojosemipibu" # Escolher entre saojosemipibu ; tubarao ; tabocas
pasta_dados = "../Dados Climaticos/Solar-$localidade/"



# x = convert(Matrix,readtable("$(pasta_dados)x.csv", header = false))
# y = convert(Matrix,readtable("$(pasta_dados)y.csv", header = false))
#x = convert(Matrix,readtable("x_nonlin.csv", header = false))
#y = convert(Matrix,readtable("y_nonlin.csv", header = false))

############### Retirada de um certo número de pontos #########
# n = length(x)

# R"selecao = as.logical(rbinom($n, 1, 1))"
# @rget selecao
# x = x[selecao]
# y = y[selecao]
#################################################################3


lambdas = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
#lambdas = [50]
lambdas = [100]



###################### Leitura a partir do R #################
caminho_arquivo = "../Dados Climaticos/Solar-$localidade/$localidade.xlsx" # caminho do arquivo contendo os dados
@rput hora caminho_arquivo
R"
  	library('dplyr')
		library('readxl')
		library('lattice')
    library('stringr')
		dados <- read_excel(path = caminho_arquivo)[,1:6]
 		dados_por_hora <- dados %>% filter(y_t0 > 0 , y_t1 > 0, hour %in% hora)
		n = nrow(dados_por_hora)
	  dados_por_hora =	dados_por_hora %>% mutate(y_t0modif = y_t0 + rnorm(n)*0.001) %>% arrange(y_t0modif) %>% select(y_t0modif, y_t1)
  	x = dados_por_hora$y_t0modif
 		y = dados_por_hora$y_t1
"
@rget x y n


# plot(x)
# scatter(x,y)



pyplot()
# gui()
# "teste $lambdas"




## Shows boxplot of months
# tmp_data = readxl("$(pasta_dados)tubarao solar-28.467_-49.005_uncorrected.xlsx","data")


p = 1 # lag

# takes a time series and transforms its lags

n = length(x)
#y = rand(Normal(), n)
I = (p+2):(n-1)
I2 = 1:n

alphas = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
# alphas = 0.05:0.05:0.95
#alphas = 0.01:0.01:0.99


#alphas = [0.05]
Alf = 1:length(alphas)
thetasoptMat = zeros(n, length(alphas))

for ilamb in 1:length(lambdas) # Itera so bos valores de lambda
	lambda = lambdas[ilamb] #  ilamb = 1
	##################### Inicio laço principal #######################
m = Model(solver = solvfunc)

	@variable(m, delta_mais[I, Alf] >= 0)
	@variable(m, delta_menos[I, Alf] >= 0)
	@variable(m, theta[1:n, Alf])
	@variable(m, D2_theta[I2, Alf])

	@objective(m, Min, sum{alphas[j] *  delta_mais[i, j] + (1-alphas[j]) *
						delta_menos[i, j], i = I, j = Alf } + lambda * sum{ D2_theta[i, j]  , i=I, j = Alf})
		# D2_theta ser igual ao valor absoluta da segunda diferenca de theta
	#
	# @addConstraint(m, absolut_posit[i = 3:n, j = Alf], D2_theta[i, j] >=
	# 		( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) )  ))
	# @addConstraint(m, absolut_negat[i = 3:n, j = Alf], D2_theta[i,j] >=
	# 	-( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) )  ) )
	#


		@addConstraint(m, absolut_posit[i = 3:n, j = Alf], D2_theta[i, j] >=
				( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) ) / (x[i]-2*x[i-1]+x[i-2]) ))
		@addConstraint(m, absolut_negat[i = 3:n, j = Alf], D2_theta[i,j] >=
			-( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) ) /  (x[i]-2*x[i-1]+x[i-2]) ))



	########## Evitar cruzamento de quantis
	@addConstraint(m, evita_cross[i = 1:n, j = 2:length(alphas)], theta[i,j] >= theta[i,j-1])

		# Dar valores ao delta_mais e delta_menos
	@addConstraint(m, deltas[i = I, j = Alf], delta_mais[i,j] - delta_menos[i,j] == y[i] - theta[i,j])

	status = solve(m)

	thetasopt = getvalue(theta)

	for i in 1:n , j in 1:length(alphas)
	  thetasoptMat[i,j] = thetasopt[i,j]
	end

	scatter(x,y, xlabel = "y_t-1", ylabel = "y_t", title = "Hora $hora - Lambda $(lambdas[1])",  leg = false)
	for alf in 1:length(alphas)
		# print(alf)
		plot!(x[I],thetasoptMat[I,alf], linewidth = 2)
	end
	# oplot(x,y, "r.")
	# savefig("icaraizinho-crossing-" * string(lambda) * ".pdf")
  return thetasoptMat
#savefig("sem ccruzar quantile.pdf")
end # Fim do laço dos lambdas



# tst()

toc()



# scatter(x,y, xlabel = "y_t", ylabel = "y_{t-1}")
# for alf in 1:length(alphas)
# 	print(alf)
# 	plot!(x[I],thetasoptMat[I,alf])
# end
#
plot!()

# itp = interpolate(thetasoptMat)

savefig("../Documento Regressao Quantilica/Figuras/npqar-solar-$localidade/hora$hora-lambda$(lambdas[1]).png")
#
# plot!(x[I],thetasoptMat[I,1])
# plot!(x[I],thetasoptMat[I,2])
# plot!(x[I],thetasoptMat[I,3])
# plot!(x[I],thetasoptMat[I,4])
# plot!(x[I],thetasoptMat[I,5])
# nomes_series = ["a", "b", "c", "d", "e", "f", "g", "h"]
# plot!(lab = nomes_series)

#
# plot(rand(5,2))
# plot!(lab = ["asdf","asdfsda"])
