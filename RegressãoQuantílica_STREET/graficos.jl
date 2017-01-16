using Winston, DataFrames, RCall

cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
# cd("C:/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
# cd("C://Users//mcruas//Dropbox//Pesquisa Doutorado//Paper NPQuantile//RegressãoQuantílica_STREET")

x = convert(Matrix,readtable("icaraizinho.csv", header = false))
n = length(x)


# takes a time series and transforms its lags


#y = rand(Normal(), n)
P = 1; I = (P+2):(n-1)
plot(x[I-P], x[I], "r.", linewidth  = 4,
	title = "Icaraizinho first lag", xlabel = "y_{t-1}", ylabel = "y_t")
savefig("icaraizinho-$P-lag.pdf")

P = 4; I = (P+2):(n-1)
plot(x[I-P], x[I], "r.", linewidth  = 4,
	title = "Icaraizinho fourth lag", xlabel = "y_{t-$P}", ylabel = "y_t")
savefig("icaraizinho-$P-lag.pdf")

P = 11; I = (P+2):(n-1)
plot(x[I-P], x[I], "r.", linewidth  = 4,
	title = "Icaraizinho eleventh lag", xlabel = "y_{t-$P}", ylabel = "y_t")
savefig("icaraizinho-$P-lag.pdf")


P = 12; I = (P+2):(n-1)
plot(x[I-P], x[I], "r.", linewidth  = 4,
	title = "Icaraizinho twelth lag", xlabel = "y_{t-$P}", ylabel = "y_t")
savefig("icaraizinho-$P-lag.pdf")

plot(x[1:12])


R"library(ggfortify)
  x = $x[,1]
  autoplot(ts(x))"



##################### Gráfico lasso #####################
tabela = readtable("table-betas-sellassonorm-alpha-05.csv", header = false)

nvar = colApply((elemApply(tabela, abs) .> 0.0000001), sum)
lambdas = e .^ (-5:0.1:5)

plot(-log(lambdas), nvar, "r", xlabel = "-log_e(\\lambda)", ylabel = "Number of variables included")



##################### Gráfico aproximacao funcao de distribuicao ##############################
using Plots
