using Winston, RCall
Pkg.add("RCall")
rcopy("1+2")

cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/Documento Regressao Quantilica/Figuras")

x = RCall.rcopy("arima.sim(list(order = c(1,0,2), ar = 0.7, ma = c(0.2, 0.6)), n = 200)")
n = size(x)[1]
plot(x, "b-", ylabel = "y_t", xlabel = "t")
savefig("exemplo-yt.pdf", height = 40)
plot(x[2:n], x[1:(n-1)], "r.", ylabel = "y_t", xlabel = "y_{t-1}")
savefig("exemplo-ar.pdf")

#plot(x,y, "m." , xlabel = "y_t", ylabel = "y_{t-1}")


