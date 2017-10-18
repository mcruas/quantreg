# include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/convergencia-dist-rodar.jl")

using JuMP, DataFrames, Distributions, Plots, RCall, Interpolations, Dierckx #, Distributions
pyplot()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
cd(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/")
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl")
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl")
include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/convergencia-dist.jl")

# ns = [300,600, 1000]; # size of data
# rhos = [0.3, 0.9];
# sigmas = [0.1, 0.3,1];
# lambda1s = [0];
# lambda2s = [1];
# non_crosss = [true];
#
# ns = [20]; # size of data
# rhos = [0.3];
# sigmas = [0.1, 0.3,1];
# lambda1s = [0];
# lambda2s = [0.01,0.1, 0.3, 1];
# non_crosss = [true];


ns = collect(50:50:3000)
# ns = collect(50:50:151)
rhos = [0.3, 0.9];
sigmas = [0.1, 0.3,1];
lambda1s = [0];
lambda2s = [1];
non_crosss = [true];

# lambda2s = [0.3];
# non_crosss = [true];
# for lambda2 in lambda2s, n in ns, rho in rhos, sigma in sigmas, lambda1 in lambda1s, non_cross in non_crosss
#   print("n:$n, lambda2:$lambda2, rho:$rho, sigma:$sigma, non_cross:$non_cross\n")
#   # include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/convergencia-dist.jl")
#   try
#     convergencia_dist(lambda2, n , rho , sigma , lambda1 , non_cross)
#   catch
#     print("error!")
#   end
# end


n_iter = 100
# n	l2	rho	sigma	l1	num_erros	erros par	erros np
results = zeros(1000, 8) # Parameters on first 5 columns; n lambda2 rho sigma lambda1 num_erros error_par error_np
i = 1
############# average error over a certain number of observations ##############
for lambda2 in lambda2s, n in ns, rho in rhos, sigma in sigmas, lambda1 in lambda1s, non_cross in non_crosss
  print("n:$n, lambda2:$lambda2, rho:$rho, sigma:$sigma, non_cross:$non_cross\n")
  # include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/convergencia-dist.jl")
  temporary_results = zeros(n_iter, 2)
  for iter in 1:n_iter
      try
          temporary_results[iter,:] = convergencia_dist(lambda2, n , rho , sigma , lambda1 , non_cross = non_cross,
                                                                iter = iter, plotar = false, methods = "par")
      catch
          print("error!")
          temporary_results[iter,:] = [-1,-1]
      end
  end
  results[i, 1:5] = [n, lambda2, rho, sigma, lambda1]
  num_error = sum(temporary_results[:,1] .== -1);
  means = (sum(temporary_results, 1) + num_error)/(n_iter-num_error)
  results[i, 6] = num_error
  results[i, 7:8] = means
  i = i + 1
  # saves files every 20 iterations
  if rem(i,20) == 0
      writecsv("Resultados convergencia linear.csv", results)
  end
end
writecsv("Resultados convergencia linear.csv", results)

sigma = 0.3
results = readcsv("Resultados convergencia linear.csv")
plot(xlab = "Length of time series", ylab = "Error")
for rho in rhos, sigma in [0.3]
# for rho in rhos, sigma in sigmas
  filtrado = results[(results[:,3] .== rho) .* (results[:,4] .== sigma), 7]
  # print("\n rho:$rho, sigma:$sigma",filtrado)
  plot!(ns,filtrado, label = "rho:$rho, sigma:$sigma")
end
plot!()
savefig("Documento Regressao Quantilica/Figuras/convergencia-dist/Convergence.pdf")
# results[i, 1:5] = [n, lambda2, rho, sigma, lambda1]
# num_error = sum(temporary_results[:,1] .== -1);
# means = (sum(temporary_results, 1) + num_error)/(n_iter-num_error)
# results[i, 6] = num_error
# results[i, 7:8] = means
#



#
# ns = [10000]; # size of data
# rhos = [0.3];
# sigmas = [0.3,1];
# # sigma = 0.3
# lambda1s = [0];
# lambda2s = [0.001, 0.003, 0.01];
# non_crosss = [false, true];
#
# include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/convergencia-dist.jl")
