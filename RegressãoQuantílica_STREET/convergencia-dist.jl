function convergencia_dist(lambda2, n , rho , sigma , lambda1 ; non_cross = true, iter = 1, plotar = true, methods = "both")
############################################
#
#
# iter: to loop over and produce many different instances
# resume: if true, iter is the number of iterations and only a metric
#         of fit is returned
# methods: "np", "par" or "both"
############################################


## include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/convergencia-dist.jl")
# using JuMP, DataFrames, Distributions, Plots, RCall, Interpolations, Dierckx #, Distributions
# gr()
# usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/")
# include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl")
# include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl")


# Verificar:
# - numero de lags
# - tipos de funcoes incluidas
# - os dados gerados
# - o nome dos arquivos

################ Parameters #############################
# n = 100; # size of data
# rho = 0.3;
# sigma = 0.3
# lambda1 = 0; lambda2 = 0.1;
# non_cross = true;
# iter = 1

srand(convert(Int64,(round(n + 1000*rho + 1000*sigma))) + iter)

Alphas = collect(0.025:0.05:0.975)
push!(Alphas, 0.98, 0.99); unshift!(Alphas, 0.01, 0.02)

x_new = 0;
y0 = 0;

##########################################################

# tic()

residuostmp = rand(LogNormal(0, sigma), n);
# mean(residuostmp)
residuos = residuostmp .- exp(sigma^2/2); # relocate lognormal distribution so it has zero mean
# mean(residuos)
################ Generate data - AR(1) ##################
serie = zeros(n); serie[1] = y0;
for i in 2:n
  serie[i] = serie[i-1] * rho + residuos[i]
end

# plot(serie)

#########################################################
y = serie[2:end]
x = serie[1:end-1]
X = zeros(length(x),1); X[:,1] = x;
alpha_plot = collect(0:0.001:1);

range_errors = collect(0.025:0.05:0.975)
true_quantile = quantile(LogNormal(0,sigma),range_errors) .- exp(sigma^2/2)

if (methods == "both") || (methods == "np")
  thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2, non_cross = non_cross)
  Q_hat_np  = Estimar_Q_hat_np2(x_ord,thetas,Alphas, x_new)
  q_plot_np = Q(alpha_plot, Q_hat_np, Alphas)
  error_np = sqrt(mean((true_quantile .- Q(range_errors, Q_hat_np, Alphas)).^2))
end

if (methods == "both") || (methods == "par")
  betas0, betas = rq_par(y,X,Alphas, non_cross = non_cross)
  Q_hat_par = Estimar_Q_hat_par2(betas0, betas, x_new)
  q_plot_par = Q(alpha_plot, Q_hat_par, Alphas)
  error_par = sqrt(mean((true_quantile .- Q(range_errors, Q_hat_par, Alphas)).^2))
end


# plot(alpha_plot, q_plot)
#




# xlabel = "y_{t-1}", ylabel = "y_t",
if plotar == true
  pyplot()

  if (methods == "both") || (methods == "np")
    p = scatter(x_ord,y_ord, leg = false,  grid = false);
    p = plot!(x_ord ,thetas)
  end

  if (methods == "both") || (methods == "par")
    r = scatter(X, y, leg = false, grid = false) #"y_{t-1}", ylabel = "y_t",
    range_x = collect(minimum(X):0.01:maximum(X))
    r = plot!(range_x, (betas0' * ones(length(range_x))'  + betas' * range_x')')
  end


  s = plot(alpha_plot[2:end-1],quantile(LogNormal(0,sigma), alpha_plot[2:end-1]) .- exp(sigma^2/2), label = "True Distribution")
  if (methods == "both") || (methods == "np")
    s = plot!(alpha_plot, q_plot_np, label = "NP $error_np", legend = :topleft)
  end
  if (methods == "both") || (methods == "par")
    s = plot!(alpha_plot, q_plot_par, label = "Par $error_par")
  end
  # plot(p,r, layout = (2,1))
  if methods == "both"
    plot(plot(p,r),s, layout = (2,1) )
  elseif methods == "np"
    plot(p,s)
  else
    plot(r,s)
  end

  savefig("Documento Regressao Quantilica/Figuras/convergencia-dist/experimentos/$n-$lambda2-$rho-$sigma-$non_cross-$methods.pdf")
end

if methods == "both"
  return [error_par, error_np]
elseif methods == "np"
  return [error_np, error_np]
else
  return [error_par, error_par]
end



end
