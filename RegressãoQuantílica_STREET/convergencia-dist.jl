function convergencia_dist(lambda2, n , rho , sigma , lambda1 , non_cross)
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
mode_distribution = "par"   # choose between "par" for parametric
                            # and "np" for nonparametric


################ Parameters #############################
n = 10; # size of data
rho = 0.3;
sigma = 0.3
lambda1 = 0; lambda2 = 10;
non_cross = true;

Alphas = collect(0.025:0.05:0.975)
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

thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2, non_cross = non_cross)
 #Q_hat =  Estimar_Q_hat_np2(x_ord, thetas, Alphas, x_new; degree_splines = 2)

# scatter(x_ord,y_ord, leg = false, xlabel = "y_{t-1}", ylabel = "y_t", grid = false)
# plot!(x_ord ,thetas)


 alpha_plot = collect(0:0.001:1);

 Q_hat  = Estimar_Q_hat_np2(x_ord,thetas,Alphas, x_new)
 q_plot = Q(alpha_plot, Q_hat, Alphas)
 # plot(alpha_plot, q_plot)
#
true_quantile = quantile(LogNormal(0,sigma),alpha_plot[2:end-1]) .- exp(sigma^2/2)

p = scatter(x_ord,y_ord, leg = false, xlabel = "y_{t-1}", ylabel = "y_t", grid = false);
p = plot!(x_ord ,thetas)
s = plot(alpha_plot, q_plot, leg = false)
s = plot!(alpha_plot[2:end-1],quantile(LogNormal(0,sigma),alpha_plot[2:end-1]) .- exp(sigma^2/2))

plot(p,s)

savefig("Documento Regressao Quantilica/Figuras/convergencia-dist/experimentos/$n-$lambda2-$rho-$sigma-$non_cross.pdf")

# Erro =


end
