# Function to estimate all quantiles at the same time, iwth a noncrossing quantile constraint
# cd("/home/mcruas/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"); usesolver = "mosek"
# y = randn(100)
# x = randn(100,3)
#
# Alphas = 0.01:0.1:0.99
# Alphas = [0.01,0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
# typeof(Alphas)
#

# include("RegressãoQuantílica_STREET/funcoes_npqar.jl")
# using JuMP, DataFrames, Plots, RCall, Interpolations #, Distributions
y
X
# This function receives a vector of dimension 1 as y and a
# matrix (n x q) of regressors (X).
# The function returns coefficients for the following model
#   q_alpha =
function rq_np(y::Array{Float64}, X::Array{Float64,2}, Alphas, lambda = NaN)

    Alf = 1:length(Alphas)

    n = size(y)[1]
    T = 1:n
    Q = 1:size(X)[2]

    m = Model(solver = solvfunc)
  	@variable(m, ɛ_tmais[T, Alf] >= 0)
  	@variable(m, ɛ_tmenos[T, Alf] >= 0)
  	@variable(m, β[Q, Alf])
  	@variable(m, β0[Alf])
  	@objective(m, Min, sum(Alphas[j] * ɛ_tmais[i, j] + (1-Alphas[j]) *ɛ_tmenos[i, j] for i = T, j = Alf ))

  	########## Evitar cruzamento de quantis
  	@constraint(m, evita_cross[i = T, j = 2:length(Alphas)], β0[j] + sum(β[q,j] * X[i,q] for q = Q) >= β0[j-1] + sum(β[q,j-1] * X[i,q] for q = Q))

  		# Dar valores ao ɛ_tmais e ao ɛ_tmenos
  	@constraint(m, epsilons[i = T, j = Alf], ɛ_tmais[i,j] - ɛ_tmenos[i,j] == y[i] - β0[j] - sum(β[q,j] * X[i,q] for q = Q))

  	status = solve(m)

  	tmp_betas0opt = getvalue(β0)
    tmp_betasopt = getvalue(β)

    ## Transform both variables into an array
    betasopt = zeros(size(X)[2], length(Alphas))
    betas0opt = zeros(length(Alphas))
    for q in 1:size(X)[2] , j in 1:length(Alphas)
      betasopt[q,j] = tmp_betasopt[q,j]
    end
    for j in 1:length(Alphas)
      betas0opt[j] = tmp_betas0opt[j]
    end


    return betas0opt', betasopt
end

# This function calculates the value from an estimated quantile function
# alpha is the probability to be calculated ; Q_hat the empirical quantiles ;
# Alphas is a vector of which probabilities
function Q(alpha, Q_hat::Array{Float64}, Alphas::Array{Float64})

  ## First step is to define the falue of Q(0) and Q(1)
  Q_0 = 2* Q_hat[1] - Q_hat[2]
  Q_fim = 2*Q_hat[end] - Q_hat[end-1]
  # scatter([0; Alphas; 1] ,
         #  [Q_0 ; Q_hat ; Q_fim], k = 1)

  ## Creates a Spline object
  sp1 = Spline1D([0; Alphas; 1] ,
                 [ Q_0 ; Q_hat ; Q_fim], k = 2)
  ponto = evaluate(sp1, alpha)
  # scatter!([valor],[ponto])
  return ponto

end

# betas0, betas = rq(y,x, Alphas)
#
#
# @rput betas0 betas x y Alphas
# R"plot(x[,1],y)
#
#  for (alf in 1:length(Alphas)) {
#       abline(a = betas0[alf],b= betas[1 , alf])
#    }
# "
