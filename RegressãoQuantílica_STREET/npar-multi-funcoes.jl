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

# tmp = [rand(5) ["a", "b", "c", "d" ,"e"]]
# sort(tmp, by)
# sortrows(tmp, by=x->(x[1]))

# scatter(x,y)


# lambda2 = 10
# This function receives a vector of dimension 1 as y and a
# vector of regressors (X).
# The function returns coefficients for the following model
#   q_alpha =
using JuMP


function rq_np(y::Array{Float64}, x::Array{Float64}, Alphas, lambda1 = NaN, lambda2 = NaN; range_data = NaN, non_cross = true)

    Alf = 1:length(Alphas)
    n = size(y)[1]
    # ensures that (x_i,y_i) are ordered in terms of x, so that the algorithm can work ; it also
    # adds a residual to ensure that there are no infinite derivatives
    tmp = [(x .+ rand(n) * 0.0001 ) y];
    tmp = sortrows(tmp, by=x->(x[1]));
    x = tmp[:,1];
    y = tmp[:,2];

    # range_data = NaN

    if any(!isnan(range_data))
      incluir = (y .< range_data[2]) .* (x .< range_data[2]) .* (range_data[1] .< y) .* (range_data[1] .< x)
      x = x[incluir]
      y = y[incluir]
    end

    n = size(y)[1] # updates n, in case some values were cut
    T = 1:n
    nQ = size(x)[1]

    p = 1 # uses only one lag
    # takes a time series and transforms its lags

    n = length(x)
    I = (p+2):(n-1)
    I2 = 1:n
    Alf = 1:length(Alphas)
    m = Model(solver = solvfunc)

    	@variable(m, delta_mais[I, Alf] >= 0)
    	@variable(m, delta_menos[I, Alf] >= 0)
    	@variable(m, theta[1:n, Alf])
      @variable(m, gamma[I2, Alf])
    	@variable(m, xi[I2, Alf])

    	@objective(m, Min, sum(Alphas[j] *  delta_mais[i, j] + (1-Alphas[j]) *
    						delta_menos[i, j] for i = I, j = Alf ) + lambda2 * sum(xi[i, j] for i=I, j = Alf) +
                lambda1 * sum(gamma[i, j] for i=I, j = Alf))


    	# Defines expressions
      @expression(m, D1_theta[i=3:n, j=Alf],
         (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) )
       @expression(m, D2_theta[i=3:n, j=Alf],
         ( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) ) / (x[i]-2*x[i-1]+x[i-2]) ))


         # D1_theta ser igual ao valor absoluta da primeira diferenca de theta
      @constraint(m, absolut_posit_1[i = 3:n, j = Alf], gamma[i,j] >= D1_theta[i, j] )
      @constraint(m, absolut_negat_1[i = 3:n, j = Alf], gamma[i,j] >= - D1_theta[i,j])


         # D2_theta ser igual ao valor absoluta da segunda diferenca de theta
  		@constraint(m, absolut_posit_2[i = 3:n, j = Alf], xi[i,j] >= D2_theta[i, j] )
  		@constraint(m, absolut_negat_2[i = 3:n, j = Alf], xi[i,j] >= - D2_theta[i,j])





    	########## Evitar cruzamento de quantis
      if non_cross
    	   @constraint(m, evita_cross[i = 1:n, j = 2:length(Alphas)], theta[i,j] >= theta[i,j-1])
       end

    		# Dar valores ao delta_mais e delta_menos
    	@constraint(m, deltas[i = I, j = Alf], delta_mais[i,j] - delta_menos[i,j] == y[i] - theta[i,j])

    	status = solve(m)

      thetasopt = getvalue(theta)

      thetasoptMat = zeros(n, length(Alphas))
      for i in 1:n , j in 1:length(Alphas)
    	  thetasoptMat[i,j] = thetasopt[i,j]
    	end

    	# scatter(x,y, xlabel = "y_{t-1}", ylabel = "y_t", title = "Lambda $(lambda2)",  leg = false)
    	# for alf in 1:length(Alphas)
    	# 	# print(alf)
    	# 	plot!(x[I],thetasoptMat[I,alf], linewidth = 2)
    	# end
      # plot!()
      return thetasoptMat, x, y
end


# Difference between rq_np and rq_np2 is that rq_np2 doesnt include the first derivative
# lambda1 is kept to
function rq_np2(y::Array{Float64}, x::Array{Float64}, Alphas, lambda1 = NaN, lambda2 = NaN; range_data = NaN, non_cross = true)

    Alf = 1:length(Alphas)
    n = size(y)[1]
    # ensures that (x_i,y_i) are ordered in terms of x, so that the algorithm can work ; it also
    # adds a residual to ensure that there are no infinite derivatives
    tmp = [(x .+ rand(n) * 0.0001 ) y];
    tmp = sortrows(tmp, by=x->(x[1]));
    x = tmp[:,1];
    y = tmp[:,2];

    if any(!isnan(range_data))
      incluir = (y .< range_data[2]) .* (x .< range_data[2]) .* (range_data[1] .< y) .* (range_data[1] .< x)
      x = x[incluir]
      y = y[incluir]
    end

    n = size(y)[1] # updates n, in case some values were cut
    T = 1:n
    nQ = size(x)[1]

    p = 1 # uses only one lag
    # takes a time series and transforms its lags

    n = length(x)
    I = (p+2):(n-1)
    I2 = 1:n
    Alf = 1:length(Alphas)
    m = Model(solver = solvfunc)

    	@variable(m, delta_mais[I, Alf] >= 0)
    	@variable(m, delta_menos[I, Alf] >= 0)
    	@variable(m, theta[1:n, Alf])
    	@variable(m, xi[I2, Alf])

    	@objective(m, Min, sum(Alphas[j] *  delta_mais[i, j] + (1-Alphas[j]) *
    						delta_menos[i, j] for i = I, j = Alf ) + lambda2 * sum(xi[i, j] for i=I, j = Alf))



    	# Defines expressions
       @expression(m, D2_theta[i=3:n, j=Alf],
         ( ( (theta[i,j]-theta[i-1,j])/(x[i]-x[i-1]) - (theta[i-1,j]-theta[i-2,j])/(x[i-1]-x[i-2]) ) / (x[i]-2*x[i-1]+x[i-2]) ))

         # D2_theta ser igual ao valor absoluta da segunda diferenca de theta
  		@constraint(m, absolut_posit_2[i = 3:n, j = Alf], xi[i,j] >= D2_theta[i, j] )
  		@constraint(m, absolut_negat_2[i = 3:n, j = Alf], xi[i,j] >= - D2_theta[i,j])





    	########## Evitar cruzamento de quantis
      if non_cross
    	   @constraint(m, evita_cross[i = 1:n, j = 2:length(Alphas)], theta[i,j] >= theta[i,j-1])
       end

    		# Dar valores ao delta_mais e delta_menos
    	@constraint(m, deltas[i = I, j = Alf], delta_mais[i,j] - delta_menos[i,j] == y[i] - theta[i,j])

    	status = solve(m)

      thetasopt = getvalue(theta)

      thetasoptMat = zeros(n, length(Alphas))
      for i in 1:n , j in 1:length(Alphas)
    	  thetasoptMat[i,j] = thetasopt[i,j]
    	end

    	# scatter(x,y, xlabel = "y_{t-1}", ylabel = "y_t", title = "Lambda $(lambda2)",  leg = false)
    	# for alf in 1:length(Alphas)
    	# 	# print(alf)
    	# 	plot!(x[I],thetasoptMat[I,alf], linewidth = 2)
    	# end
      # plot!()
      return thetasoptMat, x, y
end




# thetas = rq_np(y,x, Alphas, 10)




# This function calculates the value from an estimated quantile function
# alpha is the probability to be calculated ; Q_hat the empirical quantiles ;
# Alphas is a vector of which probabilities
function Q(alpha, Q_hat::Array{Float64}, Alphas::Array{Float64})

  ## First step is to define the value of Q(0) and Q(1)
  Q_0 = 2* Q_hat[1] - Q_hat[2]
  Q_fim = 2*Q_hat[end] - Q_hat[end-1]
  # scatter([0; Alphas; 1] ,
         #  [Q_0 ; Q_hat ; Q_fim], k = 1)

  ## Creates a Spline object
  sp1 = Spline1D([0; Alphas; 1] ,
                 [ Q_0 ; Q_hat ; Q_fim], k = 1)
  ponto = evaluate(sp1, alpha)
  # scatter!([valor],[ponto])
  return ponto

end

# Recebe como input a série, variaveis explicativas, Alphas e um novo valor de X
# Faz a estimação e calcula os valores de Q_hat na discretização dada por Alphas
function Estimar_Q_hat_np(y,x, Alphas, lambda1, lambda2, x_new; degree_splines = 2 , range_data = NaN, non_cross = true)
  thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2, non_cross = non_cross, range_data = range_data)
  Q_hat = zeros(length(Alphas))
  for col_alpha in 1:length(Alphas)
    sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = degree_splines)
    Q_hat[col_alpha] = evaluate(sp1, x_new)
  end
  return Q_hat
end

# Recebe como input os valores estimados de betas, betas0 e o novo valor de x
# Para estimar os valores dos coeficientes, utilizar a função rq_par, abaixo:
# thetas, x_ord, y_ord = rq_np(y,x,Alphas, lambda1, lambda2, non_cross = non_cross, range_data = range_data)
function Estimar_Q_hat_np2(x_ord, thetas, Alphas, x_new; degree_splines = 2)

  Q_hat = zeros(length(Alphas))
  for col_alpha in 1:length(Alphas)
    sp1 = Spline1D( x_ord , thetas[:,col_alpha] , k = degree_splines)
    Q_hat[col_alpha] = evaluate(sp1, x_new)
  end
  return Q_hat
end






# This function receives a vector of dimension 1 as y and a
# matrix (n x q) of regressors (X).
# The function returns coefficients for the following model
#   y = beta0 + X * beta
function rq_par(y::Array{Float64}, X::Array{Float64,2}, Alphas; non_cross = true)

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
    if non_cross
  	   @constraint(m, evita_cross[i = T, j = 2:length(Alphas)], β0[j] + sum(β[q,j] * X[i,q] for q = Q) >= β0[j-1] + sum(β[q,j-1] * X[i,q] for q = Q))
     end

  		# Dar valores ao ɛ_tmais e ao ɛ_tmenos
  	@constraint(m, epsilons[i = T, j = Alf], ɛ_tmais[i,j] - ɛ_tmenos[i,j] == y[i] - β0[j] - sum(β[q,j] * X[i,q] for q = Q))

  	status = solve(m)

  	tmp_betas0opt = getvalue(β0)
    tmp_betasopt = getvalue(β)
    objectiveValue = getobjectivevalue(m)
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


function rq_par_mip(y::Array{Float64}, X::Array{Float64,2}, Alphas; non_cross = true, max_K = NaN, TimeLimit = NaN, MIPGap = NaN)

    Alf = 1:length(Alphas)
    M = 2
    n = size(y)[1]
    T = 1:n
    Q = 1:size(X)[2]

    if isnan(max_K)
      max_K = size(X)[2]
    end

    m = Model(solver = GurobiSolver(MIPGap = MIPGap, TimeLimit = TimeLimit))
  	@variable(m, ɛ_tmais[T, Alf] >= 0)
  	@variable(m, ɛ_tmenos[T, Alf] >= 0)
  	@variable(m, β[Q, Alf])
  	@variable(m, β0[Alf])
    @variable(m, z[Q, Alf], Bin)


    @objective(m, Min, sum(Alphas[j] * ɛ_tmais[i, j] + (1-Alphas[j]) *ɛ_tmenos[i, j] for i = T, j = Alf ))

  	########## Evitar cruzamento de quantis
    if non_cross
  	   @constraint(m, evita_cross[i = T, j = 2:length(Alphas)], β0[j] + sum(β[q,j] * X[i,q] for q = Q) >= β0[j-1] + sum(β[q,j-1] * X[i,q] for q = Q))
     end

  		# Dar valores ao ɛ_tmais e ao ɛ_tmenos
  	 @constraint(m, epsilons[i = T, j = Alf], ɛ_tmais[i,j] - ɛ_tmenos[i,j] == y[i] - β0[j] - sum(β[q,j] * X[i,q] for q = Q))

     # Restringir a estimação a no máximo K valores
     @constraint(m, range_beta_inf[i = Q, a = Alf], - M * z[i,a] <= β[i,a])
     @constraint(m, range_beta_sup[i = Q, a = Alf], β[i,a] <= M * z[i,a])

     # Permite que no máximo max_K seja escolhido no modelo
     @constraint(m, max_var[a = Alf], sum(z[i,a] for i = Q) <= max_K)



  	status = solve(m)

  	tmp_betas0opt = getvalue(β0)
    tmp_betasopt = getvalue(β)
    objectiveValue = getobjectivevalue(m)
    ## Transform both variables into an array
    betasopt = zeros(size(X)[2], length(Alphas))
    betas0opt = zeros(length(Alphas))
    for q in 1:size(X)[2] , j in 1:length(Alphas)
      betasopt[q,j] = tmp_betasopt[q,j]
    end
    for j in 1:length(Alphas)
      betas0opt[j] = tmp_betas0opt[j]
    end


    return betas0opt', betasopt, objectiveValue, status
end



# max_K = 3; TimeLimit = 60; Grupos = 3; MIPGap = 0.0; non_cross = true


function rq_par_mip_grupos(y::Array{Float64}, X::Array{Float64,2}, Alphas; non_cross = true, max_K = NaN, TimeLimit = NaN, MIPGap = NaN, Grupos = NaN)

    Alf = 1:length(Alphas)
    M = 10
    n = size(y)[1]
    T = 1:n
    P = size(X)[2]

    if isnan(max_K)
      max_K = size(X)[2]
    end

    if isnan(Grupos)
      G = size(X)[2]
    else
      G = Grupos
    end

    m = Model(solver = GurobiSolver(MIPGap = MIPGap, TimeLimit = TimeLimit))
  	@variable(m, ɛ_tmais[T, Alf] >= 0)
  	@variable(m, ɛ_tmenos[T, Alf] >= 0)
  	@variable(m, β[1:P, Alf])
  	@variable(m, β0[Alf])
    @variable(m, z[1:P, 1:G], Bin)
    @variable(m, I[1:P, Alf], Bin)


    @objective(m, Min, sum(Alphas[j] * ɛ_tmais[i, j] + (1-Alphas[j]) *ɛ_tmenos[i, j] for i = T, j = Alf ))


  	########## Evitar cruzamento de quantis
    if non_cross
  	   @constraint(m, evita_cross[i = T, j = 2:length(Alphas)], β0[j] + sum(β[p,j] * X[i,p] for p = 1:P) >= β0[j-1] + sum(β[p,j-1] * X[i,p] for p = 1:P))
     end

  		# Dar valores ao ɛ_tmais e ao ɛ_tmenos
  	 @constraint(m, epsilons[i = T, j = Alf], ɛ_tmais[i,j] - ɛ_tmenos[i,j] == y[i] - β0[j] - sum(β[p,j] * X[i,p] for p = 1:P))

     # Restringir a estimação a no máximo K valores
     @constraint(m, range_beta_inf[p = 1:P, a = Alf, g = 1:G], β[p,a] <= M*(2 - (1-z[p,g]) - I[g,a]) )
     @constraint(m, range_beta_sup[p = 1:P, a = Alf, g = 1:G],  - M*(2 - (1-z[p,g]) - I[g,a]) <= β[p,a]  )



     # Permite que no máximo max_K seja escolhido no modelo
     @constraint(m, max_var[g = 1:G], sum(z[p,g] for p = 1:P) <= max_K)

     # Permite que um alpha pertença a apenas 1 Grupo
     @constraint(m, pertence_grupo[a = Alf], sum(I[g,a] for g = 1:G) == 1)




  	status = solve(m)

  	tmp_betas0opt = getvalue(β0)
    tmp_betasopt = getvalue(β)
    objectiveValue = getobjectivevalue(m)
    # objetivo = getobjective()
    ## Transform both variables into an array
    betasopt = zeros(size(X)[2], length(Alphas))
    betas0opt = zeros(length(Alphas))
    for q in 1:size(X)[2] , j in 1:length(Alphas)
      betasopt[q,j] = tmp_betasopt[q,j]
    end
    for j in 1:length(Alphas)
      betas0opt[j] = tmp_betas0opt[j]
    end


    return betas0opt', betasopt, objectiveValue, status
end





function rq_par_mip_grupos_rampa(y::Array{Float64}, X::Array{Float64,2}, Alphas; non_cross = true, max_K = NaN, TimeLimit = NaN, MIPGap = NaN, Grupos = NaN)

    Alf = 1:length(Alphas)
    Alfm = 1:length(Alphas)-1 # it is the same indexes of Alpha but without the last observation.
    M = 5
    
    n = size(y)[1]
    T = 1:n
    P = size(X)[2]
    M2 = P * M
    if isnan(max_K)
      max_K = size(X)[2]
    end

    if isnan(Grupos)
      G = size(X)[2]
    else
      G = Grupos
    end

    m = Model(solver = GurobiSolver(MIPGap = MIPGap, TimeLimit = TimeLimit))
  	@variable(m, ɛ_tmais[T, Alf] >= 0)
  	@variable(m, ɛ_tmenos[T, Alf] >= 0)
  	@variable(m, β[1:P, Alf])
  	@variable(m, β0[Alf])
    @variable(m, z[1:P, Alf], Bin)
    @variable(m, ϕ[1:P, Alfm] >= 0)  # m_α,p indicates the absolute value of z_α,p 
    @variable(m, r[Alfm], Bin)  # m_α,p indicates the absolute value of z_α,p 
        


    @objective(m, Min, sum(Alphas[j] * ɛ_tmais[i, j] + (1-Alphas[j]) *ɛ_tmenos[i, j] for i = T, j = Alf ))


  	########## Evitar cruzamento de quantis
    if non_cross
  	   @constraint(m, evita_cross[i = T, j = 2:length(Alphas)], β0[j] + sum(β[p,j] * X[i,p] for p = 1:P) >= β0[j-1] + sum(β[p,j-1] * X[i,p] for p = 1:P))
     end

  		# Dar valores ao ɛ_tmais e ao ɛ_tmenos
  	 @constraint(m, epsilons[i = T, j = Alf], ɛ_tmais[i,j] - ɛ_tmenos[i,j] == y[i] - β0[j] - sum(β[p,j] * X[i,p] for p = 1:P))

     # Restringir a estimação a no máximo K valores
    #  @constraint(m, range_beta_inf[p = 1:P, a = Alf, g = 1:G], β[p,a] <= M*(2 - (1-z[p,g]) - I[g,a]) )

     @constraint(m, range_beta_inf[p = 1:P, a = Alf], β[p,a] <= M*z[p,a] )
     @constraint(m, range_beta_sup[p = 1:P, a = Alf],  - M*z[p,a] <= β[p,a]  )


     # Limitar número de variáveis por α-quantil
     @constraint(m, max_var[a = Alf], sum(z[p,a] for p = 1:P) <= max_K)


     # Limitar número de rampas
     @constraint(m, n_rampa[a = Alfm], sum(ϕ[p,a] for p = 1:P) <= r[a]*M)


     # Limitar número de rampas
     @constraint(m, max_rampas, sum(r[a] for a = Alfm) <= G - 1)

     # Permite que no máximo max_K seja escolhido no modelo
     @constraint(m, max_var[a = Alf], sum(z[i,a] for i = Q) <= max_K)


     # Coloca |z_pa+1 - z_pa| em phi
     @constraint(m, phi_positivo[p = 1:P, a = Alfm], z[p,a+1] - z[p,a]  <= ϕ[p,a] )
     @constraint(m, phi_negativo[p = 1:P, a = Alf],  - z[p,a+1] + z[p,a]  <= ϕ[p,a]  )



  	status = solve(m)

  	tmp_betas0opt = getvalue(β0)
    tmp_betasopt = getvalue(β)
    objectiveValue = getobjectivevalue(m)
    # objetivo = getobjective()
    ## Transform both variables into an array
    betasopt = zeros(size(X)[2], length(Alphas))
    betas0opt = zeros(length(Alphas))
    for q in 1:size(X)[2] , j in 1:length(Alphas)
      betasopt[q,j] = tmp_betasopt[q,j]
    end
    for j in 1:length(Alphas)
      betas0opt[j] = tmp_betas0opt[j]
    end


    return betas0opt', betasopt, objectiveValue, status

end



# Recebe como input a série, variaveis explicativas, Alphas e um novo valor de X
# Faz a estimação e calcula os valores de Q_hat na discretização dada por Alphas
function Estimar_Q_hat_par(y,x, Alphas, x_new; non_cross = true)
  betas0, betas = rq_par(y,X, Alphas, non_cross = non_cross); # coeficientes do modelo linear
  Q_hat = (betas0 + x_new * betas)[1,:];
  return Q_hat
end

# Recebe como input os valores estimados de betas, betas0 e o novo valor de x
# Para estimar os valores dos coeficientes, utilizar a função rq_par, abaixo:
#     betas0, betas = rq_par(y,X, Alphas, non_cross = non_cross)
function Estimar_Q_hat_par2(betas0,betas, x_new)
  Q_hat = (betas0 + x_new * betas)[1,:];
  return Q_hat
end


# function Plotar_quantis(x,y,thetas, )
