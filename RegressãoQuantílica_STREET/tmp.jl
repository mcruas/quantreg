### This functions estimates values of beta_0 for the process to be a QAR
### As input are a range of values of X and a the sequence of coefficients beta(alpha)
# β = 0.3 - 0.6 * (Alphas .- (Alphas .^ 2)) 
using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, RCall, Dierckx, JLD
unicodeplots()
usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");


Alphas = collect(0.01:0.01:0.99)
Alphas2 = collect(0.01:0.2:0.99)


# β = 0.3 - 0.6 * (Alphas .- (Alphas .^ 2)) 
β = 0.3 - 0.6 * (Alphas .- (Alphas .^ 2)) 
β = reshape(β, 1, length(β))
betas = β
p_qr = plot(Alphas, betas', legend = false, xlab = "\$\\alpha\$", ylab = "\$\\beta_{p} (\\alpha)\$", title = "Quantile Regression")            

β2 = 0.3 - 0.6 * (Alphas2 .- (Alphas2 .^ 2)) 
β2 = reshape(β2, 1, length(β2))
betas2 = β2



seq_X = 1:0.01:8
X = seq_X[:,:]



function SimB0Qar(X::Array{Float64,2}, Alphas, β::Array{Float64,2})

  J = 1:length(Alphas)
  Jm1 = 2:length(J) # it is the same indexes of Alpha but without the last observation.
  T = 1:size(X)[1]
  P = 1:size(X)[2]

      # Defines the optimization model
  m = Model(solver = GurobiSolver(OutputFlag = 0))

      @variable(m, β0[J])
      
      @objective(m, Min,  β0[length(J)] - β0[1] )

      @constraint(m, evita_cross[i = T, j = Jm1], β0[j] + sum(β[p,j] * X[i,p] for p = P) >= β0[j-1] + sum(β[p,j-1] * X[i,p] for p = P))
      @constraint(m, β0[1] + minimum(X) * betas[1]  == minimum(X) )

  global a = time()

  @time status = solve(m; suppress_warnings=true)

  tmp_betas0opt = getvalue(β0)
  betas0 = zeros(length(Alphas))
  for j in 1:length(Alphas)
    betas0[j] = tmp_betas0opt[j]
  end


  # Transforms the data so they are all generated along the same images
  max_Im1 = betas0[length(J)] + betas[length(J)]*maximum(X)
  min_Im1 = betas0[1] + betas[1]*minimum(X)
  betas_transformed = betas .* (maximum(X) - minimum(X))/(max_Im1 - min_Im1)
  betas0_transformed =  minimum(X) .+  (betas0' + minimum(X) .* betas .- min_Im1) .* 
            ((maximum(X)-minimum(X))/(max_Im1 - min_Im1)) - (minimum(X) .* betas_transformed)

  return betas0_transformed, betas_transformed

  
  
end


### Study showing how the granularity of the Alpha set impacts very little on the range of data
# betas0 = SimB0Qar(X, Alphas, β)
# s_qr1 = plot(seq_X, (betas0' * ones(length(seq_X))'  + betas' * seq_X')', leg = false)

# betas0_2 = SimB0Qar(X, Alphas2, β2)
# s_qr2 = plot(seq_X, (betas0_2' * ones(length(seq_X))'  + betas2' * seq_X')', leg = false)

# plot(s_qr1, s_qr2)

# plot(plot(betas0'), plot(betas0_2'))

pyplot()
### Study showing how the granularity of the Alpha set impacts very little on the range of data
betas0,betas = SimB0Qar(X, Alphas, β)
s_qr1 = plot(seq_X, (betas0' * ones(length(seq_X))'  + betas' * seq_X')', leg = false)




seed = 1
n = 400
y0 = 3
# Function to generate values that follow a QAR process
function simqar(Alphas, n, seed, betas0, betas; y0 = 0) # n = 100 ; seed = 123
  n_burnin = 3000
  
  y = zeros(n + n_burnin)
  y[1] = y0
  # create quantile function
  # Alphas = collect(0.05:0.05:0.95)
  Q_hat = zeros(length(Alphas))
  # plot(betas)
  for t = 2:length(y)
      ytm1 = y[t-1]
      Q_hat = (betas0' * ones(length([ytm1]))'  + betas' * [ytm1]')'
      unif = rand(1)
      y[t] = Q(unif, Q_hat[:], Alphas)[1]
  end 
  return y[end-n:end]
end   
serie = simqar(Alphas, n, 2, betas0, betas; y0 = 2)
plot(serie)


X_lags = lagmatrix(serie,0:1)  # Makes a matrix matching y_t and its regressor y_{t-1}
y = X_lags[:,1]; X = X_lags[:, 2:end]; # Assigns y and X values

results_qr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = 0.0, non_cross = true) # estimates QR with quantile regularization
beta0_est = results_qr[:1]; betas_est = results_qr[:2]
# keep_betas[iter,i_gamma, :] = betas[1,:]

results_penqr = rq_par_lasso(y, X, Alphas; lambda = 0, gamma = 0.01, non_cross = true) # estimates QR with quantile regularization
beta0_estpen = results_penqr[:1]; betas_estpen = results_penqr[:2]


range_x = collect(minimum(X):0.01:maximum(X))


scatter(X,y, xlab = "\$y_{t-1}\$", ylab =  "\$y_{t}\$", leg = false, grid = false)
            s_qr = plot!(range_x, (beta0_est' * ones(length(range_x))'  + betas_est' * range_x')')

plot(plot(betas'), plot(betas_est'), plot(betas_estpen'))



##########################
            # Evaluates fit by testing, for each alpha, how close are the number of obs on the expected range
            for i_alpha = 1:length(Alphas)
              alpha = Alphas[i_alpha]
              keep_APD[i_gamma, i_alpha] = abs(sum(q_predicted_fold[:,i_alpha] .> y_test_fold)/length(y_test_fold) - alpha)
              # println("iter $iter    gamma $i_gamma   alpha $i_alpha     APD $(keep_APD[i_gamma, i_alpha])")
          end


          
q_hat = q_predicted_fold[:,i_alpha]
# Evaluates fit by testing, for each alpha, how close are the number of obs on the expected range
function Calculate_APD(y_test, q_hat , Alphas)
  output = zeros(length(Alphas))
  for i_alpha = 1:length(Alphas)
    alpha = Alphas[i_alpha]
    output[i_alpha] = abs(sum(q_hat[:,i_alpha] .> y_test)/length(y_test) - alpha)
  end
  return output
end

Calculate_APD(y_test_fold, q_predicted_fold, Alphas)

zeros(length(Alphas))
for i_alpha = 1:length(Alphas)
  alpha = Alphas[i_alpha]
  keep_APD[i_gamma, i_alpha] = abs(sum(q_predicted_fold[:,i_alpha] .> y_test_fold)/length(y_test_fold) - alpha)
  # println("iter $iter    gamma $i_gamma   alpha $i_alpha     APD $(keep_APD[i_gamma, i_alpha])")
end

