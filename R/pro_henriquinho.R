library(rjulia) # https://github.com/armgong/rjulia


#### Chama código em Júlia para fazer a otimização
julia_init()

n <- 30; p <- 30
P <- runif(n,1,4)
Q <- matrix(runif(n*p), 30)
r2j(matriz_correlacao, "rho")
r2j(L_mip, "L_mip")
r2j(L_lasso, "L_lasso")
# r2j(arquivo_procedimento, "arquivo_procedimento")
# Código júlia
julia_void_eval("print(convert(Array{Int64,1},L_lasso))")
julia_void_eval("@show size(L_lasso)[1]")

julia_void_eval("using JuMP, Gurobi;\\
                      @show  pwd();\\
                      include(\"R/procedure_distancia_R.jl\");\\
                      @show rho;\\
                      L_lasso = convert(Array{Int64,1},L_lasso); L_mip = convert(Array{Int64,1},L_mip);\\
                      distancia, matriz_delta = procedure_distancia(L_lasso, L_mip, rho);") #                
saidas <- j2r("(distancia, matriz_delta)")              