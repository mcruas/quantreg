function procedure_distancia(L_lasso, L_mip, rho)

    # L_mip = [1,2,3]
    # L_lasso = [1,3,4]
    K = size(L_mip)[1]

    # rho = [1.0 0.7 0.8 0.3;
    #     0.7 1.0 0.5 0.2;
    #     0.8 0.5 1.0 0.1;
    #     0.3 0.2 0.1 1.0]


    m = Model(solver = GurobiSolver())
    @variable(m, 0 <= delta[1:K, 1:K] <= 1)

    @objective(m, Min, (1/K) * sum( delta[i,j] * (1 - abs(rho[L_mip[i],L_lasso[j]]))   for i=1:K, j = 1:K  ) )

    @constraint(m, arestai[i = 1:K], sum(delta[i,j] for j = 1:K) == 1 )
    @constraint(m, arestaj[j = 1:K], sum(delta[i,j] for i = 1:K) == 1 )
    # print(m)
    status = solve(m)

    matriz_delta = getvalue(delta)
    distancia = getobjectivevalue(m)
    return distancia, matriz_delta
end


# P = size(sel_mip)[1]

# L_mip = find(sel_mip)
# L_lasso = find(sel_lasso)
# L_nao_mip = find(!sel_mip)
# L_nao_lasso = find(!sel_lasso)

# K = size(L_mip)[1]

# m = Model(solver = GurobiSolver())
# @variable(m, 0 <= delta[1:P, 1:P] <= 1)

# @objective(m, Max, (1/K)* sum( delta[i,j] * abs(rho[i,j])   for i=1:P, j = 1:P  ) )

# @constraint(m, arestai[i = L_mip], sum(delta[i,j] for j = L_lasso) == 1 )
# @constraint(m, arestaj[j = L_lasso], sum(delta[i,j] for i = L_mip) == 1 )
# @constraint(m, zerari[i = L_nao_mip, j = 1:P], delta[i,j] == 0 )
# @constraint(m, zerarj[j = L_nao_lasso, i = 1:P], delta[i,j] == 0 )
# # print(m)
# status = solve(m)

# matriz_delta = getvalue(delta)
# distancia = 1- getobjectivevalue(m)
# return distancia, matriz_delta
# end
