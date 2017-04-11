function procedure_distancia(sel_lasso, sel_mip, rho)

# sel_mip = sel_lasso = trues(12)
# rho = [1.0 0.7 0.8;
#        0.7 1.0 0.5;
#        0.8 0.5 1.0]

P = size(sel_mip)[1]

L_mip = find(sel_mip)
L_lasso = find(sel_lasso)
L_nao_mip = find(!sel_mip)
L_nao_lasso = find(!sel_lasso)



m = Model(solver = GurobiSolver())
@variable(m, 0 <= delta[1:P, 1:P] <= 1)

@objective(m, Max, (1/P)* sum( delta[i,j] * abs(rho[i,j])   for i=1:P, j = 1:P  ) )

@constraint(m, arestai[i = L_mip], sum(delta[i,j] for j = L_lasso) == 1 )
@constraint(m, arestaj[j = L_lasso], sum(delta[i,j] for i = L_mip) == 1 )
@constraint(m, zerari[i = L_nao_mip, j = 1:P], delta[i,j] == 0 )
@constraint(m, zerarj[j = L_nao_lasso, i = 1:P], delta[i,j] == 0 )
# print(m)
status = solve(m)

matriz_delta = getvalue(delta)
distancia = getobjectivevalue(m)
return distancia, matriz_delta
end
