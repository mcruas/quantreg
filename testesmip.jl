using JuMP
using Gurobi
function blau1()

    m = Model(solver = GurobiSolver())
    @variable(m, 0 <= x <= 2 )
    @variable(m, 0 <= y <= 30 , Int)

    @objective(m, Max, 5x + 3*y )
    @constraint(m, 1x + 5y <= 3.0 )

    print(m)

    status = solve(m)
    return m
end

asdfg = blau1()
print(asdfg)