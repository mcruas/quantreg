using JuMP
using Gurobi

# We will use Gurobi and disable PreSolve, Cuts, and (in-built) Heuristics so
# only our heuristic will be used
m = Model(solver=GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0))

n = 100000; V = 3
# Define our variables to be inside a box, and integer
w = rand(n)
v = rand(n)
@variable(m, x[1:n], Bin);

# Optimal solution is trying to go towards top-right corner (2.0, 2.0)
@objective(m, Max, sum(w[i]*x[i] for i = 1:n));

# @add

# We have one constraint that cuts off the top right corner
@constraint(m, sum(v[i]*x[i] for i=1:n) <= V);


##### Pegar soluções intermediárias ###########
# solutionvalues = Vector{Float64}[]
# function infocallback(cb)
#     push!(solutionvalues, JuMP.getvalue(x))
# end
# addinfocallback(m, infocallback, when = :MIPSol)
################################################


##### Pegar resultados #########################
type NodeData_tst
         time::Float64  # in seconds since the epoch
        x::Vector{Float64}
         bestbound::Float64
         obj::Float64

end
 
function infocallback(cb)
    # println(a - time())
    i = find(check .== 0)[1]
    # println ("Tempo $(time()-a) e $(times[i])")
    if time()-a > times[i]
        obj       = MathProgBase.cbgetobj(cb)
        bestbound = MathProgBase.cbgetbestbound(cb)
        solution = JuMP.getvalue(x)
        # push!(bbdata2, NodeData(obj,bestbound))
        push!(bbdata2, NodeData_tst(time()-a,solution,obj,bestbound))
        check[i] = 1
    end
end
##############################################
     # build model ``m`` up here

times = [15,20, 30, 80, Inf] # Times that a new solution must be  ; leave Inf at the end
global check = zeros(length(times)) # Whether a 

bbdata2 = NodeData_tst[];
 

addinfocallback(m, infocallback, when = :Intermediate)
################################################

a = time()

@time status = solve(m)

# tmp_x = getvalue(x)
# find(tmp_x .== 1.)
# function tmpfun(x)
#        (x .== 1.)
# end
# mapslices(tmpfun, solutionvalues,1)



# global blau = zeros(5)
# function tmp()
#     blau[1] = 2
# end
# tmp()
# blau



##################################################################



# Optimal solution of relaxed problem will be (1.5, 2.0)

# We now define our callback function that takes one argument,
# the callback handle. Note that we can access m, x, and y because
# this function is defined inside the same scope
# function myheuristic(cb)
#     x_val = getvalue(x)
#     y_val = getvalue(y)
#     println("In callback function, x=$x_val, y=$y_val")

#     setsolutionvalue(cb, x, floor(x_val))
#     # Leave y undefined - solver should handle as it sees fit. In the case
#     # of Gurobi it will try to figure out what it should be.
#     addsolution(cb)

#     # Submit a second solution
#     setsolutionvalue(cb, x, ceil(x_val))
#     addsolution(cb)
# end  # End of callback function




# # Tell JuMP/Gurobi to use our callback function
# addheuristiccallback(m, myheuristic)

# function solution2matrix(x)
#     sizes = JuMP.size(x)
#     tmp = zeros(sizes)
#     for q in 1:sizes[1] , j in 1:sizes[2]
#       tmp[q,j] = x[q,j]
#     end
#     return tmp
# end



# type NodeData
#     time::Float64  # in seconds since the epoch
#     node::Int
#     obj::Float64
#     bestbound::Float64
# end

# # build model ``m`` up here

# bbdata = NodeData[]

# function infocallback(cb)
#     node      = MathProgBase.cbgetexplorednodes(cb)
#     obj       = MathProgBase.cbgetobj(cb)
#     bestbound = MathProgBase.cbgetbestbound(cb)
#     push!(bbdata, NodeData(time(),node,obj,bestbound))
# end
# addinfocallback(m, infocallback, when = :Intermediate)

# solve(m)

# # Save results to file for analysis later
# open("bbtrack.csv","w") do fp
#     println(fp, "time,node,obj,bestbound")
#     for bb in bbdata
#         println(fp, bb.time, ",", bb.node, ",",
#                     bb.obj, ",", bb.bestbound)
#     end
# end



# # Print our final solution
# println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")

