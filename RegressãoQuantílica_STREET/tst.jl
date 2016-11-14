X_reduc = [85	96
118	105
97	131
117	133
76	149
99	110
84	132
113	146
90	160]

y = [1
1
1
1
1
0
0
0
0]


### MODEL ###

m = Model(solver = solvfunc)
# m = Model(solver = MosekSolver())

A = 1:3
B = 1:4


@defVar(m, a[A,B])

# Objective
@setObjective(m, Max, sum{a[i,j], i = A, j = B})
# Constraints

tmp = 1:3
typeof(tmp)

i=1

@addConstraints(m, begin
                  #class0[i = nM0, k = K],  sum{X_reduc[M0[i],f] * p[f,k,M0[i]], f = F} + q[k,M0[i]]  <= - δ
                  class0[i = nM0], q[1,M0[i]] <= 3
                  #class1[i = M0, k = K, ]
                  #igualwp3[j = P], -a[j] * M <= w[j]
                  #igualwp4[j = P], w[j] <= a[j] * M
                end)
