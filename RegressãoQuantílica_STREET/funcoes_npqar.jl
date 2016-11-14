# Quando um array é variável de decisão, usar essa função já retorna um vetor
function getValueVector(a)
    x = getValue(a)
    n = size(x)[1]
    tmp = zeros(n)
    for i = 1:n
        tmp[i] = x[i]
    end
    return(tmp)
end


# Quando um array é variável de decisão, usar essa função já retorna um vetor
function getValueVectorCoco(a, dim1)
    x = getValue(a)
    tmp = zeros(dim1)
    for i = 1:dim1
        tmp[i] = x[i]
    end
    return(tmp)
end


# Quando um array é variável de decisão, usar essa função já retorna um vetor
function getValueMatrix(a)
    x = getValue(a)
    tamanho = size(x)
    tmp = zeros(tamanho)
    for i = 1:tamanho[1], j = 1:tamanho[2]
        tmp[i,j] = x[i,j]
    end
    return(tmp)
end


# Quando um array é variável de decisão, usar essa função já retorna um vetor
# We need to inform both dimensions of matrix
# dim1 is the number of lines, dim2 of columns
function getValueMatrixCoco(a, dim1, dim2)
    x = getValue(a)
    tmp = zeros(dim1, dim2)
    for i = 1:dim1, j = 1:dim2
        tmp[i,j] = x[i,j]
    end
    return(tmp)
end


function printvector(a)
  n = size(a)[1]
  for i = 1:n
    print(i, ": ", a[i], "     ")
  end
end

# Função que recebe os valores previstos e os valores de y
#t

function avaliacao(previsoes, y_test)
  y0 = y_test .== 0
  y1 = y_test .== 1
  # Sobre as avaliações, é muito pior o ETII, em que um paciente com Cancêr
  # recebe o diagnóstico que não tem, do que o contrário.

  # Avaliação
   # Erros Tipo I: Paciente sem cancer (y=0) recebe o diagnóstico de que tem (y=1)
  erros_T1 = sum(previsoes[y0])/sum(y0)
  # Erros Tipo II: Paciente com cancer (y=1) recebe o diagnóstico de que não tem (y=0)
  erros_T2 = 1 - sum(previsoes[y1])/sum(y1)
  return erros_T1, erros_T2
end

function colMeans(x)
  p = size(x)[2]
  output = zeros(p)
  for i in 1:p
    output[i] = mean(x[:,i])
  end
  return output
end

# apply function f over the columns of x
function colApply(x, f)
  p = size(x)[2]
  output = zeros(p)
  for i in 1:p
    output[i] = f(x[:,i])
  end
  return output
end

# apply function f over the rows of x
function rowApply(x, f)
  p = size(x)[1]
  output = zeros(p)
  for i in 1:p
    output[i] = f(x[i,:])
  end
  return output
end

# apply function f over all elements of x
function elemApply(x, f)
  p = size(x)
  output = zeros(p)
  for i in 1:p[1], j in 1:p[2]
    output[i,j] = f(x[i,j])
  end
  return output
end



function dict2vec(a)
  tmp = zeros(P)
  for i = P
    tmp[i] = a[i]
  end
end

# receives a vector x and a number of lags do make de matrix
function lagmatrix(x, lags)
  n = length(x)
  I = (lags+1):(n)
  x_estim = zeros(n-lags, lags)
  for i in 1:lags
    x_estim[:,i] = x[I - i]
  end
  return x_estim
end


function normalize_matrix(x)
  n = size(x)
  col_mean = colApply(x, mean)
  col_sd = colApply(x, std)
  output = zeros(n)
  for i in 1:n[1]
    output[i,:] = (x[i,:] - col_mean') ./ col_sd'
  end
  return output, col_mean, col_sd
end




if usesolver == "mosek"
  using Mosek
  solvfunc = MosekSolver()
elseif usesolver == "gurobi"
  using Gurobi
  solvfunc = GurobiSolver()
else # usar o GLPK
  using GLPKMathProgInterface
  solvfunc = GLPKSolverMIP()
end
