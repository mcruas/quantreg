@everywhere function thuener2(N::Int64)
  const N_quebra = 10000
  N_iter = convert(Int64,ceil(N / N_quebra))
  n_circulo = @parallel (+) for i= 1:N_iter
    sum((rand(N_quebra).^2 + rand(N_quebra).^2) .< 1)
    # rand()
  end
  n_amostra = N_iter * N_quebra
  pi_amostra = 4*n_circulo/n_amostra
end

@everywhere function thuener3(N::Int64)
  N_quebra = 10000
  N_iter = convert(Int64,ceil(N / N_quebra))
  n_circulo = @parallel (+) for i= 1:N_iter
    sum((rand(N_quebra).^2 + rand(N_quebra).^2) .< 1)
    # rand()
  end
  n_amostra = N_iter * N_quebra
  pi_amostra = 4*n_circulo/n_amostra
end

# addprocs(7)
N = 1000000000
tic(); @show thuener2(N); toc()
