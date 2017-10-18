## 
## OMP_NUM_THREADS=5 julia ~/Dropbox/Pesquisa\ Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/Procedimentos_experimentos/experimentos_grupos_3.jl
## Experimentos em que as simulações são executadas até atingida a convergência ou 18000s (5 horas)

using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx, JLD, RCall #, Distributions


usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
n = 100;

# para icaraizinho
serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"; x_new = 30

n_tmp = length(serie); # keeps size of series before cutting them

### Inicialização dos dados

tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X =zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]
y = serie[3:n_tmp];
n = length(y);55
T = 1:n;


################################################################################################
################################################################################################
# Experimentos com o seleção inteira
X_lags = lagmatrix(serie,0:12)
Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
TimeLimit = 157
vetor_TimeLimit = [200, 600, 1800]


# Recolocar os grupos
vetor_Grupos = [1,2,3,10, length(Alphas)]
vetor_max_K = [2] 
 vetor_TimeLimit = [18000] # Máximo de dois dias de simulação


# TEMP
vetor_Grupos = [10]
vetor_max_K = [2] 
 vetor_TimeLimit = [157] # Máximo de dois dias de simulação



max_K = 2
Grupos = 36
# gr()
# Guarda_Registros = Registros[];
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");

for max_K in vetor_max_K
    results = zeros(length(vetor_TimeLimit)*(length(vetor_Grupos)*2+1), 3); i = 1 # i keeps the row to write on results
    pasta_escrita = "RegressãoQuantílica_STREET/BetasMIP/K$max_K"

      ##  MIP com Grupos
    for TimeLimit in vetor_TimeLimit, Grupos in vetor_Grupos  
        println("******************************************")
        println("$max_K - G$Grupos")
        println("******************************************")
        betas0, betas, obj, status, tempo = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
            non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = Grupos)
        writecsv("$pasta_escrita/Betas $TimeLimit G$Grupos - $tempo.csv", [betas0 ; betas])
      # betas_tmp = betas
         # Limpa valores quase zero, restando apenas os K maiores
        for coluna in 1:size(betas)[2]
              betas[sortperm(abs(betas[:,coluna]), rev = true)[max_K+1:end], coluna] = 0
        end
        # push!(Guarda_Registros ,Registros(betas0, betas, tempo, "$status", obj, "G$(Grupos)", max_K));
        results[i,:] = [length(Alphas), TimeLimit, obj] ; i += 1
        titulo = "G$Grupos $TimeLimit - $obj"
        nome_arquivo = "$pasta_escrita/Heatmap - $(round(TimeLimit)) - G$grupos - $obj - $tempo.png"
        # heatmap(betas .!= 0, title = titulo, xaxis = "Índice Alpha", yaxis = "Lag", legend = false)
        betas_diff = (betas .!= 0)+0
        R"
        tmp_betas = $betas_diff
        colnames(tmp_betas) = $Alphas
        rownames(tmp_betas) = 1:12
        png($nome_arquivo)
        heatmap(tmp_betas, Rowv=NA, Colv=NA, col = c(8,4) , scale='column', margins=c(5,10), xlab = 'Probability', ylab = 'Lags', main = $titulo)
        dev.off()
        "
        # save("$pasta_escrita/Guarda_Registros.jld", "Guarda_Registros", Guarda_Registros)
     end

      ##  MIP com rampa
    for TimeLimit in vetor_TimeLimit, Grupos in vetor_Grupos[1:end-1]  
        println("******************************************")
        println("$max_K - R$Grupos")
        println("******************************************")

         betas0, betas, obj, status, tempo = rq_par_mip_grupos_rampa(X_lags[:,1], X_lags[:, 2:end], Alphas;
            non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = Grupos)
        writecsv("$pasta_escrita/Betas $TimeLimit R$Grupos.csv - $tempo", [betas0 ; betas])
      # betas_tmp = betas
        # Limpa valores quase zero, restando apenas os K maiores
        for coluna in 1:size(betas)[2]
              betas[sortperm(abs(betas[:,coluna]), rev = true)[max_K+1:end], coluna] = 0
        end
        # push!(Guarda_Registros ,Registros(betas0, betas, tempo, "$status", obj, "R$(Grupos)", max_K));
        results[i,:] = [length(Alphas), TimeLimit, obj] ; i += 1
        titulo = "R$Grupos $TimeLimit - $obj"
        nome_arquivo = "$pasta_escrita/Heatmap - $(round(TimeLimit)) - R$grupos - $obj - $tempo.png"
        # heatmap(betas .!= 0, title = titulo, xaxis = "Índice Alpha", yaxis = "Lag", legend = false)
        betas_diff = (betas .!= 0)+0
        R"
        tmp_betas = $betas_diff
        colnames(tmp_betas) = $Alphas
        rownames(tmp_betas) = 1:12
        png($nome_arquivo)
        heatmap(tmp_betas, Rowv=NA, Colv=NA, col = c(8,4) , scale='column', margins=c(5,10), xlab = 'Probability', ylab = 'Lags', main = $titulo)
        dev.off()
        "
                # save("$pasta_escrita/Guarda_Registros.jld", "Guarda_Registros", Guarda_Registros)

    end
    # writecsv("$pasta_escrita/results.csv", results)
end
