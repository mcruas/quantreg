# for max_K in [4,5,2,6,1]
#   include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/experimentos_grupo.jl")
# end
# include(homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/RegressãoQuantílica_STREET/experimentos_grupo.jl")

#######################################################
# Este script gera gráficos para o relatório, mostrando como uma variação
# em x_t afeta o quantil do período seguinte.

################# Gráficos linear #####################
vetor_x = [-1 0 1 2 3 10] ; n = 100
# X_tau = 500 # Qual o valor de X utilizado para a previsão
# limx = 1500; limy = 2000; # para ENA sudeste
limx = 55; limy = 55; # para icaraizinho
# limx = NaN; limy = NaN # quando não se sabe

using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx, JLD, RCall #, Distributions



usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
# max_K = 4
cd(pasta_trabalho)

#####################################################
############# Carregar dados em R ########################
pwd();
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
# include(pwd()*"/RegressãoQuantílica_STREET/par-multi.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
n = 100;
# @rput n pasta_trabalho;
# R"
# serie = arima.sim(n = n, list(ar = c(0.9,-0.2), ma = NULL),
#           sd = sqrt(0.1796))
# setwd(pasta_trabalho)
# source('R/biblioteca-funcoes-npquantile.R') # Carrega dados por script do R
# # source('R/biblioteca-funcoes-npquantile.R')
# # ipak(c('readxl', 'dplyr', 'lattice'))
# # dados <- read_excel(path = 'Dados Climaticos/Solar-tubarao/tubarao solar.xlsx')[,1:6]
# # dados_filtrados <- dados %>% select(yt0, yt1, hora, mes) %>% as.matrix
# # boxplot(yt0 ~ hora, data = dados.filtrados, col = 2)
# 
#
# @rget serie dados_filtrados
# serie = dados_filtrados[:,1].values;
################################################################

############# Carregar dados ###################################

# para icaraizinho
serie = readcsv("Dados Climaticos/icaraizinho.csv"); nomeserie = "icaraizinho"; x_new = 30
# serie = readcsv("Dados Climaticos/Enas/enas sudeste 1931-2014.csv"); nomeserie = "enasudeste"; serie = serie ./ 100; x_new = 500 # ENA Sudeste
# vetor_X_tau = collect(200:50:700) # para ENA sudeste


# serie = readcsv("Dados Climaticos/Solar-tabocas/tabocas.csv"); nomeserie = "solartabocas"; x_new = 0.4 # Para Tabocas


n_tmp = length(serie); # keeps size of series before cutting them
### Inicialização

# Does adaptation for a AR(2) model

# X = [serie[2:n_tmp-1] serie[1:n_tmp-2]]; # AR(2)

# X = Array{Float64,2} # AR(1)
tmp = serie[2:n_tmp-1];   # this is done so that X has two dimensions, instead of 1.
X =zeros(Float64,length(tmp),1);
X[:,1] = tmp;
x = X[:,1]

y = serie[3:n_tmp];
n = length(y);55
T = 1:n;

# X_lags = lagmatrix(serie,0:2)
# Alphas = collect(0.005:0.005:0.995);


# Y = [ y * ones(1,n_cenarios);
#       zeros(max_sim,n_cenarios)];

limx = isnan(limx) ? maximum(serie) : limx
limy = isnan(limy) ? maximum(serie) : limy

# plot(serie)


################################################################################################
################################################################################################
# Experimentos com o seleção inteira
X_lags = lagmatrix(serie,0:12)
Alphas = [0.05,0.1,0.25,0.5,0.75,0.9,0.95]
Alphas = vcat(collect(0.005:0.005:0.05), collect(0.1:0.05:0.9), collect(0.95:0.005:0.995))
TimeLimit = 600
vetor_TimeLimit = [200, 600, 1800]
vetor_Grupos = [1,2,3,10]

vetor_max_K = [1,5,6] 
# vetor_TimeLimit = [200, 600, 1800, 5400]


# vetor_Grupos = [1,2,3,10]
max_K = 2
Grupos = 37
# gr()
Guarda_Registros = Registros[];
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");

for max_K in vetor_max_K
    results = zeros(length(vetor_TimeLimit)*(length(vetor_Grupos)*2+1), 3); i = 1 # i keeps the row to write on results
    pasta_escrita = "RegressãoQuantílica_STREET/BetasMIP/K$max_K"

      ## MIP convencional       
    for TimeLimit in vetor_TimeLimit  
        betas0, betas, obj, status, tempo = rq_par_mip(X_lags[:,1], X_lags[:, 2:end], 
                                    Alphas; non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00)
        # betas_tmp = betas
        R"$betas[which(abs(apply($betas,2,order)) > $max_K)] = 0"
        writecsv("$pasta_escrita/Betas $TimeLimit GAll.csv", [betas0 ; betas])
        push!(Guarda_Registros ,Registros(betas0, betas, tempo, "$status", obj, "G$(length(Alphas))", max_K));
        results[i,:] = [length(Alphas), TimeLimit, obj] ; i += 1
        titulo = "GAll $TimeLimit - $obj"
        nome_arquivo = "$pasta_escrita/Heatmap $TimeLimit GAll.png"
        heatmap(betas .!= 0, title = titulo, xaxis = "Índice Alpha", yaxis = "Lag", legend = false)
        betas_diff = (betas .!= 0)+0
        R"
        tmp_betas = $betas_diff
        colnames(tmp_betas) = $Alphas
        rownames(tmp_betas) = 1:12
        png($nome_arquivo)
        heatmap(tmp_betas, Rowv=NA, Colv=NA, col = c(8,4) , scale='column', margins=c(5,10), xlab = 'Probability', ylab = 'Lags', main = $titulo)
        dev.off()
        "
        save("$pasta_escrita/Guarda_Registros.jld", "Guarda_Registros", Guarda_Registros)

        # savefig("$nome_arquivo")
        # close(fig)
    end
      ##  MIP com Grupos
    for TimeLimit in vetor_TimeLimit, Grupos in vetor_Grupos  
        betas0, betas, obj, status, tempo = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
            non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = Grupos)
        writecsv("$pasta_escrita/Betas $TimeLimit G$Grupos.csv", [betas0 ; betas])
      # betas_tmp = betas
        @rput betas
        R"betas[which(apply(abs(betas),2,order) <= $max_K)] = 0"
        @rget betas
        push!(Guarda_Registros ,Registros(betas0, betas, tempo, "$status", obj, "G$(Grupos)", max_K));
        results[i,:] = [length(Alphas), TimeLimit, obj] ; i += 1
        titulo = "G$Grupos $TimeLimit - $obj"
        nome_arquivo = "$pasta_escrita/Heatmap $TimeLimit G$Grupos.png"
        heatmap(betas .!= 0, title = titulo, xaxis = "Índice Alpha", yaxis = "Lag", legend = false)
        betas_diff = (betas .!= 0)+0
        R"
        tmp_betas = $betas_diff
        colnames(tmp_betas) = $Alphas
        rownames(tmp_betas) = 1:12
        png($nome_arquivo)
        heatmap(tmp_betas, Rowv=NA, Colv=NA, col = c(8,4) , scale='column', margins=c(5,10), xlab = 'Probability', ylab = 'Lags', main = $titulo)
        dev.off()
        "
        save("$pasta_escrita/Guarda_Registros.jld", "Guarda_Registros", Guarda_Registros)
     end

      ##  MIP com rampa
    for TimeLimit in vetor_TimeLimit, Grupos in vetor_Grupos  
         betas0, betas, obj, status, tempo = rq_par_mip_grupos_rampa(X_lags[:,1], X_lags[:, 2:end], Alphas;
            non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = Grupos)
        writecsv("$pasta_escrita/Betas $TimeLimit R$Grupos.csv", [betas0 ; betas])
      # betas_tmp = betas
        R"betas1 = $betas
          betas1[which(apply(abs(betas1),2,order) > $max_K)] = 0
          betas1"
        @rget betas1  
        push!(Guarda_Registros ,Registros(betas0, betas, tempo, "$status", obj, "R$(Grupos)", max_K));
        results[i,:] = [length(Alphas), TimeLimit, obj] ; i += 1
        titulo = "R$Grupos $TimeLimit - $obj"
        nome_arquivo = "$pasta_escrita/Heatmap $TimeLimit R$Grupos.png"
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
                save("$pasta_escrita/Guarda_Registros.jld", "Guarda_Registros", Guarda_Registros)

    end




    writecsv("$pasta_escrita/results.csv", results)
end




##### Para recarregar arquivos:
Auxiliar = load("$pasta_escrita/Guarda_Registros.jld") 
Guarda_Registros = Auxiliar["Guarda_Registros"]
# 
for i in 1:length(Guarda_Registros)
   println("$i:  $(Guarda_Registros[i].tempo)  $(Guarda_Registros[i].grupos)   $(Guarda_Registros[i].K)")
   end
# 
betas = Guarda_Registros[36].betas


# for Grupos in vetor_Grupos
#     betas0, betas, obj, status = rq_par_mip(X_lags[:,1], X_lags[:, 2:end], Alphas; non_cross = true, max_K = max_K, TimeLimit = 36000, MIPGap = 0.00)
    
#     betas0, betas, obj, status = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
#        non_cross = true, max_K = max_K, TimeLimit = 36000, MIPGap = 0.00, Grupos = Grupos)
# end


# y=X_lags[:,1]; X=X_lags[:, 2:end]; non_cross = true; max_K = max_K; TimeLimit = 200; MIPGap = 0.00; Grupos = Grupos

# TimeLimit = 150
# include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
# global times = [80,120, 150,190 , Inf] # Times that a new solution must be  ; leave Inf at the end
# global check = zeros(length(times)) # Whether a 
#  betas0_rampa, betas_rampa, obj_rampa, status_rampa, tempo_rampa = rq_par_mip_grupos_rampa(X_lags[:,1], X_lags[:, 2:end], Alphas;
#    non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = Grupos)
# 
    # betas0_grupos, betas_grupos, obj_grupos, status_grupos = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
    #    non_cross = true, max_K = max_K, TimeLimit = 900, MIPGap = 0.00, Grupos = Grupos)


# heatmap(Alphas, betas .== 0)
# savefig("Coeficientes MIP apos 6000s.pdf")

# status == :Blau? 2 : status == :Optimal? 1 : 0


# betas0_rampa, betas_rampa, obj_rampa, status_rampa, tempo_rampa
# push!(Guarda_Registros ,Registros(betas0_rampa, betas_rampa, tempo_rampa, status_rampa, obj_rampa, Grupos, max_K));