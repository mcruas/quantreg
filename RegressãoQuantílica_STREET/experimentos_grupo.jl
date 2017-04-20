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

using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx #, Distributions



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
TimeLimit = 200
vetor_TimeLimit = [200]
vetor_Grupos = [1]

vetor_max_K = [4,5,6,2,1] 
vetor_TimeLimit = [200, 600, 1800, 5400]
vetor_Grupos = [1,2,3,10]
pyplot()
for max_K in vetor_max_K
    results = zeros(length(vetor_TimeLimit)*(length(vetor_Grupos)+1), 3); i = 1 # i keeps the row to write on results
    pasta_escrita = "RegressãoQuantílica_STREET/BetasMIP/K$max_K"
    for TimeLimit in vetor_TimeLimit  
        betas0, betas, obj, status = rq_par_mip(X_lags[:,1], X_lags[:, 2:end], Alphas; non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00)
        writecsv("$pasta_escrita/Betas GAll $TimeLimit.csv", [betas0 ; betas])
        results[i,:] = [length(Alphas), TimeLimit, obj] ; i += 1
        heatmap(betas .== 0, title = "GAll $TimeLimit - $obj")
        savefig("$pasta_escrita/Heatmap GAll $TimeLimit.png")
        # close(fig)
    end

    for TimeLimit in vetor_TimeLimit, Grupos in vetor_Grupos  
        betas0, betas, obj = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
            non_cross = true, max_K = max_K, TimeLimit = TimeLimit, MIPGap = 0.00, Grupos = Grupos)
        writecsv("$pasta_escrita/Betas G$Grupos $TimeLimit.csv", [betas0 ; betas])
        results[i,:] = [Grupos, TimeLimit, obj] ; i += 1
        heatmap(betas .== 0, title = "GAll $TimeLimit - $obj")
        savefig("$pasta_escrita/Heatmap G$Grupos $TimeLimit.png")
    end

    
    writecsv("$pasta_escrita/results.csv", results)
end

max_K = 3

Grupos = 2
for Grupos in vetor_Grupos
    betas0, betas, obj, status = rq_par_mip(X_lags[:,1], X_lags[:, 2:end], Alphas; non_cross = true, max_K = max_K, TimeLimit = 36000, MIPGap = 0.00)
    
    betas0, betas, obj, status = rq_par_mip_grupos(X_lags[:,1], X_lags[:, 2:end], Alphas;
       non_cross = true, max_K = max_K, TimeLimit = 36000, MIPGap = 0.00, Grupos = Grupos)
end


    betas0, betas, obj, status = rq_par_mip_grupos_rampa(X_lags[:,1], X_lags[:, 2:end], Alphas;
       non_cross = true, max_K = max_K, TimeLimit = 36000, MIPGap = 0.00, Grupos = Grupos)

# heatmap(Alphas, betas .== 0)
# savefig("Coeficientes MIP apos 6000s.pdf")

# status == :Blau? 2 : status == :Optimal? 1 : 0