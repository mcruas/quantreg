using JuMP, DataFrames, Plots, Interpolations, LaTeXStrings, Dierckx, JLD, RCall #, Distributions


usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)

include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");
n = 100;





########### PARAMETROS ##########################
window_size = 24*30





#################################################

# para icaraizinho
WindData = readcsv("Dados Climaticos/Dados-kaggle/DadosKg2014Adj.csv")[2:(end-1),:];
# nomeserie = "icaraizinho"; x_new = 30

dates_input = map(x->string(x), WindData[:,1])

year = map(x->parse(Int64,x[1:4]),dates_input) # extracts a given position from a whole vector of strings
month= map(x->parse(Int64,x[5:6]),dates_input)
day = map(x->parse(Int64,x[7:8]),dates_input)
hour = map(x->parse(Int64,x[9:10]),dates_input)
# DateTime(string(WindData[:,1]),"yyyymmdd").
n = size(WindData,1)




################################## MAIN ################################


rq_par_mip(y::Array{Float64}, X::Array{Float64,2}, Alphas; non_cross = true, max_K = NaN, TimeLimit = NaN, MIPGap = NaN)
rq_par_mip_grupos(y::Array{Float64}, X::Array{Float64,2}, Alphas; non_cross = true, max_K = NaN, TimeLimit = NaN, MIPGap = NaN, Grupos = NaN)

lags = collect()
oldest_lag = 

lagmatrix(serie,0:12)


##################################

# Escolho dias para o BacktestData
day_ini = 20130117;
day_end = 20130118;

day_ini_ano = (string(day_ini))[1:4]
day_ini_mes =(string(day_ini))[5:6]
day_ini_dia = (string(day_ini))[7:8]
DiaIni = Date(Date(string(day_ini_ano,"-",day_ini_mes,"-", day_ini_dia)));

day_end_ano = (string(day_end))[1:4]
day_end_mes =(string(day_end))[5:6]
day_end_dia = (string(day_end))[7:8]
DiaEnd = Date(Date(string(day_end_ano,"-",day_end_mes,"-", day_end_dia)));

DiasAux = collect(DiaIni:Dates.Day(1):DiaEnd) ;
Dias = Array(Int64 , size(DiasAux,1)) ;
for i=1:size(Dias,1)
  Dias[i]=parse(Int64,string(string(DiasAux[i])[1:4],string(DiasAux[i])[6:7],string(DiasAux[i])[9:10]))
end

DiasTotais = size(Dias,1)
nUnitComits = convert(Int64,DiasTotais)


#################################################

for dia = 1:size(Dias,1)

 tic()

 print("Otimização para o dia: ");println(Dias[dia])

 day = Dias[dia];
 day_ano = (string(day))[1:4]
 day_mes =(string(day))[5:6]
 day_dia = (string(day))[7:8]
 DiaAtual = DateTime(Date(string(day_ano,"-",day_mes,"-", day_dia)));

 #Encontro indice para o dia
 Idday = find(WindData[:,1].==day*100);

 # Encontrando Dias recentes para o U.S. ##########
 AuxD1 = Date(DiaAtual - Dates.Day(nDiasPassados));
 AuxDN = Date(DiaAtual - Dates.Day(1));

 #IDiasPassados = collect((DiaAtual - Dates.Day(nDiasPassados) ): ( DiaAtual - Dates.Day(1)))
 AuxCmp1 = parse(Int64,string(string(AuxD1)[1:4],string(AuxD1)[6:7],string(AuxD1)[9:10]))
 AuxCmp2 = parse(Int64,string(string(AuxDN)[1:4],string(AuxDN)[6:7],string(AuxDN)[9:10]))
 IDiasRecentes=find((WindData[:,1] .>=   AuxCmp1*100)  & (WindData[:,1] .<   day*100));

 if mod(size(IDiasRecentes,1),24) != 0 ; error("Problema de datas - intradiário") ; end;

 #### Encontrando Dias de anos anteriores para o U.S. #########

 AuxD1 = Date(DiaAtual - Dates.Year(1) - Dates.Day(nDiasAnoAnterior));
 AuxDN = Date(DiaAtual - Dates.Year(1) + Dates.Day(nDiasAnoAnterior));

 AuxCmp1 = parse(Int64,string(string(AuxD1)[1:4],string(AuxD1)[6:7],string(AuxD1)[9:10]))
 AuxCmp2 = parse(Int64,string(string(AuxDN)[1:4],string(AuxDN)[6:7],string(AuxDN)[9:10]))
 IDiasAnoPassado =find((WindData[:,1] .>= AuxCmp1*100)  & (WindData[:,1] .< AuxCmp2*100));

 if mod(size(IDiasAnoPassado,1),24) != 0 ; error("Problema de datas - intradiário") ; end;

 ###### Número de dias do Uncertainty set ###########

 dDOM.nHist = trunc(size(IDiasRecentes,1)/24 + size(IDiasAnoPassado,1)/24,0)
 dDOM.w = zeros(dDH.nW,dDH.nT,dDOM.nHist);

 IDiasAux = [IDiasAnoPassado ; IDiasRecentes]
 dDOM.d = Array(Float64,dDH.nB,dDH.nT,dDOM.nHist);

 for k = 1:dDOM.nHist
   IVento = ((k-1) * 24 + 1):((k-1) * 24 + 24)

   dDOM.w[:,:,k] = round(0+(WindData[IDiasAux[IVento],2:end])',2)
   dDOM.d[:,:,k] = dDE.d̂ + 0

end;

 dDE.ŵ = round(mean(dDOM.w[:,:,convert(Int64,trunc(size(IDiasAnoPassado,1)/24,0)+1):end],3),2)[:,:,1]


 ############ JANELA MOVEL ###########
 