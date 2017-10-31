# Gerar dados p testar poder computacional do MIP vs LASSO, 1000 regressores , 200 dados e variando K
# Ver base de dados do trabalho com algoritmo genético
# R:  install.packages(c("glmnet", "mvtnorm", "clusterGeneration"))

using RCall, Plots
T = 40; N = 30; q = 2; rho = 0
sigma_e = 0.3; iid = true


############################ Gerar dados ########################
@rput T N q sigma_e rho iid;
R"
source('R/data-generation-mipvslasso.R');
dados = dgp(T,N,q, iid, rho = rho, sigma_e = sigma_e)
";
@rget dados;
X = dados[:X];
y = dados[:Y];
sigma = dados[:sigma];
beta_true = dados[:beta];
rho = dados[:rho];
if (rho > 0) 
  @show mean(y);
  # R"jb.norm.test($y)"
  @show E_yt = (q)*0.3/(1-rho);
  @show var(y);
  @show V_yt = (beta_true' * sigma * beta_true + sigma_e)[1] / (1-rho^2);
else 
  print(mean(y));
  # jb.norm.test(y)
  @show E_yt = (q)*0.3;
  @show var(y);
  @show V_yt = (beta_true' * sigma * beta_true + sigma_e)[1];
end


######################## Procedimentos para seleção automática ################
# Estimar a distribuição condicional com grupos
Alphas = collect(0.05:0.05:0.95)

usesolver = "gurobi"    # Escolher entre os valores 'mosek' ou 'gurobi'
# cd("/home/marcelo/Dropbox/Pesquisa Doutorado/Paper NPQuantile/RegressãoQuantílica_STREET")
pasta_trabalho = homedir()*"/Dropbox/Pesquisa Doutorado/Paper-NPQuantile/"
cd(pasta_trabalho)
include(pwd()*"/RegressãoQuantílica_STREET/funcoes_npqar.jl");
include(pwd()*"/RegressãoQuantílica_STREET/npar-multi-funcoes.jl");


# @show resultados

R"
library(rqPen)
x <- matrix(rnorm($T * $N),ncol=$N)
y <- 1 + x[,1] - 1*x[,2] + 1*x[,3] + rnorm($T)
cv_model <- cv.rq.pen(x,y,tau=0.3,lambda=NULL,weights=NULL,penalty='LASSO', intercept=TRUE)
"
@rget x y cv_model

resultados = rq_par_mip(y,X,Alphas; max_K = 4, MIPGap = 0.2)


#  cv_model$models[[91]]
