# Gerar dados p testar poder computacional do MIP vs LASSO, 1000 regressores , 200 dados e variando K
# Ver base de dados do trabalho com algoritmo genético
# R:  install.packages(c("glmnet", "mvtnorm", "clusterGeneration"))

using RCall, Plots
T = 1000; N = 20; q = 20; rho = 0.5
autoregressive = true ; sigma_e = 0.3; iid = false


############################ Gerar dados ########################
@rput T N q autoregressive sigma_e rho iid;
R"
source('R/data-generation-mipvslasso.R');
dados = dgp(T,N,q, iid, rho = rho, autoregressive =  autoregressive,sigma_e = sigma_e)
";
@rget dados E_yt V_yt;
X = dados[:X];
y = dados[:Y];
sigma = dados[:sigma];
beta_true = dados[:beta];
rho = dados[:rho];
if (autoregressive) 
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