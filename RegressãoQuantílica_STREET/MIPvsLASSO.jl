# Gerar dados p testar poder computacional do MIP vs LASSO, 1000 regressores , 200 dados e variando K
# Ver base de dados do trabalho com algoritmo genético
# R:  install.packages(c("glmnet", "mvtnorm", "clusterGeneration"))

using RCall, Plots
T = 100000; N = 20; q = 20


@rput T N  q
R"
source('R/data-generation-mipvslasso.R');
dados = dgp(T,N,q);
"
@rget dados
X = dados[:X]
y = dados[:Y]
beta_true = dados[:beta]

# gr()
# histogram(y)

# R"glm($y ~ $X)"

mean(y)