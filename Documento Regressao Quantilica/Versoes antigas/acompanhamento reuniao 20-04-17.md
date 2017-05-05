---
title: "Acompanhamento Reunião 20-04-17"
author: "Marcelo Ruas"
date: "April 20, 2017"
output: pdf_document
---

- 1. Introduzir $L_\lambda$ antes do modelo, na equação 3.14

- 2. Colocar equação 3.22 na 3.23. Eliminar e substituir por 
$$
\lambda^*_K = \underset{\lambda}{\text{arg min}} \left\lbrace \left.  obj_{\lambda}^{*} \quad  \right| \, \| \beta^*_\lambda \|_0 = K \right\rbrace,
$$
Esta nova equação está como equação 3.19.
Defini a função objetivo como $obj^*_\lambda$ ao invés de $\hat{\sigma}_\lambda^*$.

- 2.1. Definir a distância como um problema de otimização:
\begin{eqnarray}
d(\beta^*_{MILP(K)}, \beta^*_{\lambda^*_K}) =	1 - \max_{0\leq\delta_{ij}\leq1} & \sum\sum_{j} \delta_{ij} |\rho_{ij}| \label{eq:metricad0} \\
\text{s.t.} & \sum_{j}\delta_{ij}=1 & \forall i\in L_{K}^{MILP},\\
& \sum_{i}\delta_{ij}=1 & \forall j\in L_{K}^{LASSO},\\
& \delta_{i,j} = 0, & \forall i \in \overline{L}_{K}^{MILP}, \forall j \in \{1,\dots,P\},\\
& \delta_{i,j} = 0, & \forall j \in \overline{L}_{K}^{LASSO}, \forall i \in \{1,\dots,P\},\label{eq:metricad4}
\end{eqnarray}

- 3. Indicar o modelo selecionado pelo LASSO e pelo MIP, indicando graficamente. Feito na figura 3.2.
- 3.1. Incluir seção sobre o modelo MILP com divisão em grupos. O problema ficou definido como: 
\begin{eqnarray}
 \underset{\beta_{0\alpha},\beta_\alpha,z_{p \alpha} \varepsilon_{t \alpha}^{+},\varepsilon_{t \alpha}^{-}}{\text{min}} & \sum_{\alpha \in A} \sum_{t\in T}\left(\alpha\varepsilon_{t \alpha}^{+}+(1-\alpha)\varepsilon_{t\alpha}^{-}\right) \label{eq:mip0} \\
\mbox{s.t } & \varepsilon_{t \alpha}^{+}-\varepsilon_{t \alpha}^{-}=y_{t}-\beta_{0 \alpha}-\sum_{p=1}^{P}\beta_{p \alpha}x_{t,p},& \qquad\forall t \in T ,\forall \alpha \in A, \label{eq:mip1}\\
& \varepsilon_{t \alpha}^{+},\varepsilon_{t \alpha}^{-}\geq0,&\qquad\forall t \in T ,\forall \alpha \in A, \label{eq:mip2}\\
& - M z_{p \alpha} \leq \beta_{p \alpha} \leq M z_{p \alpha},&\qquad \forall \alpha \in A, \forall p\in\{1,\dots,P\}, \label{eq:mip3}\\
& \sum_{p=1}^P z_{p \alpha} \leq K, & \qquad \forall \alpha \in A, \label{eq:mip4}\\
& z_{p \alpha} \in \{0,1\},&\qquad \forall \alpha \in A, \forall p\in\{1,\dots,P\}, \label{eq:mip5}\\
& \beta_{0\alpha} + \beta_{\alpha}^T x_{t} \leq \beta_{0\alpha'} + \beta_{\alpha'}^T x_{t}, & \qquad \forall t \in T, \forall (\alpha, \alpha') \in A \times A,  \alpha < \alpha',\nonumber\\ \label{eq:mip6}
\end{eqnarray}
Ver seção 3.1

- 4. Gerar dados para testar poder computacional do MIP vs LASSO. 1000 regressores, 200 dados e variando $K$

- 5. Pesquisar literatura sobre avaliar a distribuição completa em séries temporais.
    - Observar se o tema de determinar $F_{Y_t|X_t}$ (distribuição) está presente
  
  
- 6. Pensar na seleção do $\lambda_2$ como uma minimização do erro da distribuição incondicional. 


