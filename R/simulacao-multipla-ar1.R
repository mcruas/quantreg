library(forecast)
library(quantreg)
library(reshape2)
library(dplyr)
library(xtable)
library(stringr)
library(ggplot2)

# 1) tamanho da amostra, sugiro utilizarmos 4 tamanhos (20, 50, 100 e 1000
# observações), e 2) número de realizações: sugiro reportarmos não o resultado
# para uma série simulada mas sim a média para 10,000 (dez mil) simulações desse
# processo para cada parametrização.
# 
# 1 e 2 implicam em colocar dois for no codigo, um por fora, para considerar os 
# diferentes tamanhos de série e outro por dentro, que repete o experimento para
# cada parametrização 10000 vezes e acumula um resultado médio (não guarde o 
# resultado individual, vá acumulando os valores das razoes e tirando a média ao
# mesmo tempo).
# 
# Dê prioridade a isso do que à métrica.


# Simulate AR(1) process --------------------------------------------------
n.sim <- 1000
vector.n <- c(20, 50, 100, 1000)
vector.alpha <- c(0.05, 0.1, 0.5, 0.9, 0.95)
vector.phi <- c(0.25, 0.5, 0.7, 0.9)
vector.RSN <- c(0.01, 0.05, 0.1, 0.5, 1)
n.rows <- length(vector.alpha) * length(vector.phi) * length(vector.RSN) * length(vector.n)
n = 20; alpha = 0.05; i = 1; phi = 0.9; RSN = 0.1; phi0  = 1; serie = rep(0,n); iter = 1

results <- data.frame(n = rep(0,n.rows), alpha = rep(0,n.rows), phi = rep(0,n.rows), RSN = rep(0,n.rows), sigma2 = rep(0,n.rows),
                      phi.hat.mean = rep(0,n.rows), phi.hat.sd = rep(0,n.rows), 
                      sigma2.hat.mean = rep(0,n.rows),sigma2.hat.sd = rep(0,n.rows), 
                      phi0.hat.mean = rep(0,n.rows), phi0.hat.sd = rep(0,n.rows),  
                      beta.hat.mean = rep(0,n.rows), beta.hat.sd = rep(0,n.rows),
                      beta0.hat.mean = rep(0,n.rows), beta0.hat.sd = rep(0,n.rows),
                      RMSE.arima = rep(0,n.rows), RMSE.qr = rep(0, n.rows),
                      diffrmse.mean = rep(0,n.rows), diffrmse.sd = rep(0,n.rows), 
                      diff.intercept.mean = rep(0,n.rows), diff.intercept.sd = rep(0,n.rows),
                      diff.lag.mean = rep(0,n.rows), diff.lag.sd = rep(0,n.rows)
                      )
tmp.result <- data.frame(phi.hat = rep(0,n.sim), phi0.hat = rep(0,n.sim), sigma.e = rep(0,n.sim),
                         beta = rep(0,n.sim), beta0 = rep(0,n.sim), 
                         RMSE.arima = rep(0,n.sim), RMSE.qr = rep(0, n.sim), intercept = rep(0,n.sim))

for (n in vector.n) {
  for (alpha in vector.alpha) {
    for (phi in vector.phi) {
      for (RSN in vector.RSN) {
        set.seed(n + 100*alpha + 100*phi + 100*RSN)
        sigma.e = (phi0/(1-phi)) * RSN
        # y_t = phi0 + phi * y_t-1 + epsilon
        for (iter in 1:n.sim) {
          
          serie.gen <- rep(0,2*n); serie[1] = (phi0/(1-phi))
          for (j in 2:(2*n)) {
            serie.gen[j] <- serie.gen[j-1] * phi + phi0 + rnorm(n=1,sd = sqrt(sigma.e))
          }
          testing <- serie.gen[(n+1):(2*n)]
          serie <- serie.gen[1:n]
          true.quantile <- testing[-1] + qnorm(alpha) * sqrt(sigma.e)
          
          # Estimate AR(1) model with normal ARIMA ----------------------------------
          
          # ar(serie, order.max = 1, method = "ols")
          # fit1 <- Arima(serie, order = c(1,0,0), include.constant = TRUE)
          # phi0.hat <- coef(fit1)[2]*(1-phi)
          # sigma2.hat <- fit1$sigma2
          # phi.hat <- coef(fit1)[1]
          # 
          
          fit.tmp <- lm(formula = serie[2:n] ~ serie[1:(n-1)] + 1)
          phi.hat <- coef(fit.tmp)[2]
          phi0.hat <- coef(fit.tmp)[1]
          sigma2.hat <- var(fit.tmp$residuals)
          
          
          
          # lm(serie ~ lag(serie) + 1)
          
          prediction.arima <- (testing[-n] * phi.hat + phi0.hat + qnorm(alpha) * sqrt(sigma2.hat))
          # lines(prediction.arima, col = 2)
          
          fit2 <- rq(serie[2:n] ~ serie[1:(n-1)] + 1, tau = alpha)
          beta0.hat <- coef(fit2)[1]
          beta.hat <- coef(fit2)[2]
          
          prediction.qr <- (testing[-n] * beta.hat + beta0.hat)
          
          error.arima <- sqrt(mean((prediction.arima - true.quantile)^2))
          error.qr <- sqrt(mean((prediction.qr - true.quantile)^2))
          
          intercept.ar <- phi0.hat + qnorm(alpha) * sqrt(sigma2.hat)
          diff.intercept <- (intercept.ar) - beta0.hat
          
          tmp.result[iter, ] <- c(phi.hat, phi0.hat, sigma2.hat,
                                  beta.hat, beta0.hat, 
                                  error.arima, error.qr, intercept.ar)
          
        }
        
        # View(tmp.result)
                
        # plot(density(tmp.result[, 1]), main = str_c("N = ",n, ", Phi = ", phi,", RSN = ", RSN, 
        #                                                    " and Alpha = ", alpha))
        # lines(density(tmp.result[,"beta"]), col = 2)
        # abline(v = phi)
        
        # plot.ts(true.quantile)
        # write.table(serie, file = "R/serie_ar1.csv", row.names = FALSE, col.names = FALSE, sep = ",")
        
        # plot.ts(serie)
        
        # tmp.result[1,]
        
        diff.intercept <- tmp.result[,"intercept"] - tmp.result[,"beta0"]
        diff.lag <- tmp.result[,"phi.hat"] - tmp.result[,"beta"]
        diff.rmse.mean <- sum(tmp.result[,"RMSE.qr"])/sum(tmp.result[,"RMSE.arima"])
        diff.rmse.sd <- sd((tmp.result[,"RMSE.qr"] / tmp.result[,"RMSE.arima"]))
        results[i, ] <-   c(n, alpha, phi, RSN, sigma.e,
                            mean(tmp.result[,"phi.hat"]), sd(tmp.result[,"phi.hat"]),  
                            mean(tmp.result[,"sigma.e"]), sd(tmp.result[,"sigma.e"]),
                            mean(tmp.result[,"phi0.hat"]), sd(tmp.result[,"phi0.hat"]),
                            mean(tmp.result[,"beta"]), sd(tmp.result[,"beta"]),
                            mean(tmp.result[,"beta0"]), sd(tmp.result[,"beta0"]),
                            mean(tmp.result[,"RMSE.arima"]), mean(tmp.result[,"RMSE.qr"]),
                            diff.rmse.mean, diff.rmse.sd, 
                            mean(diff.intercept),  sd(diff.intercept),
                            mean(diff.lag), diff.rmse.sd)
        
        
        # plot(prediction.qr, prediction.arima, main = str_c("Phi = ", phi,", RSN = ", RSN, 
        #                                                    " and Alpha = ", alpha), 
        #                       # xlim = c(0, max(prediction.qr)), ylim = c(0,max(prediction.arima)),
        #      pch = 16 )
        #     abline(a = 0, b = 1, lty = 2, col = 2)
        i = i + 1
        # plot(prediction.qr[1:200],main = str_c("Phi = ", phi,", RSN = ", RSN, 
        #                                       " and Alpha = ", alpha))
        # points(prediction.arima[1:200],pch = 3, col = 2)
      }
    }
  }
}

# View results and prints to files ---------------------------------------------------------

# View(results)
write.csv(results, file = "results-multi-sim-tmp.csv", row.names = FALSE)
# names(results)
# pasta.destino <- "Documento Regressao Quantilica/Tabelas/"
# 
# vector.analysis <- c("diffrmse","diff.intercept","diff.lag")
# 
# vector.texts <- c("RMSE ratio ($RMSE^{QR} / RMSE^{AR} $) ", 
#                   "Coefficient difference between the non-autoregressive part (($\\hat{\\phi}_0 + z_\\alpha  \\hat{\\sigma}^2_\\varepsilon) - \\hat{\\beta}_0)$ ",
#                   "Diffeference between the autoregressive coefficients ($\\hat{\\phi} - \\hat{\\beta}$) ")
# file.name <- c("rmse", "intercept", "auto")
# for (alpha in vector.alpha) {
#   b.alpha <- alpha # dplyr doesnt understant if put (alpha == alpha)
#   # print(b.alpha)
#   rest.of.caption <- str_c("for estimating quantile
#                            $\\alpha = $ ", alpha, ". $\\phi$ stands for the autoregressive coefficient 
#                            and RSN is the signal to noise ratio. Details for these experiments can 
#                            be found on section \\ref{sec:simulation-ar1}")
#   # for latex table
#   alpha_out <- str_replace(alpha, "\\.", "")
#   
#   for (i in 1:length(vector.analysis)) {
#     f.name <- file.name[i]
#     var.analyse <- vector.analysis[i]
#     text.display <- vector.texts[i]
#     tmp <- results %>% filter(alpha == b.alpha) %>% dcast(phi ~ RSN, value.var = var.analyse)
#     names(tmp)[1] <- "$\\phi \\backslash RSN$"
#     row.names(tmp) <- NULL
#     tabela.latex <- xtable(tmp, caption = str_c(
#       text.display, rest.of.caption ), include.rownames=FALSE, 
#       digits = c(2,2,8,8,8,8,8), label = str_c("tab:sim-",f.name,"-",alpha_out))
#     print(tabela.latex, type = "latex", include.rownames = FALSE, file = str_c(pasta.destino,
#                                                                                "simulation-ar-", f.name, "-", alpha_out, ".tex"), sanitize.text.function=function(x){x})
#   }
#   # row.names(tmp) <- tmp[, 1]
#   # tmp <- tmp[, -1]
#   # levelplot(as.matrix(tmp), xlab = "phi", ylab = "RSN", main = str_c("alpha = ", alpha))
#   # heatmap(as.matrix(tmp))
#   # rownames(tmp) <- str_c("$\\beta_{", 0:12, "}$")
#   # colnames(tmp) <- str_c("K=",1:12)
# }  
# 

# pasta.destino <- "Documento Regressao Quantilica/Figuras/Simulations-Multi/"
# for (phi in vector.phi) {
#   for (RSN in vector.RSN) {
#     b.phi <- phi # dplyr doesnt understant if put (alpha == alpha)
#     b.rsn <- RSN
#     phi_out <- str_replace(phi, "\\.", "")
#     rsn_out <- str_replace(RSN, "\\.", "")
#     selected.results <- results %>% filter(phi == b.phi) %>% filter(RSN == b.rsn)
#     # plot(selected.results$diffrmse, type = "b", axes = FALSE)
#     pdf(str_c(pasta.destino,"phi", phi_out,"-rsn", rsn_out,".pdf"))
#     
#     plot(selected.results$alpha, selected.results$diffrmse, main = str_c("Phi = ", phi," and RSN = ", RSN),
#          xlab = "Alpha", ylab = "RMSE^QR / RMSE^AR", type ="b")
#     abline(h = 1, lty = 2, col = 2)
#     dev.off()
#     # as.character(selected.results$alpha)
#   }
# }
# 

# library(lattice)
# View(results)

