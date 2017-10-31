library(forecast)
library(quantreg)
library(reshape2)
library(dplyr)
library(xtable)
library(stringr)
library(ggplot2)
# Simulate AR(1) process --------------------------------------------------
vector.alpha <- c(0.05, 0.1, 0.5, 0.9, 0.95)
vector.phi <- c(0.25, 0.5, 0.7, 0.9)
vector.RSN <- c(0.01, 0.05, 0.1, 0.5, 1)
n.rows <- length(vector.alpha) * length(vector.phi) * length(vector.RSN)
n = 10000; alpha = 0.05; i = 1; phi = 0.9; RSN = 0.1; phi0  = 1; serie = rep(0,n)

results <- data.frame(alpha = rep(0,n.rows), phi = rep(0,n.rows), RSN = rep(0,n.rows),
                      phi.hat = rep(0,n.rows), phi0.hat = rep(0,n.rows), sigma.e = rep(0,n.rows),
                      beta = rep(0,n.rows), beta0 = rep(0,n.rows), 
                      RMSE.arima = rep(0,n.rows), RMSE.qr = rep(0, n.rows),
                      diffrmse = rep(0,n.rows), diff.intercept = rep(0,n.rows), 
                      diff.lag = rep(0,n.rows))
                      
for (alpha in vector.alpha) {
  for (phi in vector.phi) {
    for (RSN in vector.RSN) {
      set.seed(123+i)
      
      # y_t = phi0 + phi * y_t-1 + epsilon
      
      
      sigma.e = (phi0/(1-phi)) * RSN
      serie.gen <- rep(0,2*n); serie[1] = (phi0/(1-phi))
      for (j in 2:(2*n)) {
        serie.gen[j] <- serie.gen[j-1] * phi + phi0 + rnorm(n=1,sd = sqrt(sigma.e))
      }
      testing <- serie.gen[(n+1):(2*n)]
      serie <- serie.gen[1:n]
      true.quantile <- testing[-1] + qnorm(alpha) * sqrt(sigma.e)
      # plot.ts(true.quantile)
      # write.table(serie, file = "R/serie_ar1.csv", row.names = FALSE, col.names = FALSE, sep = ",")
      
      # plot.ts(serie)
      
      # Estimate AR(1) model with normal ARIMA ----------------------------------
      
      # ar(serie, order.max = 1, method = "ols")
      fit1 <- Arima(serie, order = c(1,0,0), include.constant = TRUE)
      phi.hat <- coef(fit1)[1]
      phi0.hat <- coef(fit1)[2]*(1-phi)
      sigma2.hat <- fit1$sigma2
      # lm(serie ~ lag(serie) + 1)
      
      prediction.arima <- (testing[-n] * phi.hat + phi0.hat + qnorm(alpha) * sqrt(sigma2.hat))
      # lines(prediction.arima, col = 2)
      
      fit2 <- rq(serie[-1] ~ serie[-n] + 1, tau = alpha)
      beta0.hat <- coef(fit2)[1]
      beta.hat <- coef(fit2)[2]
      
      prediction.qr <- (testing[-n] * beta.hat + beta0.hat)
      
      error.arima <- sqrt(mean((prediction.arima - true.quantile)^2))
      error.qr <- sqrt(mean((prediction.qr - true.quantile)^2))
      
      diff.intercept <- (phi0.hat + qnorm(alpha) * sqrt(sigma2.hat)) - beta0.hat
      
      results[i, ] <- c(alpha, phi, RSN, phi.hat, phi0.hat, sigma2.hat, beta.hat, 
                        beta0.hat, error.arima, error.qr, error.qr/error.arima , 
                        diff.intercept , phi.hat - beta.hat)
      
      # plot(prediction.qr, prediction.arima, main = str_c("Phi = ", phi,", RSN = ", RSN, 
      #                                                    " and Alpha = ", alpha), 
      #                       # xlim = c(0, max(prediction.qr)), ylim = c(0,max(prediction.arima)),
      #      pch = 16 )
      #     abline(a = 0, b = 1, lty = 2, col = 2)
      # i = i + 1
      # plot(prediction.qr[1:200],main = str_c("Phi = ", phi,", RSN = ", RSN, 
      #                                       " and Alpha = ", alpha))
      # points(prediction.arima[1:200],pch = 3, col = 2)
    }
  }
}


# View results and prints to files ---------------------------------------------------------

# View(results)
write.csv(results, file = "results.csv", row.names = FALSE)
# names(results)
pasta.destino <- "Documento Regressao Quantilica/Tabelas/"

vector.analysis <- c("diffrmse","diff.intercept","diff.lag")

vector.texts <- c("RMSE ratio ($RMSE^{QR} / RMSE^{AR} $) ", 
                  "Coefficient difference between the non-autoregressive part (($\\hat{\\phi}_0 + z_\\alpha  \\hat{\\sigma}^2_\\varepsilon) - \\hat{\\beta}_0)$ ",
                  "Diffeference between the autoregressive coefficients ($\\hat{\\phi} - \\hat{\\beta}$) ")
file.name <- c("rmse", "intercept", "auto")
for (alpha in vector.alpha) {
  b.alpha <- alpha # dplyr doesnt understant if put (alpha == alpha)
  # print(b.alpha)
  rest.of.caption <- str_c("for estimating quantile
$\\alpha = $ ", alpha, ". $\\phi$ stands for the autoregressive coefficient 
and RSN is the signal to noise ratio. Details for these experiments can 
be found on section \\ref{sec:simulation-ar1}")
  # for latex table
  alpha_out <- str_replace(alpha, "\\.", "")
  
  for (i in 1:length(vector.analysis)) {
    f.name <- file.name[i]
    var.analyse <- vector.analysis[i]
    text.display <- vector.texts[i]
    tmp <- results %>% filter(alpha == b.alpha) %>% dcast(phi ~ RSN, value.var = var.analyse)
    names(tmp)[1] <- "$\\phi \\backslash RSN$"
    row.names(tmp) <- NULL
      tabela.latex <- xtable(tmp, caption = str_c(
            text.display, rest.of.caption ), include.rownames=FALSE, 
            digits = c(2,2,8,8,8,8,8), label = str_c("tab:sim-",f.name,"-",alpha_out))
    print(tabela.latex, type = "latex", include.rownames = FALSE, file = str_c(pasta.destino,
                                                   "simulation-ar-", f.name, "-", alpha_out, ".tex"), sanitize.text.function=function(x){x})
  }
    # row.names(tmp) <- tmp[, 1]
  # tmp <- tmp[, -1]
  # levelplot(as.matrix(tmp), xlab = "phi", ylab = "RSN", main = str_c("alpha = ", alpha))
  # heatmap(as.matrix(tmp))
  # rownames(tmp) <- str_c("$\\beta_{", 0:12, "}$")
  # colnames(tmp) <- str_c("K=",1:12)
}  
  

pasta.destino <- "Documento Regressao Quantilica/Figuras/Simulations/"
for (phi in vector.phi) {
  for (RSN in vector.RSN) {
    b.phi <- phi # dplyr doesnt understant if put (alpha == alpha)
    b.rsn <- RSN
    phi_out <- str_replace(phi, "\\.", "")
    rsn_out <- str_replace(RSN, "\\.", "")
    selected.results <- results %>% filter(phi == b.phi) %>% filter(RSN == b.rsn)
    # plot(selected.results$diffrmse, type = "b", axes = FALSE)
    pdf(str_c(pasta.destino,"phi", phi_out,"-rsn", rsn_out,".pdf"))
    
    plot(selected.results$alpha, selected.results$diffrmse, main = str_c("Phi = ", phi," and RSN = ", RSN),
         xlab = "Alpha", ylab = "RMSE^QR / RMSE^AR", type ="b")
    abline(h = 1, lty = 2, col = 2)
    dev.off()
    # as.character(selected.results$alpha)
  }
}

  
# library(lattice)
View(results)
