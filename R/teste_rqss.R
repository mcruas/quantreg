library(quantreg)
n <- 200
x <- sort(rchisq(n,4))
z <- x + rnorm(n)
y <- log(x)+ .1*(log(x))^2 + log(x)*rnorm(n)/4 + z

vector_tau <- seq(0.05,0.95,0.05);

for (tau in vector_tau) {
  plot(x, y-z, main = tau)
  
  f.N  <- rqss(y ~ qss(x, constraint= "N") + z, tau = tau)
  f.I  <- rqss(y ~ qss(x, constraint= "I") + z, tau = tau)
  f.CI <- rqss(y ~ qss(x, constraint= "CI") + z, tau = tau)

  lines(x[-1], f.N $coef[1] + f.N $coef[-(1:2)])
  lines(x[-1], f.I $coef[1] + f.I $coef[-(1:2)], col="blue")
  lines(x[-1], f.CI$coef[1] + f.CI$coef[-(1:2)], col="red")

}










  for (tau in vector_tau) {
  i = which(tau ==vector_tau)
  qs_alpha_N[i] <- rqss(y ~ qss(x, constraint= "N") + z, tau = tau)
  qs_alpha_I[i] <- rqss(y ~ qss(x, constraint= "I") + z, tau = tau)
  qs_alpha_CI[i] <- rqss(y ~ qss(x, constraint= "CI") + z, tau = tau)
}

cbind(qs_alpha_N, qs_alpha_I, qs_alpha_CI)

lines(x[-1], f.N $coef[1] + f.N $coef[-(1:2)])
lines(x[-1], f.I $coef[1] + f.I $coef[-(1:2)], col="blue")
lines(x[-1], f.CI$coef[1] + f.CI$coef[-(1:2)], col="red")



## A bivariate example
data(CobarOre)
fCO <- rqss(z ~ qss(cbind(x,y), lambda= .08), data=CobarOre)
plot(fCO)






