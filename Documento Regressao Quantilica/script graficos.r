library(ggplot2)
library(ggfortify)
x1 = seq(0,1,0.01)
x2 = seq(1,2,0.01)
x = c(x1,x2)
y1 = 1 + rnorm(101,sd = 0.05)
y2 = x2*1 + rnorm(101,sd = 0.05)
y = c(y1,y2)
plot(x, y)
library(quantreg)
fit = rq(y1 ~ x1,tau=0.2)
lines(x1, coef(fit)[1] + coef(fit)[2] * x1, type = "l", col = 2)
fit = rq(y2 ~ x2,tau=0.2)
lines(x2, coef(fit)[1] + coef(fit)[2] * x2, type = "l", col = 2)

fit = rq(y ~ x,tau=0.2)
qplot(x, y) + geom_abline(intercept = coef(fit)[1], slope  =  coef(fit)[2], colour = "red")

write.table(x, file = "tmp/x_nonlin.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(y, file = "tmp/y_nonlin.csv", row.names = FALSE, col.names = FALSE, sep = ",")

n <- 100000
samp.t <- rt(n, df = 5)
samp.n <- rnorm(n)
# plot(density(samp.t), xlim = c(-4,4), ylim = c(0,0.5))
# lines(density(samp.n),col = 2, xlim = c(-4,4), ylim = c(0,0.5))
media.t <- mean(samp.t)
media.n <- mean(samp.n)
var.t <- sd(samp.t)^2
var.n <- sd(samp.n)^2
(k.t <- 1/n*sum((samp.t-media.t)^4) / (var.t)^2)
(k.n <- 1/n*sum((samp.n-media.n)^4) / (var.n)^2)
plot(ecdf(samp.t))
lines(ecdf(samp.n), col=2)

x = seq(-4,4,0.01)
plot(x,dt(x, df = 5), type = "l", ylim = c(0,0.4),col=4)
lines(x, dnorm(x), col = 2)
abline(h = 0)
