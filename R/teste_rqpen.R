library(rqPen)
# results_r = rq.lasso.fit.mult($X,$y,0.5, 1)
set.seed(123)
X <- matrix(rnorm(8000),ncol=10)
y <- 1 + X[,1] - 3*X[,5] + rnorm(800)
colMeans(X)
lambda = 0.2
debug(rq.lasso.fit)
rq.lasso.fit(X,y,lambda=lambda, tau = 0.7)

rho <- function(alpha, x) {
  if (x > 0) {
    alpha * x
  } else {
    (1- alpha) * x
  }
}

L <- function(alpha, y, q) {
  sum(rho(alpha, y - q))
}

l1_norm <- function(beta) sum(abs(beta))

