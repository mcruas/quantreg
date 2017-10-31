install.packages("devtools")
library("devtools")
install_github("hadley/l1tf")

library(l1tf)
plot(sp500$log, type='l', col='blue')
lines(l1tf(sp500$log, lambda = 50), col = "red", type='l')
lines(l1tf(sp500$log, prop = 0.01), col = "green3", type='l')
