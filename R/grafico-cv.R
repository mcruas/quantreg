par(mar=c(0,0,0,0))
plot(0,0,xlim=c(0,28),ylim=c(0,1),
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",type="n", main = "Cross-validation scheme")
i <- 1
total_folds <- 5
for(j in 1:total_folds)
{
  test <- (6+j):26
  train <- 1:(5+j)
  arrows(0,1-j/total_folds,27,1-j/total_folds,0.05)
  points(train,rep(1-j/total_folds,length(train)),pch=19,col="blue")
  if(length(test) >= i)
    points(test[i], 1-j/total_folds, pch=19, col="red")
  if(length(test) >= i)
    points(test[-i], rep(1-j/total_folds,length(test)-1), pch=19, col="gray")
  else
    points(test, rep(1-j/total_folds,length(test)), pch=19, col="gray")
}
text(28,.95,"time")