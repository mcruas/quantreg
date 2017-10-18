library(tidyverse); library(readr); library(lattice)


lagmatrix <- function(x, lags) {
  n = length(x)
  I = (lags+1):(n)
  x_estim = matrix(0,n-lags, lags)
  for (i in 1:lags) {
    x_estim[ ,i] = x[I - i]
  }
  colnames(x_estim) <- paste0("lag",1:lags)
  return(x_estim)
}



WindData <- read_csv("Dados Climaticos/Dados-kaggle/DadosKg2014Adj_edited.csv", col_types = "iiiiidddddddddd")
n = nrow(WindData)
WindData_melt <- WindData %>% gather(`W1`,`W2`,`W3`,`W4`,`W5`,`W6`,`W7`,`W8`,`W9`,`W10`,key = "location", value = "wind")
month_mean <- WindData_melt %>% group_by(MONTH1, YEAR1,location) %>% summarise(mean = mean(wind))

n_lags <- 48

# Adds lags 
WindData_oneserie <- WindData_melt %>% filter(location == "W1") # get only values from data.frame
WindData_oneserie_withlags <- cbind(WindData_oneserie[-(1:n_lags),],lagmatrix(WindData_oneserie$wind, n_lags))


# Creates subset of data to use for estimation purposes
n_mmlags <- 6 # defines the number of averages to include in the end of the dataset 
Data_estimation <- WindData_oneserie_withlags %>% filter(MONTH1 >= 3, MONTH1 <= 10, YEAR1 == 2013) %>% mutate(MM0 = 0, MM1 = 0, MM2 = 0, MM11 = 0, MM12 = 0, MM13 = 0) # creates new columns to the end of the data to put statistics, which will be used for estimation, as well as define the 
Initial_column_MM <- which(colnames(Data_estimation) == "MM0") # defines the initial column to write Month Statistics

# For each observation, the lag is added to the last columns
for (i in 1:nrow(Data_estimation)) { 
  yeari = Data_estimation[i,"YEAR1"]; monthi = Data_estimation[i,"MONTH1"]; locationi = Data_estimation[i,"location"]
  
  # makes data.frame of lags to use
  
  info_lags <- data.frame(monthi = rep(0,n_mmlags), yeari = rep(0,n_mmlags), locationi = rep("a",n_mmlags), stringsAsFactors = FALSE)
  
  info_lags[1,1] <- monthi;       info_lags[1,2] <-  yeari;       info_lags[1,3] <-  locationi  # MM0
  info_lags[2,1] <- monthi - 1;   info_lags[2,2] <-  yeari;       info_lags[2,3] <-  locationi      # MM1
  info_lags[3,1] <- monthi - 2;   info_lags[3,2] <-  yeari;       info_lags[3,3] <-  locationi    # MM2
  info_lags[4,1] <- monthi + 1;   info_lags[4,2] <-  yeari - 1;   info_lags[4,3] <-  locationi    # MM11
  info_lags[5,1] <- monthi;      info_lags[5,2] <-  yeari - 1;   info_lags[5,3] <-  locationi    # MM12
  info_lags[6,1] <- monthi -1;   info_lags[6,2] <-  yeari - 1;   info_lags[6,3] <-  locationi    # MM13
  for (mmlags in 1:n_mmlags) { # For each MMlag, finds statistic on table and puts it into its column
    Data_estimation[i,Initial_column_MM - 1 + mmlags] <-  (month_mean %>% filter(MONTH1 == info_lags[mmlags,1], YEAR1 == info_lags[mmlags,2], location == info_lags[mmlags,3]) %>% ungroup() %>% select(mean))[[1]] # Finds statistic and puts into the dataset
  }
  
  # Data_estimation[i,]
  print(i)
}
# write_csv(x = Data_estimation, path = "Dados Climaticos/Dados-kaggle/Data_kaggle2014_oneserie_estimation.csv")
# getwd()


###### VARIABLE SELECTION ############################
Data_estimation <- read_csv(file = "Dados Climaticos/Dados-kaggle/Data_kaggle2014_oneserie_estimation.csv")
tau = 0.95

library(rqPen); library(stringi)
col_X <- stri_detect_fixed(colnames(Data_estimation), "MM") | stri_detect_fixed(colnames(Data_estimation), "lag")
X = as.matrix(Data_estimation[,col_X])[,-1] # The matrix of predictors
y = Data_estimation$wind # Vector of response values

# cross-validation estimation
fit <- cv.rq.pen(x = X,
          y = as.matrix(Data_estimation$wind),
          tau=tau) # debug(cv.rq.pen)

#gets all values of lambda used on cross-validation
values_cv <- fit$cv[,1]
keep_betas_cv <- matrix(0,ncol = ncol(X)+1, nrow = length(values_cv))
rownames(keep_betas_cv) <- values_cv; colnames(keep_betas_cv) <- c("intercept",colnames(Data_estimation)[col_X])

# Reruns the lasso estimates to every lambda used on the cross-validation procedure, keeping its variable values
for (i in 1:length(values_cv)) {
  lambda <- values_cv[i] 
  fit_i <- rq.lasso.fit(x = as.matrix(Data_estimation[,col_X]),
                        y = as.matrix(Data_estimation$wind),
                        tau=tau, lambda = lambda)
  keep_betas_cv[i, ] <- fit_i$coefficients
}
  
View(keep_betas_cv)

cross




fit_i$coefficientsfit_i <- rq.lasso.fit(x = as.matrix(Data_estimation[,col_X]),
              y = as.matrix(Data_estimation$wind),
              tau=tau, lambda = 0.0013530478)
fit_i$coefficients

plot(density(diff(y)))
media = mean(diff(y))
desv_pad = sd(diff(y))
lines(seq(-2,2,0.01),dnorm(seq(-2,2,0.01), mean = media, sd = desv_pad))

###############################33
# library(quantregGrowth)
# data(growthData)
# taus <- c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
# model_gcrq <- gcrq(y ~ ps(x, monotone = TRUE, lambda = 0) + z, tau = taus, data = growthData)
# plot(model_gcrq)
densityplot(~W1 | factor(MONTH1),data = WindData, main="Wind power density comparison across different months", xlab = "")
