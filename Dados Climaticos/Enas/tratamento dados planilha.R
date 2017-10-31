library(psych)
tmp <- read.table(file = "Dados Climaticos/Enas/enas sudeste matriz 1931-2014.csv", sep = ",", header=TRUE)
tmp2 <- as.vector(t(as.matrix(tmp)))
write.table(tmp2, file="Dados Climaticos/Enas/enas sudeste 1931-2014.csv", 
            sep = ",", row.names = FALSE, col.names = FALSE)
