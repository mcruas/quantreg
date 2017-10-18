#install.packages("data.table")
library(data.table)
library(dplyr)

users = data.frame(fread("users.dat", sep=":")[, c(1,3,5,7,9), with=FALSE])
colnames(users) = c("UserID","Gender","Age","Occupation","Zip-code")
head(users)

ratings = data.frame(fread("ratings.dat", sep=":")[, c(1,3,5,7), with=FALSE])
colnames(ratings) = c("UserID","MovieID","Rating","Timestamp")
head(ratings)

movies <- as.vector(readLines("movies.dat"))
xx = gsub("\\|","::",movies)
head(xx)
yy = strsplit(xx,"::")
head(yy)

# máximo de elementos em uma única linha
aaa = NULL
for (i in 1:length(xx)){
  aaa[i] = length(yy[[i]])
}
max(aaa)

matrix.filmes = matrix(NA,length(movies),8)
for (i in 1:nrow(matrix.filmes)){
  matrix.filmes[i,seq(1,length(yy[[i]]))] = unlist(yy[i])
}
head(matrix.filmes)


# Criando a matriz para análise com dummy filmes
matrix.dados = data.frame(cbind(ratings[,-4],matrix(NA,nrow(ratings),18)))
colnames(matrix.dados)=c("UserID","MovieID","Rating","Action",
                         "Adventure","Animation","Children's",
                         "Comedy","Crime","Documentary","Drama",
                         "Fantasy","Film-Noir","Horror","Musical",
                         "Mystery","Romance","Sci-Fi","Thriller",
                         "War","Western")
                         


head(matrix.dados)
head(ratings)
head(matrix.filmes)


####################################################################################################################################
# Fazendo o crossing para as dummies
####################################################################################################################################
for(lin in 1:nrow(matrix.dados)){
  ind = length(which(is.na(matrix.filmes[which(matrix.filmes[,1]==matrix.dados$MovieID[lin]),])))
  aux = NULL
  for(i in 3:(max(aaa)-ind)){
    aux[i]=which(colnames(matrix.dados)==matrix.filmes[which(matrix.filmes[,1]==matrix.dados$MovieID[lin]),][i])
  }
  matrix.dados[lin,aux[c(-1,-2)]]=1
}

# Ao final
matrix.dados[is.na(matrix.dados)]=0

####################################################################################################################################
# Fazendo o crossing as idades
####################################################################################################################################
matrix.dados$UserID <- factor(matrix.dados$UserID)

zipcode <- users$`Zip-code`;userid <- users$UserID; Gender <- users$Gender; Age <- factor(users$Age); occupation <- factor(users$Occupation)
dummies = model.matrix(users$UserID ~ Gender + Age + occupation - 1)
tmp <- cbind(userid, dummies, zipcode)

library(dplyr)
db <- inner_join(x = matrix.dados, y = tmp, by = c("UserID" = "userid"), copy=TRUE )
save(db, file = "database.dad")
db <- merge(x = matrix.dados, y = tmp, by.x = "UserID", by.y = "userid")


head(ratings)
head(users)

