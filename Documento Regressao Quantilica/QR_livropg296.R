library(quantreg)

example(rq)

data(engel)

# renda e gastos com alimentação

fit1 <- rq(foodexp~income, tau=0.5, data=engel)
fit1
# para obter os ICs dos coeficientes! Esses ICs são calculados pelo método de inversão de rank descrito na seção 3.5.5.
summary(fit1)

# para extrair os resíduos (r1) ou os coeficientes (c1) da relação ajustada podemos escrever:
r1 = resid(fit1)
c1 = coef(fit1)


# 1) INFERÊNCIA
# Existem muitos métodos alternativos de conduzir inferências sobre os coeficientes de regressão quantílica.
# Como uma alternativa aos intervalos de confiança de inversão de rank, podemos obter uma tabela com os coeficientes,
# erros padrão, estatísticas t e p-valores

summary(fit1, se="nid")

# Os erros calculados nessa tabela são computados como descritos na Seção 3.2.3 para a formula de regressão quantílica Sanduíche,
# e usando a regra de largura de banda de Hall-Sheather. 

# Para obter a versão de Powell kernel da matrix de variância covariância estimada, deve-se especificar se="ker"
# Também é possível controlar a largura da banda. E também é possível computar os erros padrão via bootstrap se="boot". Existem 3 metodologias para isso conforme fim da pág299.

xx = income-mean(income)
fit1 = summary(rq(foodexp~xx,tau=2:98/100))


attach(engel)
plot(income,foodexp,cex=0.25,type="n",xlab="Household Income",ylab="Food Expenditure")
points(income,foodexp,cex=0.5,col="blue")
abline(rq(foodexp~income,tau=0.5),col="blue")
abline(lm(foodexp~income),lty=2,col="red")
taus=c(0.05,0.1,0.25,0.75,0.9,0.95)
for (i in 1:length(taus)){
  abline(rq(foodexp~income,tau=taus[i]),col="gray")
}
fit2 = summary(rq(foodexp~income,tau=c(0.05,0.1,0.25,0.75,0.9,0.95)))

postscript("engelcoef.ps", horizontal=FALSE, width=6.5, height=3.5)
plot(fit1)
dev.off()

latex(fit2,caption= "Engel's Law", transpose=TRUE)

z = rq(foodexp~income,tau=-1)
z$sol
z$dsol


# Estimar o quantil condicional de y em um valor específico de x é bem fácil. Vamos plotar o quantil empírico condicional estimado de gastos com comida para famílias que estejam nos percentis 0.1 e 0.9
# Posteriormente iremos plotar as densidades estimatimadas para os dois grupos usando o método de kernel adaptativo proposto por Silverman (1986). Essa função é particulamente conveniente porque permite
# massa desigual ser associada com observações tais como produzidas pelo processo de regressao quantilica.

x.pobre = quantile(income,0.1)
x.rico = quantile(income,0.9)
ps = z$sol[1,]
qs.pobre = c(c(1,x.pobre) %*% z$sol[4:5,])
qs.rico = c(c(1,x.rico) %*% z$sol[4:5,])
par(mfrow=c(1,2))
plot(c(ps,ps),c(qs.pobre,qs.rico),type="n",xlab=expression(tau),ylab="quantile")
plot(stepfun(ps,c(qs.pobre[1],qs.pobre)),do.points=FALSE,add=TRUE)
plot(stepfun(ps,c(qs.pobre[1],qs.rico)),do.points=FALSE,add=TRUE,col.hor="gray",col.vert="gray")
ap=akj(qs.pobre[1:270],qs.pobre[1:270],diff(ps))
ar=akj(qs.rico[1:270],qs.rico[1:270],diff(ps))
plot(c(qs.pobre[1:270],qs.rico[1:270]),c(ap$dens,ar$dens),type="n",xlab="Food Expenditure",ylab="Density")
lines(qs.rico[1:270],ar$dens,col="gray")
lines(qs.pobre[1:270],ap$dens,col="black")
legend(750,0.006,c("pobre","rico"),lty=c(1,1),col=c("black","gray"))



# Considerando que temos heterocedasticidade, vamos modelar com o log!
plot(income,foodexp,log="xy",xlab="Household Income", ylab="Food Expenditure")
taus=c(0.05,0.1,0.25,0.75,0.9,0.95)
abline(rq(log10(foodexp) ~ log10(income), tau=0.5),col="blue")
abline(lm(log10(foodexp) ~ log10(income)), lty=3, col="red")

for (i in 1:length(taus)){
  abline(rq(log10(foodexp) ~ log10(income),tau=taus[i]),col="gray")
}

# Vamos considerar algumas outras formas de teste. A primeira questão que devemos levantar é a seguinte: 
# A regressão quantilica estimada adapta-se a hipotese de mudança de locação que assume que todos os quantis condicionais possuem os mesmos parâmetros de inclinação?
# Para iniciar, suponha que nos estimamos os quartis ajustados para os dados de Engel conforme:
y = foodexp
x = income 

fit1 = rq(y~x, tau=0.25)
fit2 = rq(y~x, tau=0.5)
fit3 = rq(y~x, tau=0.75)

# suponha que queremos testar se o parâmetro de inclinação é o mesmo para as 3 regressões, isso é feito pela ANOVA

anova(fit1,fit2,fit3)





# previsão
data(airquality)
airq <- airquality[143:150,]
f <- rq(Ozone ~ ., data=airquality)
predict(f,newdata=airq)
f <- rq(Ozone ~ ., data=airquality,tau=1:3/4)
predict(f,newdata=airq)

