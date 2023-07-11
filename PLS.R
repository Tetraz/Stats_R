rm(list=ls())


#laod packages
library(pls)
library(plsRglm)
library(matlib)




# Define the data
x1 <- c(1, 3, 8,7,3,5)
x2 <- c(6, 2, 5,1,3,0)

x1 <- (x1-mean(x1))/sd(x1)
x2 <- (x2- mean(x2))/sd(x2)
y <- c(3, 4, 5,6,0,-2)

# Combine the data into a data frame
data <- data.frame(x1, x2, y)

# Fit the pls model
model <- plsr(y ~ x1 + x2, data=data, ncomp = 2)
model2 <- plsR(y ~ x1 + x2, data=data, ncomp = 2)


# Print the model summary
summary(model)
model$coefficients
model2$Coeffs

#base 1 c'est la base des x
model$loadings#c'est les coefficients de regression entre les x1 et les t1 dans l'espace des observations de dimension n
model2$pp
.

#wiegths
base2 <- model$loading.weights#coordonnées des x dans la base 2 cette base est orthogonal. inverse= transposee
model2$wwnorm#coordonnée base 2

#la base 3 c'est le w* de simca: c'est une matrice de changeemnt de base tq T=XW*. Cette matrice n'est pas orthogonale.
wstar <- model2$wwetoile#ok


#scores
score <- model$scores##coordonnées des obs dans la base 3
model2$tt#score ok


#loadings Y
model$Yloadings


#propriétes sur les composantes
t1 <- score[,1]
t2 <- score[,2]
cor(t1,t2)

dt1 <- sqrt(sum(t1*t1))
dt2 <- sqrt(sum(t2*t2))
ps <- sum(t1*t2)/(dt1*dt2)
#les composantes sont orthogonal dans l'espace a n dimension et leur covariane est nulle.





#base 2
m <- matrix(base2,ncol=2)#base 2 est orthogonal donc l'inverse=transposée
inv_base2 <- t(m)#coordonnées des composantes dans la base 1
inv(m)

#base 3
inv_wstar <- inv(wstar)#


#verification th=X*w*
data_m <- as.matrix(data[,-3])
data_m%*%wstar
(t1 <- model2$tt[,1])



#Deflation
norme <- as.vector((t1%*%t1))
projection <- matrix(t1)%*%t(matrix(t1))/norme
xdeflated <- data_m-projection%*%data_m


#on verifie que Th=Xdeflatedh*compo_h de la base 2
(t2_deflated <- xdeflated%*%base2[,2])#on regarde comment xdelated se projete sur la 2eme composante 
t2








##################################################################################""
#GRAPH
###################################################################################





#dans la base d'origine
plot(x1,x2)#données avant delfation
points(xdeflated[,1],xdeflated[,2],col="red")#données après delfation
arrows(0,0,1,0,col="grey80")#abse1
arrows(0,0,0,1,col="black")#base1
arrows(0,0,inv(model2$wwnorm)[1,1],inv(model2$wwnorm)[1,2],col="red")
arrows(0,0,inv(model2$wwnorm)[2,1],inv(model2$wwnorm)[2,2],col="red")
arrows(0,0,inv_wstar[1,1],inv_wstar[1,2],col="lightblue")
arrows(0,0,inv_wstar[2,1],inv_wstar[2,2],col="blue")
# prepare "circle data"
radius = 1
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
# draw the circle
lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y)



#EN UTILISANT LE THEOREME DE THALES? ON RETORUVE BIEN QUE LES COORDONN2ES DE x1 SUR C1 ET COORDON2ES DE X2 SUR c1 SONT LES MEMES QUELQUE SOIT LA BASE 5 ABSE 2 OU BASE 3°


#dans la base 3
plot(score[,1],score[,2])#données avant delfation
arrows(0,0,1,0,col="lightblue")#base3
arrows(0,0,0,1,col="blue")#base3
arrows(0,0,wstar[1,1],wstar[1,2],col="grey80")
arrows(0,0,wstar[2,1],wstar[2,2],col="black")



#biplot

#pour les variables, corrélations dans l'espace des observations à n dimension.
cor(x1,t1)
cor(x1,t2)


#pour les scores, ils sont ramené dans le cercle de rayon 1 en divisant tous les scores par la plus grande valeur absolue.






