Xpre <- diag(rep(1/sqrt(2),4))

Xpre <- Xpre - apply(Xpre,1,mean)

a <- 1
b <- 2
t <- pi/3

Xpost <- rbind(c(a,b),c(a+1,b),c(a+1+cos(t), b+sin(t)),c(a+cos(t),b+sin(t)))

plot(Xpost)
segments(Xpost[,1],Xpost[,2],Xpost[c(2:4,1),1],Xpost[c(2:4,1),2])

segments(Xpost[1,1],Xpost[1,2],Xpost[2,1],Xpost[2,2])
segments(Xpost[2,1],Xpost[2,2],Xpost[3,1],Xpost[3,2])
segments(Xpost[3,1],Xpost[3,2],Xpost[4,1],Xpost[4,2])
segments(Xpost[4,1],Xpost[4,2],Xpost[1,1],Xpost[1,2])

# Xpost ‚ð4ŽŸŒ³‹óŠÔÀ•W‚É‚·‚é‚×‚­A‘æ‚RE‘æ‚SŽ²‚ÌÀ•W‚ð‚O‚Æ‚µ‚Ä’Ç‰Á‚·‚é
Xpost. <- matrix(0,4,4)
Xpost.,[1:2] <- Xpost

# “K“–‚È‚SŽŸŒ³‰ñ“]s—ñ‚ðì‚é
library(GPArotation)

R <- Random.Start(4)

# M %*% Xpre = R %*% Xpost
M <- R %*% Xpost. %*% solve(Xpre)

eigen(M) # 2ŒÂ‚ÌŒÅ—L’l‚Í0

# X1,X2 •½–Ê‚Ì‚QƒxƒNƒgƒ‹
v1 <- c(runif(2),0,0)
v2 <- c(runif(2),0,0)

Mv1 <- M %*% v1
Mv2 <- M %*% v2

sum(v1 * v2) - sum(Mv1 * Mv2)

###############

Xpre <- diag(rep(1/sqrt(2),4))

Xpre <- Xpre - Xpre[,1]

a <- 0
b <- 0
t <- pi/3

Xpost <- rbind(c(a,b),c(a+1,b),c(a+1+cos(t), b+sin(t)),c(a+cos(t),b+sin(t)))

plot(Xpost)
segments(Xpost[1,1],Xpost[1,2],Xpost[2,1],Xpost[2,2])
segments(Xpost[2,1],Xpost[2,2],Xpost[3,1],Xpost[3,2])
segments(Xpost[3,1],Xpost[3,2],Xpost[4,1],Xpost[4,2])
segments(Xpost[4,1],Xpost[4,2],Xpost[1,1],Xpost[1,2])

# Xpost ‚ð4ŽŸŒ³‹óŠÔÀ•W‚É‚·‚é‚×‚­A‘æ‚RE‘æ‚SŽ²‚ÌÀ•W‚ð‚O‚Æ‚µ‚Ä’Ç‰Á‚·‚é
Xpost. <- matrix(0,4,4)
Xpost.[,1:2] <- Xpost

# “K“–‚È‚SŽŸŒ³‰ñ“]s—ñ‚ðì‚é
library(GPArotation)

R <- Random.Start(4)

# M %*% Xpre = R %*% Xpost
M <- R %*% Xpost. %*% solve(Xpre)

eigen(M) # 2ŒÂ‚ÌŒÅ—L’l‚Í0

# X1,X2 •½–Ê‚Ì‚QƒxƒNƒgƒ‹
v1 <- c(runif(2),0,0)
v2 <- c(runif(2),0,0)

Mv1 <- M %*% v1
Mv2 <- M %*% v2

sum(v1 * v2) - sum(Mv1 * Mv2)

#######

theta <- runif(1) * pi

K <- matrix(0,4,4)
diag(K) <- cos(theta)
K[1,2] <- K[2,3] <- K[3,4] <- -sin(theta)
K[2,1] <- K[3,2] <- K[4,3] <- sin(theta)

K
Xpre <- diag(rep(1/sqrt(2),4))

K %*% Xpre

########

A <- matrix(0,4,4)
diag(A) <- 1
A[1,4] <- A[2,1] <- A[3,2] <- A[4,3] <- -1
A <- A / sqrt(2)

thetas <- runif(2) * pi
thetas <- c(thetas,thetas+pi)

B <- rbind(cos(thetas),sin(thetas))
B <- rbind(B,matrix(0,2,4))
plot(t(B[1:2,]))

# K %*% A = B
library(matlib)
K <- B %*% Ginv(A)

out <- K %*% Xpre

plot(t(out[1:2,]))


F <- 1/sqrt(2) * matrix(c(1,0,-1,-1,1,0,0,-1,1),byrow=TRUE,3,3)

det(F)
F %*% t(F)

###

n <- 10

costhetas <- cos(2*pi/n * (1:n))
thetas <- 2*pi/n * (1:n)

X <- cbind(cos(thetas),sin(thetas))


cumsumX <- apply(X,2,cumsum)

my.right.npoly <- function(n){
	thetas <- 2*pi/n * (1:n)
	edges <- cbind(cos(thetas),sin(thetas))
	X <- apply(edges,2,cumsum)
	return(list(X=X,edges=edges))
}
out <- my.right.npoly(n)

plot(out$X)

# ‘S•Ó‚Ì’·‚³‚Í‚P
dist(out$X)

my.zitter.npoly <- function(x,s=0.5){
	thetas <- 2*pi/n * (1:n)
	thetas <- thetas + rnorm(n) * s

	edges <- cbind(cos(thetas),sin(thetas))
	X <- apply(edges[1:(n-2),],2,cumsum)
	xn <- c(0,0)
	xn_2 <- X[n-2,]
	
	d <- sqrt(sum((xn-xn_2)^2))
	
	perp <- c(xn_2[2],-xn_2[1])
	
	k <- sqrt(1-(d/2)^2) / sqrt(sum(perp^2))
	
	xn_1 <- xn_2/2 + k * perp
	
	X <- rbind(X,xn_1,xn)
	edges <- apply(X[c(1:n,1),],2,diff)
	return(list(X=X,edges=edges))
}

outz <- my.zitter.npoly(n,s=0.1)

plot(outz$X)

# ‘S•Ó‚Ì’·‚³‚Í‚P
dist(outz$X)
segments(outz$X[,1],outz$X[,2],outz$X[c(2:n,1),1],outz$X[c(2:n,1),2])



A <- matrix(0,n,n)
diag(A) <- 1
for(i in 1:n){
	if(i==1){
		A[i,n] <- -1
	}else{
		A[i,i-1] <- -1
	}
}
A <- A / sqrt(2)

B <- t(outz$edges)
B <- rbind(B,matrix(0,n-2,n))
plot(t(B[1:2,]))

# K %*% A = B
library(matlib)
K <- B %*% Ginv(A)

Xpre <- 1/sqrt(2) * diag(rep(1,n))
out <- K %*% Xpre

plot(t(out))

range(as.matrix(dist(t(out))) - as.matrix(dist(outz$X)))

