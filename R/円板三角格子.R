# Triangle lattice
my.tri.lattice <- function(a){
	v1 <- c(1,0)
	v2 <- c(cos(pi/3),sin(pi/3))
	v3 <- c(cos(2*pi/3),sin(2*pi/3))

	v <- rbind(v1,v2,v3)
	z <- (-a):a
	z3 <- as.matrix(expand.grid(z,z,z))
	
	unique(z3 %*% v)
}

# ‰~”Â
my.tri.lattice.circle <- function(b){
	a <- b*sqrt(2) + 3
	X <- my.tri.lattice(a)

	xlen <- sqrt(apply(X^2,1,sum))

	X. <- X[which(xlen < b),]
	return(X.)
}


len <- 2

xlen <- sqrt(apply(X^2,1,sum))

X.<- my.tri.lattice.circle(len)

plot(X.)

theta <- pi/6 / len

Rot <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)

X.2 <- X. %*% Rot

z <- sqrt(2* cos(pi/6) -1)

X0 <- cbind(X.,rep(0,length(X.[,1])))
X1 <- cbind(X.2,rep(z,length(X.2[,1])))

X <- rbind(X0,X1)
library(rgl)


d <- as.matrix(dist(X))
d1 <- which(abs(d-1) < 10^(-2),arr.ind=TRUE)

image(abs(d-1) < 10^(-3))

plot3d(X)

for(i in 1:length(d1[,1])){
	segments3d(rbind(X[d1[i,1],],X[d1[i,2],]))
}

