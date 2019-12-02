N <- 50

psi <- 2 * pi * (1:N)/N
x <- cbind(cos(psi),sin(psi))
x. <- apply(x,2,cumsum)

plot(x.)
for(i in 1:(length(psi)-1)){
	segments(x.[i,1],x.[i,2],x.[i+1,1],x.[i+1,2])
}
segments(x.[length(psi),1],x.[length(psi),2],x.[1,1],x.[1,2])


X0 <- c(0,0)
X1 <- c(0.7,0)
theta <- seq(from=0,to=1,length=101) * 2 * pi
theta <- theta[-1]

Z <- cbind(X1[1] + cos(theta),X1[2] + sin(theta))

plot(rbind(X0,X1,Z))
points(Z,type="l")
points(rbind(X0,X1),pch=20,cex=1,col=2)
for(i in 1:length(theta)){
	segments(X0[1],X0[2],Z[i,1],Z[i,2])
	#segments(X1[1],X1[2],Z[i,1],Z[i,2],col=2)
}
L <- sqrt(apply(Z^2,1,sum))
plot(theta,L,type="l",ylim=c(0,2))
abline(h=0)



k <- 0.4

x <- sqrt(1 + 2 * k * cos(theta) + k^2)

plot(theta,x)

n <- 3
ks <- seq(from=0,to=1,length = 101) * n

X <- matrix(0,length(theta),length(ks))

for(i in 1:length(ks)){
	k <- ks[i]
	x <- sqrt(1 + 2 * k * cos(theta) + k^2)
	X[,i] <- x
}

matplot(theta,X,type="l")


# 3ü•ªŠÔ‹——£‚ª L ‚Ì‚Æ‚«AŠÔ‚Ì‚Q“_‚ð“®‚©‚µ‚½‚¢
L <- runif(3)

