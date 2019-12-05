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
L <- runif(1)*3
#L <- 1.1
theta1 <- runif(1) * pi/2
#theta1 <- pi/2 * 0.95
A <- cos(theta1)
B <- sin(theta1)

x1 <- 1/2 * (L-A + B * sqrt((3+2*L*A-L^2)/(L^2-2*L*A+1)))
x2 <- 1/2 * (L-A - B * sqrt((3+2*L*A-L^2)/(L^2-2*L*A+1)))

y1 <- sqrt(1-x1^2)
y2 <- sqrt(1-x2^2)

C1 <- L-x1
D1 <- y1
C2 <- L-x2
D2 <- y2


(A-C1)^2 + (B-D1)^2
(A-C2)^2 + (B-D2)^2


#############

L <- 2 + runif(1)
thetas <- seq(from = 0, to = 1, length=101) * pi/2

V1 <- V2 <- V3 <- matrix(0,length(thetas),2)

thetas2_1 <- thetas2_2 <- rep(0,length(thetas))

ones2 <- ones3 <- rep(FALSE,length(thetas))

for(i in 1:length(thetas)){
	theta1 <- thetas[i]
	A <- cos(theta1)
	B <- sin(theta1)

	x1 <- 1/2 * (L-A + B * sqrt((3+2*L*A-L^2)/(L^2-2*L*A+1)))
	x2 <- 1/2 * (L-A - B * sqrt((3+2*L*A-L^2)/(L^2-2*L*A+1)))

	y1 <- sqrt(1-x1^2)
	y2 <- sqrt(1-x2^2)

	C1 <- L-x1
	D1 <- y1
	C2 <- L-x2
	D2 <- y2
	
	V1[i,] <- c(A,B)
	V2[i,] <- c(C1,D1)
	V3[i,] <- c(C2,D2)
	
	thetas2_1[i] <- acos(x1)
	thetas2_2[i] <- acos(x2)
	
	# ’†‰›‚Ìü•ª‚Ì’·‚³‚ª‚P‚Å‚ ‚é‚±‚Æ‚ÌŠm”F
	if(abs((A-C1)^2+(B-D1)^2-1)<10^(-6)){
		ones2[i] <- TRUE
	}
	if(abs((A-C2)^2+(B-D2)^2-1)<10^(-6)){
		ones3[i] <- TRUE
	}
}


plot(rbind(c(0,0),c(L,0),V1,V2,V3),xlim = c(-1,L),ylim=c(-1,L),pch=20,cex=0.1)
for(i in 1:length(thetas)){
	if(ones2[i]){
		segments(V1[i,1],V1[i,2],V2[i,1],V2[i,2],col=2)
		segments(V2[i,1],V2[i,2],L,0,col=1)
		segments(V1[i,1],V1[i,2],0,0,col=1)
	}
	if(ones3[i]){
		segments(V1[i,1],V1[i,2],V3[i,1],V3[i,2],col=3)
		segments(V3[i,1],V3[i,2],L,0,col=1)
		segments(V1[i,1],V1[i,2],0,0,col=1)
	}
}

plot(rbind(c(0,0),c(L,0),V1,V2,V3),xlim = c(-1,L),ylim=c(-1,L),pch=20,cex=0.1)
for(i in 1:length(thetas)){
	if(ones2[i]){
		segments(V1[i,1],V1[i,2],V2[i,1],V2[i,2],col=i)
		segments(V2[i,1],V2[i,2],L,0,col=i)
		segments(V1[i,1],V1[i,2],0,0,col=i)
	}
	if(ones3[i]){
		segments(V1[i,1],V1[i,2],V3[i,1],V3[i,2],col=i)
		segments(V3[i,1],V3[i,2],L,0,col=i)
		segments(V1[i,1],V1[i,2],0,0,col=i)
	}
}

plot(thetas[ones2],thetas2_1[ones2],type="l",col=2,ylim=c(0,pi/2))
points(thetas[ones3],thetas2_2[ones3],type="l",col=3)


#matplot(thetas,cbind(thetas2_1,thetas2_2),type="l")

plot(thetas + thetas2_1) # ˆê’è‚È‚í‚¯‚Å‚Í‚È‚¢