my.rSn <- function(n,d=3){
	X <- matrix(rnorm(n*d),ncol=d)
	X/sqrt(apply(X^2,1,sum))
}
n1 <- 1000
n2 <- 980
s1 <- my.rSn(n1)
s2 <- my.rSn(n2)
nr <- 100
rq <- my.rSn(nr,4)

my.rotmat.q <- function(q){
	s <- Mod(q)^{-2}
	x11 <- 1-2*s*(j(q)^2+k(q)^2)
	x12 <- 2*s*(i(q)*j(q) - k(q)*Re(q))
	x13 <- 2*s*(i(q)*k(q)+j(q)*Re(q))
	x21 <- 2*s*(i(q)*j(q)+k(q)*Re(q))
	x22 <- 1- 2*s (i(q)^2+k(q)^2)
	x23 <- 2*s*(j(q)*k(q)-i(q)*Re(q))
	x31 <- 2*s*(i(q)*k(q)-j(q)*Re(q))
	x32 <- 2*s*(j(q)*k(q)+i(q)*r(q))
	x33 <- 1-2*s*(i(q)^2+j(q)^2)
	
	ret <- matrix(c(x11,x21,x31,x12,x22,x32,x13,x23,x33),3,3)
	return(ret)
}
my.rotmat.v4 <- function(v){
	q <- v[1] + v[2] * Hi + v[3] * Hj +v[4] * Hk
	s <- Mod(q)^{-2}
	x11 <- 1-2*s*(j(q)^2+k(q)^2)
	x12 <- 2*s*(i(q)*j(q) - k(q)*Re(q))
	x13 <- 2*s*(i(q)*k(q)+j(q)*Re(q))
	x21 <- 2*s*(i(q)*j(q)+k(q)*Re(q))
	x22 <- 1- 2*s*(i(q)^2+k(q)^2)
	x23 <- 2*s*(j(q)*k(q)-i(q)*Re(q))
	x31 <- 2*s*(i(q)*k(q)-j(q)*Re(q))
	x32 <- 2*s*(j(q)*k(q)+i(q)*Re(q))
	x33 <- 1-2*s*(i(q)^2+j(q)^2)
	
	ret <- matrix(c(x11,x21,x31,x12,x22,x32,x13,x23,x33),3,3)
	return(ret)
}



library(onion)
rot.list <- list()
for(i in 1:length(rq[,1])){
	rot.list[[i]] <- my.rotmat.v4(rq[i,])
}

my.gaus.k.ip <- function(x,y,a=1){
	sumx <- apply(x^2,1,sum)
	sumy <- apply(y^2,1,sum)
	sumxy <- x %*% t(y)
	rep.sumx <- rep(sumx,length(y[,1]))
	rep.sumy <- rep(sumy,each=length(x[,1]))
	tmp <- rep.sumx + rep.sumy - 2 * c(sumxy)
	sum(1/(length(x[,1])*length(y[,1])) * exp(-a*tmp))
}


my.gaus.k.ip(s1,s2)

out <- rep(0,length(rq[,1]))
for(i in 1:length(rq[,1])){
	rot.s2 <- s2 %*% rot.list[[i]]
	out[i] <- my.gaus.k.ip(s1,rot.s2)
}