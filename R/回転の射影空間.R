# Nは乱点の数を決める引数
# f1 は黄金比の値をデフォルトとする値。この値を黄金比から変えると格子としての良い性格がなくなる
# k=1をデフォルトとし、フィボナッチ格子を描くためにはk=1。ただし、その背景にあるらせんを描くためにはこの値をいじって、らせん上の点の数を増やす必要がある
fib.lattice.S2 <-function(N,f1=(sqrt(5)+1)/2,k=1){
	f2 <- f1-1 # 黄金比はx^2-x-1=0の１つの根。もう１つの根を取り出す
	#P <- 2*N+1
	i <- seq(from=-N, to=N,by=k)
	theta <- asin(2*i/(2*N+1))
	phi <- 2*pi*i*f2
	x <- cos(theta)*cos(phi)
	y <- cos(theta)*sin(phi)
	z <- sin(theta)
	return(cbind(x,y,z))
}

library(rgl)


# フィボナッチ格子点の座標
N <- 100
fl <- fib.lattice.S2(N)


plot3d(fl)

plot3d(rbind(fl,-fl))

library(geometry)

my.elFromConvhull <- function(h){
	n <- length(h[1,])
	ret <- matrix(0,0,2)
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			ret <- rbind(ret,h[,c(i,j)])
		}
	}
	ret <- t(apply(ret,1,sort))
	return(unique(ret))
}

my.projective.SO3 <- function(N,n){
	fl <- fib.lattice.S2(N)
	h1 <- convhulln(fl)
	flfl <- rbind(fl,-fl)
	h2 <- convhulln(flfl)
	elh1 <- my.elFromConvhull(h1)
	elh2 <- my.elFromConvhull(h2)
	
	theta <- seq(from=0,to=2*pi,length(n+1))
	theta <- theta[-1]
	
	n.S2 <- length(fl[,1])
	
	# 玉ねぎ格子の球面三角メッシュエッジ
	el1 <- elh1
	for(i in 2:(n-1)){
		el1 <- rbind(el1,elh1 + n.S2*(i-1))
	}
	# 玉ねぎ格子の半径方向エッジ
	el2 <- matrix(0,0,2)
	for(i in 2:n){
		el2 <- rbind(el2,cbind((1:n.S2)+n.S2 * (i-2),(1:n.S2)+n.S2 * (i-1)))
	}
	# 玉ねぎ格子の一番外側の球面のエッジ
	tmp <- which(elh2 > n.S2)
	elh2[tmp] <- elh2[tmp] - n.S2
	el3 <- elh2 + n.S2 * (n-1)
	el3 <- t(apply(el3,1,sort))
	el3 <- unique(el3)
	# 原点と一番内側の球面点との間のエッジ
	ori <- n.S2 * n + 1
	el4 <- cbind(rep(ori,n.S2),1:n.S2)
	
	el <- rbind(el1,el2,el3,el4)
	
	X <- fl * theta[1]
	for(i in 2:n){
		X <- rbind(X,fl*theta[i])
	}
	X <- rbind(X,rep(0,3))
	
	return(list(X=X,el=el,fl=fl,elfl=elh1))
}

out <- my.projective.SO3(200,30)

library(rgl)
plot3d(out$X)

segments3d(out$X[t(out$el),])

plot3d(out$fl)
segments3d(out$fl[t(out$elfl),])


library(igraph)

g <- graph.edgelist(out$el,directed=FALSE)
plot(g,vertex.label="",vertex.size=0.1)

adj <- get.adjacency(g)

degree(g)
table(degree(g))

#image(as.matrix(adj))

