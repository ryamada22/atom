library(igraph)

my.eigen.G.decomposition <- function(g){
	A <- get.adjacency(g)
	L <- diag(degree(g)) - A
	eigen.out <- eigen(L)
	V <- eigen.out[[2]]
	
	V.s <- sign(V)
	V.s[which(V.s==0)] <- 1
	n <- length(V[1,])
	V.s2 <- V.s[,n:1]
	cl <- matrix(0,n,n)
	cl[,1] <- V.s2[,1]
	for(i in 2:n){
		cl[,i] <- V.s2[,1:i] %*% 2^(0:(i-1))
	}
	
	X <- matrix(0,n,n)
	for(i in 1:n){
		tmp <- table(cl[,i])
		X[i,1:length(tmp)] <- sort(tmp,decreasing=TRUE)
	}
	
	A.list <- list()
	A.list[[1]] <- A
	for(i in 2:n){
		tmp <- matrix(1,n,n)
		tmp2 <- outer(cl[,i], cl[,i],"-")
		tmp[which(tmp2!=0)] <- 0
		A.list[[i]] <- A * tmp
	}
	n.delE <- rep(0,n)
	for(i in 2:n){
		n.delE[i] <- sum(A.list[[i-1]] - A.list[[i]])/2
	}
	curve <- cbind(X/n,1-cumsum(n.delE)/sum(n.delE))
	return(list(curve = curve,X=X, n.E=n.delE,cl=cl,A.list=A.list))
}

edge.list <- matrix(c(1,2,2,3,1,3,2,4,4,5,5,6,4,6,6,7,7,8,6,8,8,9,7,9),byrow=TRUE,ncol=2)
g <- graph.edgelist(edge.list,directed=FALSE)
plot(g)

out <- my.eigen.G.decomposition(g)
out

#g <- as.undirected(make_kautz_graph(4,4))
#plot(g)
#out <- my.eigen.G.decomposition(g)
#image(out$curve)



A1 <- diag(11)
A1 <- A1[c(2:11,1),]
A1 <- A1 + t(A1)
A1
out <- my.eigen.G.decomposition(graph.adjacency(A1))

A1 <- diag(11)

s1 <- sample(11)
s2 <- sample(11)
B1 <- A1[s1,]
B2 <- A1[s2,]
B <- B1 + t(B1) + B2 + t(B2)
#B <- sign(B)
g <- graph.adjacency(B)
plot(g)
out <- my.eigen.G.decomposition(g)
image(out$curve)
