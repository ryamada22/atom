# â·‚Ì—LŒü•Ó‚P–{‚ÅŒ‹‚Î‚ê‚é‚Q’¸“_ vi, vj‚ÌŒÝŠ·‚Í
# i,j,i,j,i

my.B.mut <- function(B,k){
  new.B <- B
  n <- length(B[,1])
  new.B[k,] <- - B[k,]
  new.B[,k] <- - B[,k]
  for(i in 1:n){
    for(j in 1:n){
      if(i != k){
        if(j != k){
          new.B[i,j] <- B[i,j] + sign(B[i,k]) * max(B[i,k]*B[k,j],0)
        }
      }
    }
  }
  return(new.B)
}
my.mat.mut.B <- function(B,k){
	 b <- B[,k]
	b. <- cbind(b,rep(0,length(b)))
	b. <- apply(b.,1,max)
	b.[k] <- -1
	ret <- diag(rep(1,length(b)))
	ret[,k] <- b.
	return(ret)
}

library(igraph)

gs <- list()

d <- 10
A1 <- diag(rep(1,d))

A1 <- A1[,c(2:d,1)]
A1[1,d] <- 0

B1 <- A1 - t(A1)
A1
B1

gs[[1]] <- graph.adjacency(A1)
As <- list()
As[[1]] <- A1
Bs <- list()
Bs[[1]] <- B1
Ws <- list()

plot(gs[[1]])


k <- 4
#n.mut <- 5
#ks <- c(k,k+1,k,k+1,k)
k <- 4
k1 <- 5
k2 <- 6
k3 <- 7
#ks <- c(k,k1,k,k2,k1,k3,k2,k3,k2,k)
ks <- c(k,k1,k,k2,k3,k2,k3,k2,k1,k)
n.mut <- length(ks)



#n.mut <- 9*100
#ks <- sample(1:d,n.mut,replace=TRUE)
for(i in 1:n.mut){
	k <- ks[i]
	B1 <- my.B.mut(B1,k)
	Ws[[i]] <- my.mat.mut.B(B1,k)
	A1 <- B1
	A1[which(A1 < 0)] <- 0
	gs[[i+1]] <- graph.adjacency(A1)
	As[[i+1]] <- A1
	Bs[[i+1]] <- B1
}

par(mfrow=c(2,3))

for(i in 1:(n.mut+1)){
	plot(gs[[i]],edge.arrow.size=1)
}

ev.mat <- matrix(sapply(As,function(x){eigen(x)[[1]]}),byrow=TRUE,ncol=d)

matplot(t(ev.mat),type="l")


B <- matrix(sample(c(0,1,2,3),d^2,replace=TRUE),d,d)
#B <- matrix(sample(c(0,2),d^2,replace=TRUE),d,d)
#B <- matrix(sample(c(0,1),d^2,replace=TRUE),d,d)
B <- B - t(B)
B1 <- B
A1 <- B1
A1[which(A1 < 0)] <- 0

gs[[1]] <- graph.adjacency(A1)
As <- list()
As[[1]] <- A1
Bs <- list()
Bs[[1]] <- B1
Ws <- list()

plot(gs[[1]])


add <- which(Bs[[1]] == 1,arr.ind=TRUE)
k <- add[1,2]
k2 <- add[1,1]

#tmp <- sample(1:d,2)
#k <- tmp[1]
#k2 <- tmp[2]
#k <-4
#k2 <- 5
n.mut <- 5
ks <- c(k,k2,k,k2,k)




#n.mut <- 9*100
#ks <- sample(1:d,n.mut,replace=TRUE)
for(i in 1:n.mut){
	k <- ks[i]
	B1 <- my.B.mut(B1,k)
	Ws[[i]] <- my.mat.mut.B(B1,k)
	A1 <- B1
	A1[which(A1 < 0)] <- 0
	gs[[i+1]] <- graph.adjacency(A1)
	As[[i+1]] <- A1
	Bs[[i+1]] <- B1
}

par(mfrow=c(2,3))

for(i in 1:(n.mut+1)){
	plot(gs[[i]],edge.arrow.size=1)
}

P <- diag(rep(1,d))
s <- 1:d
s[k2] <- k
s[k] <- k2
P <- P[s,]

P %*% Bs[[1]] %*% t(P) - Bs[[6]]
