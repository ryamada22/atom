V <- matrix(0,13,18)
F <- matrix(0,15,18)
E <- matrix(0,18,18)

V[1,c(12,14,7,3,4,5,6,10,8)] <- 1
V[2,c(14,7,16,11,5,6,8,10,12)] <- 1
V[3,c(14,7,16,9,8,11,10,15,12)] <- 1
V[4,c(14,2,16,11,9,13,12,10,15)] <- 1
V[5,c(14,2,16,11,15,13,6,18,12)] <- 1
V[6,c(14,2,16,4,15,6,18,17,13)] <- 1
V[7,c(1,2,16,15,4,6,17,8,18)] <- 1
V[8,c(1,2,3,4,5,6,17,8,18)] <- 1
V[9,c(1,2,3,4,6,5,10,12,8)] <- 1
V[10,c(1,2,16,4,5,6,8,12,18)] <- 1
V[11,c(14,2,16,4,5,6,8,10,12)] <- 1
V[12,c(14,2,16,4,15,6,12,18,8)] <- 1
V[13,c(14,2,16,11,10,15,6,12,8)] <- 1

theta <- (1:18) / 18 * 2 * pi
x <- cbind(cos(theta),sin(theta))

Vx <- V %*% x

plot(Vx,pch=20,col=1,asp=TRUE)

F[1,c(12,2,14,10,3,4,5,6,8)] <- 1
F[2,c(14,7,16,6,4,5,8,12,10)] <- 1
F[3,c(14,16,11,6,10,7,15,8,12)] <- 1
F[4,c(14,2,11,16,9,15,8,10,12)] <- 1
F[5,c(14,2,11,16,15,6,13,12,10)] <- 1
F[6,c(14,13,16,2,4,6,15,18,12)] <- 1
F[7,c(14,2,16,4,8,17,15,6,18)] <- 1
F[8,c(1,2,16,4,5,17,6,8,18)] <- 1
F[9,c(1,2,3,4,5,6,8,18,12)] <- 1
F[10,c(1,2,5,16,4,6,8,10,12)] <- 1
F[11,c(1,2,16,12,15,4,6,8,18)] <- 1
F[12,c(14,2,4,16,12,6,8,5,18)] <- 1
F[13,c(14,12,2,4,16,10,15,6,8)] <- 1
F[14,c(14,2,16,10,11,5,6,8,12)] <- 1
F[15,c(14,6,2,11,18,15,6,8,12,16)] <- 1

Fx <- F %*% x

points(Fx,pch=20,col=2)

E[1,c(2,10,4,3,16,5,6,8,12)] <- 1
E[2,c(14,2,7,6,8,4,5,10,12)] <- 1
E[3,c(14,16,4,11,5,6,8,10,12)] <- 1
E[4,c(10,12,14,6,2,7,11,16,8)] <- 1
E[5,c(14,6,8,16,9,11,15,12,10)] <- 1
E[6,c(14,2,15,16,13,11,8,10,12)] <- 1
E[7,c(14,16,11,2,10,6,15,18,12)] <- 1
E[8,c(14,12,2,8,16,15,13,6,18)] <- 1
E[9,c(14,2,17,16,4,12,15,18,6)] <- 1
E[10,c(1,2,16,6,14,4,15,8,18)] <- 1
E[11,c(17,1,2,16,4,6,8,12,18)] <- 1
E[12,c(1,2,3,4,5,6,8,16,18)] <- 1
E[13,c(1,2,4,12,18,5,6,8,10)] <- 1
E[14,c(1,14,16,2,4,5,6,8,12)] <- 1
E[15,c(2,16,5,15,4,6,18,8,12)] <- 1
E[16,c(14,2,12,4,16,10,6,8,18)] <- 1
E[17,c(14,10,16,2,5,6,15,8,12)] <- 1
E[18,c(14,16,2,4,11,15,12,6,8)] <- 1

Ex <- E %*% x

points(Ex, pch=20,col=3)

VFE <- rbind(V,F,E)
col.v <- c(rep(4,length(V[,1])),rep(5,length(F[,1])),rep(6,length(E[,1])))

X <- VFE %*% x


distmat <- as.matrix(dist(VFE, method="manhattan"))




el <- which(distmat==2,arr.ind=TRUE)


# ３角形は、V,F,Eを１つずつ頂点とするので、適切なトリオを選ぶ


tris <- matrix(0,0,3)

for(i in 1:length(V[,1])){
	for(j in 1:length(F[,1])){
		for(k in 1:length(E[,1])){
			tmp0 <- c(i,j+13,k+13+15)
			tmp <- distmat[tmp0,tmp0]
			diag(tmp) <- 2
			if(sum((tmp-2)^2) == 0){
				tris <- rbind(tris,tmp0)
			}
		}
	}
}

#tris <- rbind(c(1,1,2),c(1,2,2),c(2,2,3),c(2,14,3),c(2,3,4),c(3,3,5),c(3,4,5),c(4,4,6),c(4,5,6),c(5,5,7),c(5,15,7),c(5,6,8),c(6,6,9),c(6,7,9),c(7,7,10),c(7,11,10),c(7,11,11),c(7,8,11),c(8,8,12),c(8,9,12),c(9,9,13),c(9,10,13),c(9,1,1))

#tris[,2] <- tris[,2] + 13
#tris[,3] <- tris[,3] + 13 + 15

tri.col <- rep(0,length(tris[,1]))
for(i in 1:length(tris[,1])){
	tmp <- VFE[tris[i,],]
	if(prod(apply(tmp,2,sum) -1) == 0){
		tri.col[i] <- 2
	}else{
		tri.col[i] <- 3
	}
	
}


########
# Plabic graph G*
par(mfcol=c(2,2))

plot(X, pch=20,cex=2,col=col.v,asp=TRUE)
segments(X[el[,1],1], X[el[,1],2],X[el[,2],1],X[el[,2],2])

for(i in 1:length(tris[,1])){
	polygon(X[tris[i,c(1,2,3,1)],1],X[tris[i,c(1,2,3,1)],2],col=tri.col[i])
}

points(X,pch=20,col=col.v,cex=2)


#########
# Vの３角形はうまく描ける
plot(X, pch=20,cex=2,col=col.v,asp=TRUE)
segments(X[el[,1],1], X[el[,1],2],X[el[,2],1],X[el[,2],2])

for(i in 1:length(tris[,1])){
	polygon(X[tris[i,c(1,2,3,1)],1],X[tris[i,c(1,2,3,1)],2],col=tri.col[i])
}


distmat.v <- as.matrix(dist(V,method="manhattan"))
el.v <- which(distmat.v==4,arr.ind=TRUE)

points(Vx,pch=20)
segments(Vx[el.v[,1],1],Vx[el.v[,1],2],Vx[el.v[,2],1],Vx[el.v[,2],2],lw=2)

##########
# Fの3-valentグラフはこの方法ではうまく描けない
# 組み合わせ的距離では、Vを介して多数のFペアの距離が同じになるから
plot(X, pch=20,cex=2,col=col.v,asp=TRUE)
segments(X[el[,1],1], X[el[,1],2],X[el[,2],1],X[el[,2],2])

for(i in 1:length(tris[,1])){
	polygon(X[tris[i,c(1,2,3,1)],1],X[tris[i,c(1,2,3,1)],2],col=tri.col[i])
}


distmat.f <- as.matrix(dist(F,method="manhattan"))
el.f <- which(distmat.f==4,arr.ind=TRUE)

points(Fx,pch=20)
segments(Fx[el.f[,1],1],Fx[el.f[,1],2],Fx[el.f[,2],1],Fx[el.f[,2],2])

##########
# Eをつないだ箙はこの方法ではうまく描けない
# 組み合わせ的距離では、Vを介して多数のEペアの距離が同じになるから

plot(X, pch=20,cex=2,col=col.v,asp=TRUE)
segments(X[el[,1],1], X[el[,1],2],X[el[,2],1],X[el[,2],2])

for(i in 1:length(tris[,1])){
	polygon(X[tris[i,c(1,2,3,1)],1],X[tris[i,c(1,2,3,1)],2],col=tri.col[i])
}


distmat.e <- as.matrix(dist(E,method="manhattan"))
el.e <- which(distmat.e==4,arr.ind=TRUE)

points(Ex,pch=20)
segments(Ex[el.e[,1],1],Ex[el.e[,1],2],Ex[el.e[,2],1],Ex[el.e[,2],2])

