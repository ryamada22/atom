my.Epair <- function(el){
	s.el <- apply(el,1,sort)
	mv <- max(el)+1
	V <- s.el[1,] * mv + s.el[2,]
	tmp <- outer(V,V,"-")
	diag(tmp) <- 1
	pairs <- which(tmp==0,arr.ind=TRUE)
	return(pairs)
}

my.Einfo <- function(f){
	el <- rbind(f[,1:2],f[,2:3],f[,c(3,1)])
	E.V <- matrix(0,length(el[,1]),max(f))
	for(i in 1:length(el[,1])){
		E.V[i,el[i,1]] <- -1
		E.V[i,el[i,2]] <- 1
	}
	
	ne <- length(el[,1])
	edge.trio <- matrix(1:length(el[,1]),ncol=3)
	
	e.pair <- my.Epair(el)
	# 順方向か逆方向か
	e.dir <- rep(1,ne)
	e.dir[which(el[,1] < el[,2])] <- -1
	# エッジのシークエンスID
	e.id <- rep(0,ne)
	e.id[which(el[,1] > el[,2])] <- 1:(ne/2)
	for(i in 1:ne){
		if(e.id[i]==0){
			e.id[i] <- e.id[e.pair[which(e.pair[,1]==i),2]]
		}
	}
	# eid.dir
	# エッジID(第一カラム)とその順逆情報(第二カラム)
	eid.dir <- cbind(e.id,e.dir)
	
	# eidの始点・終点情報
	eid.stend <- matrix(0,ne/2,2)
	for(i in 1:length(eid.stend[,1])){
		tmp <- which(eid.dir[,1] * eid.dir[,2] == i)
		eid.stend[i,] <- el[tmp,]
	}
	eid.E.V <- matrix(0,length(eid.stend[,1]),max(f))
	for(i in 1:length(eid.E.V[,1])){
		eid.E.V[i,eid.stend[i,1]] <- -1
		eid.E.V[i,eid.stend[i,2]] <- 1
	}
	regM <- diag(rep(1,ne/2))
	regM.v <- diag(rep(1,ne/2))
	for(i in 1:length(edge.trio[,1])){
		tmp <- edge.trio[i,]
		tmp <- c(tmp,tmp[1])
		for(j in 1:3){
			e1 <- tmp[j]
			e2 <- tmp[j+1]
			
			eid1 <- eid.dir[e1,1]
			eid2 <- eid.dir[e2,1]
			
			edir1 <- eid.dir[e1,2]
			edir2 <- eid.dir[e2,2]
			
			regM[eid1,eid2] <- regM[eid2,eid1] <- 1
			regM.v[eid1,eid2] <- regM.v[eid2,eid1] <- -0.5
			if(edir1 * edir2 == -1){
				regM.v[eid1,eid2] <- regM.v[eid2,eid1] <- 0.5
			}
		}
	}
	# eid は3k本のエッジ、e2は6k本のエッジ
	return(list(eid.stend = eid.stend,eid.EV = eid.E.V,regM = regM, regM.v = regM.v, e2.stend = el,e2.EV = E.V,e.pair = e.pair,e2.iddir = eid.dir,e2.trio = edge.trio))
}

e.info <- my.Einfo(f)

my.rsvd.tri <- function(X3,f,k=3,eps=10^(-10),maxiter=10000){
	e.info <- my.Einfo(f)
	
	regM <- e.info$regM
	regM.v <- e.info$regM.v
	
	add <- which(regM==1)
	ipv <- regM.v[add]
	
	e.vec <- matrix(0,length(e.info$eid.stend[,1]),3)
	for(i in 1:length(e.vec[,1])){
		e.vec[i,] <- X3[e.info$eid.stend[i,2],] - X3[e.info$eid.stend[i,1],]
	}
	M <- e.vec %*% t(e.vec)
	# 辺の長さの平均を1にそろえる
	scaler <- mean(diag(M))
	X3. <- X3/sqrt(scaler)
	X3. <- t(t(X3.) - apply(X3.,2,mean))
	M. <- M/scaler
	
	K <- M.
	n <- length(M.[,1])
	ret <- list()
	ret[[1]] <- K
	tmp.rsvd <- rsvd(K,k=k)
	edges <- diag(sqrt(tmp.rsvd$d)) %*% t(tmp.rsvd$u)
	retX <- list()
	retX[[1]] <- edges
	iter.cnt <- 1
	while(iter.cnt < maxiter){
		# 内積制約
	  check.ip <- max(abs(K * regM - regM.v))
	  #check.rowsum <- max(abs(apply(K,1,sum)))
	  if(check.ip < eps){
	    return(list(edges = retX, IPmat = ret))
	    break
	  }
	  # 内積制約の調整
	  K[add] <- ipv
	  # rsvd
	  tmp.rsvd <- rsvd(K,k=k)
	  K <- tmp.rsvd$v %*% diag(tmp.rsvd$d) %*% t(tmp.rsvd$u)
		ret[[iter.cnt+1]] <- K
	  edges <- diag(sqrt(tmp.rsvd$d)) %*% t(tmp.rsvd$u)
	  retX[[iter.cnt+1]] <- edges
		
		iter.cnt <- iter.cnt + 1
	}
	# 内積制約の調整
	  K[add] <- ipv
	  # rsvd
	  tmp.rsvd <- rsvd(K,k=k)
	  K <- tmp.rsvd$v %*% diag(tmp.rsvd$d) %*% t(tmp.rsvd$u)
		ret[[iter.cnt+1]] <- K
	  edges <- diag(sqrt(tmp.rsvd$d)) %*% t(tmp.rsvd$u)
	  ret[[iter.cnt+1]] <- K
	  retX[[iter.cnt+1]] <- edges
	#edges <- diag(sqrt(tmp$d)) %*% t(tmp$u)
	final.K <- ret[[iter.cnt+1]]
	final.edges <- retX[[iter.cnt+1]]
	
	nv <- max(f)
	Xn <- diag(rep(1,nv))
	En <- t(e.info$eid.EV)
	
	final.edges. <- matrix(0,ne/2,ne/2)
	final.edges.[1:k,] <- final.edges
	R <- final.edges. %*% Ginv(En)
	
	X3.est <- R %*% Xn
	X3.est. <- t(X3.est[1:k,])
	return(list(edges = final.edges,IPmat=final.K,X3.est = X3.est., X3.estmat = X3.est,R=R,edges.hx = retX, IPmat.hx = ret,e.info=e.info,niter=iter.cnt,X3.ori.st=X3.,scaler=scaler))
}

icosa <- rbind(c(1,2,3),c(1,3,4),c(1,4,5),c(1,5,6),c(1,6,2),c(3,2,7),c(4,3,8),c(5,4,9),c(6,5,10),c(2,6,11),c(3,7,8),c(4,8,9),c(5,9,10),c(6,10,11),c(2,11,7),c(12,8,7),c(12,9,8),c(12,10,9),c(12,11,10),c(12,7,11))

f <- icosa
theta1 <- (1:5)/5 * 2 * pi
theta2 <- theta1 + theta1[1]/2
x <- rbind(c(0,0,4),cbind(cos(theta1),sin(theta1),rep(3,5)),cbind(cos(theta2),sin(theta2),rep(2,5)),c(0,0,1)) * sqrt(2)

x <- x + rnorm(length(x),0,0.2)



out <- my.rsvd.tri(x,f)
el <- out$e.info$e2.stend

tmp <- c(out$X3.ori.st,out$X3.est)
plot3d(rbind(rep(min(tmp),3),rep(max(tmp),3)))

#plot3d(out$X3.ori.st)
spheres3d(out$X3.ori.st,col=1,radius=0.05)
segments3d(out$X3.ori.st[c(t(el)),])
spheres3d(out$X3.est,col=2,radius=0.1)
segments3d(out$X3.est[c(t(el)),])

plot3d(rbind(rep(min(tmp),3),rep(max(tmp),3)))

#plot3d(out$X3.ori.st)
#spheres3d(out$X3.ori.st,col=1,radius=0.05)
#segments3d(out$X3.ori.st[c(t(el)),])
spheres3d(out$X3.est,col=2,radius=0.1)
segments3d(out$X3.est[c(t(el)),])


plot(out$IPmat[which(out$e.info$regM==1)])

