# 3列の行列 f : 行が三角形、列がその頂点の反時計回り

my.planar.triangulation <- function(f){
	# kは、|V|=k+2, |E|=3k, |F|=2k
	k <- length(f[,1])/2
	# fごとに構成vを反時計回りに
	fv <- f 
	# e ごとに有向エッジの始点と終点
	ev <- rbind(f[,1:2],f[,2:3],f[,c(3,1)])
	# eの向き違いペア
	ev.tmp <- t(apply(ev,1,sort))
	tmp.d <- as.matrix(dist(ev.tmp))
	ee <- rep(0,6*k)
	for(i in 1:(6*k)){
		tmp <- which(tmp.d[i,]==0)
		if(tmp[1]==i){
			ee[i] <- tmp[2]
		}else{
			ee[i] <- tmp[1]
		}
	}
	# fごとに構成eを反時計回りに
	fe <- matrix(1:(6*k),ncol=3)
	# eが帰属するfはひとつ。ペアeを取ると２つのfと対応
	ef <- rep(1:(2*k),3)
	# ペアeのfは、ef[ee] で得られる
	# vf と veはリスト,f,eの順序情報はなし
	vf <- ve <- list()
	for(i in 1:(k+2)){
		vf[[i]] <- which(f==i,arr.ind=TRUE)[,1]
		ve[[i]] <- which(ev==i,arr.ind=TRUE)[,1]
	}
	return(list(k=k,fv=fv,fe=fe,ef=ef,ev=ev,vf=vf,ve=ve,ee=ee))
}

# quiver mutationはエッジ(双方向)の付け替え
# 付け替えにあたり、ev/ve,ef/fe,vf/fvが変わる。k,eeは変わらない

# eを指定して、関係するv,e,fを取り出す
my.planar.mut.invlv <- function(pt,e){
	e1 <- e
	e2 <- pt$ee[e1]
	f1 <- pt$ef[e1]
	f2 <- pt$ef[e2]
	fe1 <- pt$fe[f1,]
	fe2 <- pt$fe[f2,]
	v1 <- pt$ev[e1,]
	v2 <- pt$ev[e2,]
	
	loc.e1 <- which(fe1==e1)
	loc.e2 <- which(fe2==e2)
	
	pre.loc.e1 <- (c(3,1,2))[loc.e1]
	pre.loc.e2 <- (c(3,1,2))[loc.e2]
	post.loc.e1 <- (c(2,3,1))[loc.e1]
	post.loc.e2 <- (c(2,3,1))[loc.e2]
	
	pre.e1 <- fe1[pre.loc.e1]
	pre.e2 <- fe2[pre.loc.e2]
	post.e1 <- fe1[post.loc.e1]
	post.e2 <- fe2[post.loc.e2]
	
	tri.v1 <- pt$fv[f1,]
	tri.v2 <- pt$fv[f2,]
	
	new.v1 <- tri.v1[which(! (tri.v1 %in% v1))]
	new.v2 <- tri.v2[which(! (tri.v2 %in% v2))]
	
	# e1 : v1 -> c(new.v1,new.v2), e2 : v2 -> c(new.v2,new.v1)
	new.ev1 <- c(new.v1,new.v2)
	new.ev2 <- new.ev1[2:1]
	
	
	new.fv1 <- c(new.v1,v1[1],new.v2)
	new.fv2 <- c(new.v1,new.v2,v1[2])
	
	new.fe1 <- c(pre.e1,post.e2,e2)
	new.fe2 <- c(pre.e2,post.e1,e1)
	
	return(list(e1=e1,e2=e2,f1=f1,f2=f2,fe1=fe1,fe2=fe2,v1=v1,v2=v2,loc.e1=loc.e1,loc.e2=loc.e2,pre.loc.e1=pre.loc.e1,pre.loc.e2=pre.loc.e2,post.loc.e1=post.loc.e1,post.loc.e2=post.loc.e2,pre.e1=pre.e1,pre.e2=pre.e2,post.e1=post.e1,post.e2=post.e2,tri.v1=tri.v1,tri.v2=tri.v2,new.v1=new.v1,new.v2=new.v2,new.ev1=new.ev1,new.ev2=new.ev2,new.fv1=new.fv1,new.fv2=new.fv2,new.fe1=new.fe1,new.fe2=new.fe2))
}

my.planar.mut.invlv.check <- function(pt.e){
	invlv <- my.planar.mut.invlv(pt,e)
	
}
my.planar.mut <- function(pt,e){
	invlv <- my.planar.mut.invlv(pt,e)
	e1 <- invlv$e1
	e2 <- invlv$e2
	f1 <- invlv$f1
	f2 <- invlv$f2
	fe1 <- invlv$fe1
	fe2 <- invlv$fe2
	v1 <- invlv$v1
	v2 <- invlv$v2
	
	loc.e1 <- invlv$loc.e1
	loc.e2 <- invlv$loc.e2
	
	pre.loc.e1 <- invlv$pre.loc.e1
	pre.loc.e2 <- invlv$pre.loc.e2
	post.loc.e1 <- invlv$post.loc.e1
	post.loc.e2 <- invlv$post.loc.e2
	
	pre.e1 <- invlv$pre.e1
	pre.e2 <- invlv$pre.e2
	post.e1 <- invlv$post.e1
	post.e2 <- invlv$post.e2
	
	tri.v1 <- invlv$tri.v1
	tri.v2 <- invlv$tri.v2
	
	new.v1 <- invlv$new.v1
	new.v2 <- invlv$new.v2
	
	# e1 : v1 -> c(new.v1,new.v2), e2 : v2 -> c(new.v2,new.v1)
	new.ev1 <- invlv$new.ev1
	new.ev2 <- invlv$new.ev2
	
	
	new.fv1 <- invlv$new.fv1
	new.fv2 <- invlv$new.fv2
	
	new.fe1 <- invlv$new.fe1
	new.fe2 <- invlv$new.fe2
	
	pt$ef[e1] <- f2
	pt$ef[e2] <- f1
	pt$ef[post.e2] <- f1
	pt$ef[post.e1] <- f2
	
	pt$ev[e1,] <- new.ev1
	pt$ev[e2,] <- new.ev2
	
	pt$fe[f1,] <- new.fe1
	pt$fe[f2,] <- new.fe2
	
	pt$fv[f1,] <- new.fv1
	pt$fv[f2,] <- new.fv2
	
	#pt$vf[[new.v2]][which(pt$vf[[new.v2]]==f2)] <- f1
	pt$vf[[new.v1]] <- c(pt$vf[[new.v1]],f2)
	#pt$vf[[new.v1]][which(pt$vf[[new.v1]]==f1)] <- f2
	pt$vf[[new.v2]] <- c(pt$vf[[new.v2]],f1)
	
	pt$vf[[v1[1]]] <- pt$vf[[v1[1]]][-which(pt$vf[[v1[1]]] == f2)]
	pt$vf[[v1[2]]] <- pt$vf[[v1[2]]][-which(pt$vf[[v1[2]]] == f1)]
	
	pt$ve[[new.v1]] <- c(pt$ve[[new.v1]],c(e1,e2))
	pt$ve[[new.v2]] <- c(pt$ve[[new.v2]],c(e1,e2))
	pt$ve[[v1[1]]] <- pt$ve[[v1[1]]][-which(pt$ve[[v1[1]]] %in% c(e1,e2))]
	pt$ve[[v1[2]]] <- pt$ve[[v1[2]]][-which(pt$ve[[v1[2]]] %in% c(e1,e2))]
	
	return(pt)
}
# ６面体
hex.tri <- rbind(c(1,2,3),c(1,3,4),c(1,4,2),c(5,3,2),c(5,4,3),c(5,2,4))
hex.pla <- my.planar.triangulation(hex.tri)
pt <- hex.pla
e <- 8

mut.pt2 <- my.planar.mut(pt,e)


# サッカーボール対応三角メッシュ
football.tri <- rbind(c(1,7,4),c(4,7,8),c(1,9,7),c(1,2,9),c(2,23,9),c(2,3,23),c(3,29,23),c(3,19,29),c(3,13,19),c(2,13,3),c(2,10,13),c(1,10,2),c(1,5,10),c(1,4,5),c(4,6,5),c(4,8,6),c(6,8,17),c(8,26,17),c(17,26,31),c(17,31,16),c(6,17,16),c(11,6,16),c(5,6,11),c(5,11,12),c(10,5,12),c(10,12,13),c(13,12,14),c(13,14,19),c(8,21,26),c(8,7,21),c(7,22,21),c(7,9,22),c(9,23,22),c(22,23,25),c(23,29,25),c(16,31,32),c(16,32,15),c(16,15,11),c(11,15,12),c(12,15,14),c(15,32,18),c(15,18,14),c(14,18,19),c(20,19,18),c(29,19,20),c(26,27,31),c(27,26,24),c(26,21,24),c(21,22,24),c(24,22,25),c(27,24,28),c(24,25,28),c(25,29,28),c(28,29,20),c(31,27,30),c(27,28,30),c(28,20,30),c(32,31,30),c(18,32,30),c(20,18,30))

foot.pla <- my.planar.triangulation(football.tri)
str(foot.pla)

mut.foot <- my.planar.mut(foot.pla,5)