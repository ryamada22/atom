# fは3列の行列
# 行数は三角形数(2k)
# 列には、反時計回りに構成頂点id(1,2,...,(k+2))
# 返り値
# kは|V|=k+2,|E|=3k,|F|=2kのk
# fvは、行F、3列でvid
# feは、行F、3列でeid(このeは向き考慮、全部で6k本のエッジ)
# efは、長さ6kのベクトル。ベクトルの番地はエッジId。反時計回りで考えて、帰属するfidを値として持つ
# evは、6k行、2列の行列。行はeid、第1列は始点vid、第２列は終点vid
# vfは、リスト、[[i]]番要素はi-th頂点が含まれるfid
# veは、リスト、[[i]]番要素はi-th頂点が接続するeid
# eeは、6k本の向きアリエッジの逆向き対応関係を表すベクトル
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

# my.planar.triangulation()の出力オブジェクトと、eidを指定して、変異の際に考慮するべき、v,e,fを取り出す
# e1,e2は、選ばれて付け替えれあれるeid(両方向に対応して２つのid)
# f1,f2は、対角線を付け替える四角形を構成する２つの三角形のfid
# fe1,fe2はf1,f2の周回エッジid3つずつ
# v1.v2はe1,e2の始点・終点vid
# loc... pre...,post...などは、処理用のオブジェクト
# pre.e1,pre.e2,post.e1,post.e2は、f1,f2の３つのエッジのうち、選ばれたエッジe1,e2と接続するエッジid。向きがあるので、pre,postの区別がある
# new.v1,new.v2は引き直されるエッジの両端vid
# new.ev1,new.ev2は、引き直されるエッジの両方向に対応する始点・終点vid
# new.fv1,new.fv2は、新たに作られる三角形の周囲vid(反時計回り)
# new.fe1,new.fe2は、新たに作られる三角形の周囲eid(反時計回り)
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
