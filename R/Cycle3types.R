library(igraph)
# 向きの違うエッジのペアを、my.EfromF()関数の出力での行番をエッジＩＤと見なして返す
# The return is two-column matrix and each row is the number of edges paring in mutually opposite direction.
my.Epair <- function(el){
	s.el <- apply(el,1,sort)
	mv <- max(el)+1
	V <- s.el[1,] * mv + s.el[2,]
	tmp <- outer(V,V,"-")
	diag(tmp) <- 1
	pairs <- which(tmp==0,arr.ind=TRUE)
	return(pairs)
}
# 三角メッシュの各面の接続に相当する行列Ｒと
# 逆向きエッジ対応に相当する行列Ｌとを返す

# This function returns two square matrices, R and L
# R %*% L is the edge-connection matrix of clockwise rotation
my.permMat <- function(F,EL,E.pair){
	n <- length(E.pair[,1])
	R <- matrix(0,n,n)
	trio.mat <- matrix(0,3,3)
	trio.mat[1,2] <- trio.mat[2,3] <- trio.mat[3,1] <- 1
	for(i in 1:length(F[,1])){
		R[F[i,],F[i,]] <- trio.mat
	}
	L <- matrix(0,n,n)
	for(i in 1:length(E.pair[,1])){
		el1 <- EL[E.pair[i,1],]
		el2 <- EL[E.pair[i,2],]
		tmp1 <- which(apply((t(EL) - el1)^2,2,sum)==0)
		tmp2 <- which(apply((t(EL) - el2)^2,2,sum)==0)
		
		L[tmp1,tmp2] <- 1
	}
	return(list(R=R,L=L))
}

# ３角形の頂点情報から、エッジリストを作る
# my.rtri()関数の出力の$fから作る
# エッジの並び順は、my.rtri()関数の出力の$Fとうまく合致するように並べ替えて出力してある
# The function makes an edgelist with direction
my.EfromF <- function(f){
	el <- rbind(f[,1:2],f[,2:3],f[,c(3,1)])
	tmp <- matrix(1:length(el[,1]),ncol=3)
	el <- el[c(t(tmp)),]
	return(el)
}
# 向き付け３角メッシュ(３列)、頂点IDを入力に、辺接続行列とその時計回り・反時計回り分解を作る
my.W <- function(f){
  F <- matrix(1:(length(f[,1])*3),byrow=TRUE,ncol=3)
  EL <- my.EfromF(f)
  E.pair <- my.Epair(EL)
  RL <- my.permMat(F,EL,E.pair)
  # RL出力から、双対グラフの時計回り接続行列を作る
  Wc <- RL$L %*% RL$R
# 双対グラフの反時計回り接続行列を作る
  Wd <- RL$L %*% t(Wc) %*% RL$L
	return(list(We = Wc + Wd,Wc = Wc, Wd =Wd, F=F, RL=RL,EL = EL, E.pair = E.pair))
}
# 置換行列から巡回の取り出し
my.permMat2cycles <- function(W){
	g <- graph.adjacency(t(W))
	cl <- clusters(g)
	n.cl <- cl$no
	cycles <- list()
	for(i in 1:n.cl){
		members <- which(cl$members==i)
		D <- distances(g,members[1],members,mode="out")
		cycles[[i]] <- members[order(D)]
	}
	return(cycles)
}
# 木下さんのサイクル取り出し
my.cycles.2 <- function(Wc,Wd){
	W <- Wd %*% Wc
	g <- graph.adjacency(t(W))
	cl <- clusters(g)
	n.cl <- cl$no
	two.step.cycles <- list()
	for(i in 1:n.cl){
		members <- which(cl$members==i)
		D <- distances(g,members[1],members,mode="out")
		two.step.cycles[[i]] <- members[order(D)]
	}
	Wc.neighbors <- apply(Wc,2,function(x){which(x==1)})
	one.step.cycles <- list()
	for(i in 1:n.cl){
		mem <- two.step.cycles[[i]]
		Wc.nb <- Wc.neighbors[mem]
		one.step.cycles[[i]] <- c(t(cbind(mem,Wc.nb)))
	}	
	return(one.step.cycles)
}
# 箙用直線サイクル
my.cycles.Quiver2 <- function(Wc,Wd){
	W <- Wd %*% Wc
	g <- graph.adjacency(t(W))
	cl <- clusters(g)
	n.cl <- cl$no
	two.step.cycles <- list()
	for(i in 1:n.cl){
		members <- which(cl$members==i)
		D <- distances(g,members[1],members,mode="out")
		two.step.cycles[[i]] <- members[order(D)]
	}
	Wc.neighbors <- apply(Wd,2,function(x){which(x==1)})
	one.step.cycles <- list()
	for(i in 1:n.cl){
		mem <- two.step.cycles[[i]]
		Wc.nb <- Wc.neighbors[mem]
		one.step.cycles[[i]] <- c(t(cbind(mem,Wc.nb)))
	}	
	#return(one.step.cycles)
	return(two.step.cycles)
}


# 向き付け３角メッシュ行列(３列)、頂点IDを入力に、箙を作る
my.F2Quiver <- function(f){
	# 以下のelは箙グラフの有向グラフのすべて
	el <- rbind(f[,1:2],f[,2:3],f[,c(3,1)])
	# 三角形サイクルと頂点周囲サイクルと擬直線サイクルとを取り出す
	tri.cycles <- matrix(1:(max(f)*3),byrow=TRUE,ncol=3)
	# 辺接続行列の時計回り・反時計回り
	Wout <- my.W(f)
	Wc <- Wout$Wc
	Wd <- Wout$Wd
	E.pairs <- Wout$E.pair
	v.cycles <- my.permMat2cycles(Wd)
	two.step.cycles <- my.cycles.Quiver2(Wd,Wc)
	one.step.cycles <- list()
	done <- rep(0,length(two.step.cycles))
	cnt <- 1
	for(i in 1:length(two.step.cycles)){
		if(done[i] ==1){
			next
		}
		this.e <- two.step.cycles[[i]][1]
		this.e.pair <- E.pairs[which(E.pairs[,1]==this.e),2]
		for(j in 1:length(two.step.cycles)){
			tmp <- which(two.step.cycles[[j]]==this.e.pair)
			print(tmp)
			if(length(tmp)>0){
				len <- length(two.step.cycles[[j]])
				if(tmp==1){
					ord <- c(1,len:2)
				}else if(tmp==len){
					ord <- len:1
				}else{
					ord <- c(tmp:1,len:(tmp+1))
				}
				print(ord)
				print(two.step.cycles[[j]][ord])
				tmp2 <- rbind(two.step.cycles[[j]][ord],two.step.cycles[[i]])
				one.step.cycles[[cnt]] <- c(tmp2)
				cnt <- cnt + 1
				done[i] <- done[j] <- 1
			}
		}
	}
	
	return(list(tri.cycles = tri.cycles,v.cycles=v.cycles,two.step=two.step.cycles,one.step.cycles=one.step.cycles,E.pairs=E.pairs))
}

my.F2Quiver(f)
