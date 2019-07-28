
# 頂点IDを1,2,3,...,k+2,としたときに
# 三角形メッシュの各三角形の頂点IDを時計回りに並べたものを行とした
# 3列行列 fを引数とし、
# その双対グラフの４通りの行列WR4R, WR4L, WL4R,WL4Lを返す関数


# f is a 3-column matrix
# Each row of f has 3 node ids in clockwise rotation
# node id should start from 1 and end with the number of nodes
# L.sort is the option to make the edges in two direction of the same edge should have their id numbers as 2i-1 and 2i

# Wc4c.cpx = WR4R, Wc4d.cpx = WR4L, Wd4c.cpx = WL4R, Wd4d.cpx = WL4L

my.W.cpx <- function(f,L.sort=TRUE){
  F <- matrix(1:(length(f[,1])*3),byrow=TRUE,ncol=3)
  EL <- my.EfromF(f)
  E.pair <- my.Epair(EL)
  RL <- my.permMat(F,EL,E.pair)
  # RL出力から、双対グラフの時計回り接続行列を作る
  Wc <- RL$L %*% RL$R
# 双対グラフの反時計回り接続行列を作る
  Wd <- RL$L %*% t(Wc) %*% RL$L
# 行・列入れ替え
  s <- sample(1:length(Wc[,1]))
  S <- diag(length(s))[,s]
  Wc. <- S %*% Wc %*% t(S)
  Wd. <- S %*% Wd %*% t(S)
# 双対グラフのゼータ関数向けの接続行列を作る
  We <- Wc + Wd
  We. <- Wc. + Wd. 
  
  # 複素行列
  # このWc.cpx,Wd.cpxは、時計回りに部分を切り取るときに、時計回りエッジを選択する場合と
  # 反時計回りに部分を切り取るときに、反時計回りエッジを選択する場合
  ret.steps <- my.return.step(Wc)
  Wc.cpx <- my.complex.W(Wc,counterclockwise=FALSE)
  #Wd.cpx <- my.complex.W(Wd,counterclockwise=FALSE)
  Wd.cpx <- RL$L %*% t(Conj(Wc.cpx)) %*% RL$L
  
  opp.cpx <- my.complex.W.opp(Wc.cpx,Wd.cpx,RL$L)
  
  # 4行列を追合関係の２行列にする
  ZeroMat <- matrix(0,length(Wc.cpx[,1]),length(Wc.cpx[,1]))
  W <- rbind(cbind(Wc.cpx,ZeroMat),cbind(ZeroMat,opp.cpx$Wc4d))
  W.star <- rbind(cbind(opp.cpx$Wd4c,ZeroMat),cbind(ZeroMat,Wd.cpx))
  if(L.sort){
    SS <- my.L.sort(RL$L)
    WR4R=t(SS) %*% Wc.cpx %*% SS
    WL4L=t(SS) %*% Wd.cpx %*% SS
    WR4L = t(SS) %*% opp.cpx$Wc4d %*% SS
    WL4R = t(SS) %*% opp.cpx$Wd4c %*% SS
    # 4行列を追合関係の２行列にする
  ZeroMat <- matrix(0,length(Wc.cpx[,1]),length(Wc.cpx[,1]))
  W <- rbind(cbind(WR4R,ZeroMat),cbind(ZeroMat,WR4L))
  W.star <- rbind(cbind(WL4R,ZeroMat),cbind(ZeroMat,WL4L))
    return(list(Wc=t(SS) %*% Wc %*% SS,Wd=t(SS) %*% Wd %*% SS,We=t(SS) %*% We %*% SS,Wc.=t(SS) %*% Wc. %*% SS, Wd.=t(SS) %*% Wd. %*% SS,We.=t(SS) %*% We. %*% SS,WR4R=WR4R, WL4L=WL4L, WR4L=WR4L,WL4R = WL4R,W = W, W.star = W.star , R = t(SS) %*% RL$R %*% SS , L = t(SS) %*% RL$L %*% SS,F=F,EL=EL,E.pair=E.pair))
  }else{
    return(list(Wc=Wc,Wd=Wd,We=We,Wc.=Wc.,Wd.=Wd.,We.=We.,WR4R=Wc.cpx,WL4L=Wd.cpx,WR4L = opp.cpx$Wc4d, WL4R = opp.cpx$Wd4c, W=W, W.star = W.star,R = RL$R, L = RL$L,F=F,EL=EL,E.pair=E.pair))
  }
  
}


# 雑な関数ではあるが
# 行列のn乗を返す関数

# power of Matrix
my.powM <- function(M,n){
	d <- length(M[,1])
	ret <- diag(d)
	for(i in 1:n){
		ret <- ret %*% M
	}
	return(ret)
}
# 固有値 -1 については、-1 + 0iと-1 - 0iとが同数になるように一工夫してある
my.powM2 <- function(M,n){
  eigen.out <- eigen(M)
  eval.WcWd <- eigen.out[[1]]
  close2one <- which(Arg(eval.WcWd)/pi > 1- (10)^(-10))
  eval.WcWd[close2one] <- rep(c(exp(1i * (-1)*pi),exp(1i * pi)),length(close2one)/2)
  eigen.out[[1]] <- eval.WcWd
  #eigen.out <- eigen(M)
  tmp <- solve(eigen.out[[2]])
  eigen.out[[2]] %*% diag(eigen.out[[1]]^n) %*% tmp
  
}
# 行列のトレースを返す
# Trace of matrix
my.trace <- function(M){
	sum(diag(M))
}

# 適当に三角メッシュ平面グラフの情報を返す
# fは、頂点のトリオを時計回り(反時計回りと見なしても良い)に並べた行からなる３列行列
# Fの各行に向きのある辺の番号を与えたものとしてい作った３列行列

# generates a random triangular planar graph information
# F is a matrix with 3 columns; each row has three integers that represent three edges that are in the order of clockwise rotation (or you can consider the order as counterclockwise.)
# Returned values = F and f
# F's row is just a three numbers of each triangle.
# f is a matrix whose each row is consisted of three vertex ids in clockwise rotation order.
my.rtri <- function(n){
	f <- matrix(1:3,ncol=3)
	f <- rbind(f,3:1)
	for(i in 1:n){
		s <- sample(1:length(f[,1]),1)
		v <- max(f) + 1
		sf <- f[s,]
		f <- f[-s,]
		f <- rbind(f,c(sf[1:2],v))
		f <- rbind(f,c(sf[2:3],v))
		f <- rbind(f,c(sf[c(3,1)],v))
	}
	F <- matrix(1:(length(f[,1])*3),byrow=TRUE,ncol=3)
	return(list(f=f,F=F))
}
# すべての3角形をn回ずつ分割する
my.div.tri <- function(f,n){
  el <- rbind(f[,1:2],f[,2:3],f[,c(3,1)])
  max.v <- max(c(el))
  el. <- t(apply(el,1,sort))
  el.v <- el.[,1] + (max.v+1) * el.[,2]
  ord <- rank(el.v,ties.method="max")/2 + max.v # 該当エッジに発生させる新頂点のID番号
  for(i in 1:n){
    len <- length(f[,1])
    for(j in 1:len){
      v1 <- ord[j]
      v2 <- ord[len+j]
      v3 <- ord[2*len+j]
      
      sf <- f[1,]
      f <- f[-1,]
      f <- rbind(f,c(sf[1],v1,v3))
		  f <- rbind(f,c(sf[2],v2,v1))
		  f <- rbind(f,c(sf[3],v3,v2))
		  f <- rbind(f,c(v1,v2,v3))
    }
  }
  F <- matrix(1:(length(f[,1])*3),byrow=TRUE,ncol=3)
	return(list(f=f,F=F))
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

# 何歩で戻るかを列挙する関数
my.return.step <- function(Wc){
  returned <- list()
  for(i in 1:length(Wc[,1])){
    returned[[i]] <- 0
  }
  for(i in 1:length(Wc[,1])){
    tmp <- diag(my.powM(Wc,i))
    if(sum(tmp)>0){
      for(j in which(tmp==1)){
        returned[[j]] <- c(returned[[j]],i)
      }    
    }
  }
  returned. <- lapply(returned,function(x){x[-1]})
  returned.  
}
# 何歩で戻るかの情報で、Wc行列を複素行列化
# この関数は、時計回りに部分を切り取るときに、時計回りエッジを選択する場合と
# 反時計回りに部分を切り取るときに、反時計回りエッジを選択する場合
my.complex.W <- function(Wc,counterclockwise=TRUE){
  returned <- my.return.step(Wc)
  tmp <- sapply(returned,min)
  ret <- Wc
  pm <- 1
  if(!counterclockwise){
    pm <- -1
  }
  for(i in 1:length(tmp)){
    val <- exp(pm * 1i * 1/tmp[i]*2*pi)
    loc <- which(ret[i,]==1)
    ret[i,loc] <- val 
  }
  return(ret)
}
# この関数は、時計回りに部分を切り取るときに、反時計回りエッジと選ぶ場合と
# 反時計回りに部分を切り取るときに、時計回りエッジを選ぶ場合
# Wc,Wdは時計回り・反時計回りの複素行列,Pはエッジ対応行列
my.complex.W.opp <- function(Wc,Wd,P){
  Wd4c <- Wd # clockwiseに切り取るときのdエッジ選択
  for(i in 1:length(Wc[,1])){
    tmpd <- which(Wd[i,] != 0)
    tmpc <- which(Wc[i,] != 0)
    counter <- which(P[tmpc,] != 0)
    tmpc2 <- which(Wc[counter,] != 0)
    #print(Wc[i,tmpc])
    #print(Wc[counter,tmpc2])
    #print("0000")
    Wd4c[i,tmpd] <- Wc[i,tmpc] * Wc[counter,tmpc2] / exp(1i * pi)
  }
  Wc4d <- P %*% t(Conj(Wd4c)) %*% t(P)
  return(list(Wd4c=Wd4c,Wc4d=Wc4d))
}

# 有向辺のペアをあらわした行列Lを引数にして、(1,2),(2,3),...とそろえるために、エッジの番号付けを変換するための置換行列を返す関数
my.L.sort <- function(L){
	pairs <- which(L==1,arr.ind=TRUE)
	pairs <- pairs[which(pairs[,1]<pairs[,2]),]
	sfl <- c(t(pairs))
	SS <- diag(length(sfl))[,sfl]
	return(SS) # t(SS) %*% L %*% SS is the sorting calculation
}

# ランダムに三角形planar graphを作り、それに関する辺接続行列等を返す関数
my.rTriRL <- function(n.tri,L.sort=TRUE){
  F <- my.rtri(n.tri)
  EL <- my.EfromF(F$f)
  E.pair <- my.Epair(EL)
  RL <- my.permMat(F$F,EL,E.pair)
  # RL出力から、双対グラフの時計回り接続行列を作る
  Wc <- RL$L %*% RL$R
# 双対グラフの反時計回り接続行列を作る
  Wd <- RL$L %*% t(Wc) %*% RL$L
# 行・列入れ替え
  s <- sample(1:length(Wc[,1]))
  S <- diag(length(s))[,s]
  Wc. <- S %*% Wc %*% t(S)
  Wd. <- S %*% Wd %*% t(S)
# 双対グラフのゼータ関数向けの接続行列を作る
  We <- Wc + Wd
  We. <- Wc. + Wd. 
  
  # 複素行列
  # このWc.cpx,Wd.cpxは、時計回りに部分を切り取るときに、時計回りエッジを選択する場合と
  # 反時計回りに部分を切り取るときに、反時計回りエッジを選択する場合
  ret.steps <- my.return.step(Wc)
  Wc.cpx <- my.complex.W(Wc,counterclockwise=FALSE)
  #Wd.cpx <- my.complex.W(Wd,counterclockwise=FALSE)
  Wd.cpx <- RL$L %*% t(Conj(Wc.cpx)) %*% RL$L
  
  opp.cpx <- my.complex.W.opp(Wc.cpx,Wd.cpx,RL$L)
  
  # 4行列を追合関係の２行列にする
  ZeroMat <- matrix(0,length(Wc.cpx[,1]),length(Wc.cpx[,1]))
  W <- rbind(cbind(Wc.cpx,ZeroMat),cbind(ZeroMat,opp.cpx$Wc4d))
  W.star <- rbind(cbind(opp.cpx$Wd4c,ZeroMat),cbind(ZeroMat,Wd.cpx))
  
  if(L.sort){
    SS <- my.L.sort(RL$L)
    Wc4c.cpx=t(SS) %*% Wc.cpx %*% SS
    Wd4c.cpx=t(SS) %*% Wd.cpx %*% SS
    Wc4d.cpx = t(SS) %*% opp.cpx$Wc4d %*% SS
    Wd4c.cpx = t(SS) %*% opp.cpx$Wd4c %*% SS
    
    # 4行列を追合関係の２行列にする
    W <- rbind(cbind(Wc4c.cpx,ZeroMat),cbind(ZeroMat,Wc4d.cpx))
  W.star <- rbind(cbind(Wd4c.cpx,ZeroMat),cbind(ZeroMat,Wd4c.cpx))
    
    return(list(Wc=t(SS) %*% Wc %*% SS,Wd=t(SS) %*% Wd %*% SS,We=t(SS) %*% We %*% SS,Wc.=t(SS) %*% Wc. %*% SS, Wd.=t(SS) %*% Wd. %*% SS,We.=t(SS) %*% We. %*% SS,Wc4c.cpx=Wc4d.cpx ,Wd4d.cpx=Wd4d.cpx,Wc4d.cpx = Wc4d.cpx, Wd4c.cpx = Wdec.cpx, W = W, W.star = W.star,R = t(SS) %*% RL$R %*% SS , L = t(SS) %*% RL$L %*% SS,F=F,EL=EL,E.pair=E.pair))
  }else{
    return(list(Wc=Wc,Wd=Wd,We=We,Wc.=Wc.,Wd.=Wd.,We.=We.,Wc4c.cpx=Wc.cpx,Wd4d.cpx=Wd.cpx,Wc4d.cpx = opp.cpx$Wc4d, Wd4c.cpx = opp.cpx$Wd4c, W=W,W.star=W.star,R = RL$R, L = RL$L,F=F,EL=EL,E.pair=E.pair))
  }
  
}

n.tri <- 3
F1 <- my.rtri(n.tri)
f1 <- F1$f

out1 <- my.W.cpx(f1)

W <- out1$We
Wvec <- c(W)
n <- length(W[,1])
