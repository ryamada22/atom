my.non.negative <- function(a){
	ret <- max(a,0)
	return(ret)
}
my.x.mut <- function(X,B,k=1:length(B[,1])){
	ret <- list()
	n <- length(B[,1])
	nk <- length(k)
	for(i in 1:nk){
		ret[[i]] <- X
		tmp1 <- 1
		tmp2 <- 1
		for(j in 1:n){
			if(j != k[i]){
				b1 <- B[j,k[i]]
				b2 <- B[k[i],j]
				if(b1 > b2){
					tmp1 <- tmp1 * X[[j]]^my.non.negative(b1)
					tmp2 <- tmp2 * X[[j]]^my.non.negative(b2)
				}else{
					tmp1 <- tmp1 * X[[j]]^my.non.negative(b2)
					tmp2 <- tmp2 * X[[j]]^my.non.negative(b1)
				}
				
			}
		}
		tmp <- 1/X[[k[i]]] * (tmp1 + tmp2)
		ret[[i]][[k[i]]] <- tmp
		ret[[i]][[k[i]]] <- Simplify(ret[[i]][[k[i]]])
	}
	if(nk==1){
		return(ret[[1]])
	}else{
		return(ret)
	}
}

my.y.mut <- function(Y,B,k=1:length(B[,1])){
	ret <- list()
	n <- length(B[,1])
	nk <- length(k)
	for(i in 1:nk){
		ret[[i]] <- Y
		Yk <- Y[[k[i]]]
		for(j in 1:n){
			if(j == k[i]){
				ret[[i]][[j]] <- 1/Yk
			}else{
				b <- B[j,k[i]]
				ret[[i]][[j]] <- Y[[j]] * (1+Yk^sign(b))^b
			}
		}
		ret[[i]][[j]] <- Simplify(ret[[i]][[j]])
	}
	if(nk==1){
		return(ret[[1]])
	}else{
		return(ret)
	}

}
my.Ebira.mut <- function(B,ks=1:length(B[,1])){
  ret <- list()
  n <- length(ks)
  for(k in 1:n){
    ret[[k]] <- matrix(0,n,n)
    for(i in 1:n){
      for(j in 1:n){
        if(i == ks[k] | j == ks[k]){
          ret[[k]][i,j] <- (-1) * B[i,j]
        }else{
          tmp1 <- B[i,ks[k]]
          tmp2 <- -B[ks[k],j]
          if(tmp1 <= 0){
            tmp1 <- 0
          }
          if(tmp2 <= 0){
            tmp2 <- 0
          }
          ret[[k]][i,j] <- B[i,j] + tmp1 * B[ks[k],j] + tmp2 * B[i,ks[k]]
        }
      }
      
    }
  }
  if(n == 1){
    return(ret[[1]])
  }
  return(ret)
}
# PDFの箙変換則をRのベクトル演算に合わせて記載
my.Ebira.mut2 <- function(B,ks=1:length(B[,1])){
  ret <- list()
  n <- length(ks)
  for(k in 1:n){
    ret[[k]] <- B
    # 変わるのはk行・k列
    out.K <- B[ks[k],]
    in.K <- B[,ks[k]]
    # 非負のみを問題にする
    in.K.nn <- in.K * (in.K >= 0)
    out.K.nn <- out.K * (out.K >= 0)
    prod.K <- matrix(in.K.nn,ncol=1) %*% matrix(out.K.nn,nrow=1)
    
    ret[[k]] <- ret[[k]] + prod.K + (-1) * t(prod.K)
    
    ret[[k]][,ks[k]] <- (-1) * in.K
    ret[[k]][ks[k],] <- (-1) * out.K
  }
  if(n == 1){
    return(ret[[1]])
  }
  return(ret)
}
# Xの双対変数を作る
my.YfromX <- function(X,B){
	n <- length(X)
	Y <- list()
	for(i in 1:n){
		Y[[i]] <- 1
		for(j in 1:n){
			if(j != i){
				Y[[i]] <- Y[[i]] * X[[j]]^B[j,i]
			}
		}
	}
	Y
}


n <- 5
B <- matrix(sample(0:6,n^2,replace=TRUE),n,n)
B <- B + (-1) * t(B)

# Symbolic calculation
library(Ryacas)
library(Ryacas0) # remotes::install_github("mikldk/raycas0")

X <- list()
Y <- list()
for(i in 1:n){
	tmp <- paste("x",i,sep="")
	X[[i]] <- Sym(tmp)
	tmp <- paste("y",i,sep="")
	Y[[i]] <- Sym(tmp)

}

Ebira.x <- my.x.mut(X,B)
Ebira.y <- my.y.mut(Y,B)

# n=2で、変異のために選択するノードを交互にする
# 交互にしないとき＝２連続で同じノードで変異すると、それは元に戻る性質を持っている


X <- list()
for(i in 1:n){
	tmp <- paste("x",i,sep="")
	X[[i]] <- Sym(tmp)
}

B <- matrix(c(0,-1,1,0),2,2)

ks <- rep(c(1,2),2)

x.series <- list()
x.series[[1]] <- X
print(x.series[[1]])
for(i in 1:length(ks)){
	x.series[[i+1]] <- my.x.mut(x.series[[i]],B,k=ks[i])
	print(x.series[[i+1]])
}

# Y変数をx変数の双対にする

n <- 3
B <- matrix(sample(0:6,n^2,replace=TRUE),n,n)
B <- B + (-1) * t(B)


X <- list()
Y <- list()
for(i in 1:n){
	tmp <- paste("x",i,sep="")
	X[[i]] <- Sym(tmp)
}


Y <- my.YfromX(X,B)

Ebira.x <- my.x.mut(X,B,k=1)


B. <- my.Ebira.mut2(B)

Y. <- my.YfromX(Ebira.x,B.[[1]])


Ebira.y <- my.y.mut(Y,B,k=1)

# どの項も同一だが、Simplify()が不完全なので
# そのように見えにくいかも
for(i in 1:n){
	print(Simplify(Y.[[i]]/Ebira.y[[i]]))
}

## 向き付け三角メッシュから反対称行列を作る

my.EbiraB.tri <- function(h){
	k <- length(h[,1])/2 * 3
	#n <- max(h)
	ret <- matrix(0,k,k)
	for(i in 1:length(h[,1])){
		ret[h[i,1],h[i,2]] <- ret[h[i,2],h[i,3]] <- ret[h[i,3],h[i,1]] <- 1
		ret[h[i,2],h[i,1]] <- ret[h[i,3],h[i,2]] <- ret[h[i,1],h[i,3]] <- -1
	}
	return(ret)
}