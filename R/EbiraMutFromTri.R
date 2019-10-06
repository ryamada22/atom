my.EbiraB.from.tri <- function(h){
	tmp <- my.planar.triangulation(h)
	ee <- cbind(1:length(tmp$ee),tmp$ee)
	ee. <- t(apply(ee,1,sort))
	ee.unique <- unique(ee.)
	# edge idを6k有向エッジ番号から呼び出せるようにする
	ee.id.pair <- cbind(rep(1:length(ee.unique[,1]),2),c(ee.unique))
	ee.id <- ee.id.pair[order(ee.id.pair[,2]),1]
	fe <- tmp$fe
	fe. <- fe
	for(i in 1:length(ee.unique[,1])){
		tmp2 <- which(fe.==ee.unique[i,2])
		fe.[tmp2] <- ee.unique[i,1]
	}
	ret <- matrix(0,tmp$k*3,tmp$k*3)
	for(i in 1:length(fe.[,1])){
		ret[ee.id[fe.[i,1]],ee.id[fe.[i,2]]] <- ret[ee.id[fe.[i,2]],ee.id[fe.[i,3]]] <- ret[ee.id[fe.[i,3]],ee.id[fe.[i,1]]] <- 1
		ret[ee.id[fe.[i,2]],ee.id[fe.[i,1]]] <- ret[ee.id[fe.[i,3]],ee.id[fe.[i,2]]] <- ret[ee.id[fe.[i,1]],ee.id[fe.[i,3]]] <- -1
	}
	tmp$ee.id <- ee.id
	tmp$ee.unique <- ee.unique
	return(list(graph.info=tmp,B=ret))
}

B.hexa <- my.EbiraB.from.tri(hexa)
library(Ryacas)
#library(Ryacas0) # remotes::install_github("mikldk/ryacas0")

my.XY.mut.fromB <- function(B){
	n <- length(B[,1])
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
	return(list(x.mut=Ebira.x,y.mut=Ebira.y))
}

my.XY.mut.from.tri <- function(h){
	B.h <- my.EbiraB.from.tri(h)
	XY.mut <- my.XY.mut.fromB(B.h$B)
	
	return(list(graph.info=B.h$graph.info,B=B.h$B,XY.mut=XY.mut))
}

XYmut.hexa <- my.XY.mut.from.tri(hexa)
XY.mut.octa <- my.XY.mut.from.tri(octa)

