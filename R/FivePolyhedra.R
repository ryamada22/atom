library(permutations)
library(rgl)
# 各多面体の各面の置換
# エッジIDの値を1:78にするためにlapply(...,"+",...)を使っている
# 加える値は、p2,p3,,,に対して、s2,s3,...
s2 <- 12
s3 <- 12 + 18
s4 <- 12 + 2 * 18
s5 <- 12 + 3 * 18
s6 <- 2 * 12 + 3 * 18

p1.face <- cycle(list(list(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12))))

p5.face <- cycle(list(lapply(list(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12)),"+",s5)))

p2.face <- cycle(list(lapply(list(c(1,2,3),c(4,5,6,7),c(8,9,10,11),c(12,13,14),c(15,16,17,18)),"+",s2)))

p3.face <- cycle(list(lapply(list(c(1,2,3),c(4,5,6,7),c(8,9,10,11),c(12,13,14),c(15,16,17,18)),"+",s3)))

p4.face <- cycle(list(lapply(list(c(1,2,3),c(4,5,6,7),c(8,9,10,11),c(12,13,14),c(15,16,17,18)),"+",s4)))

p6.face <- cycle(list(lapply(list(c(1,2,3),c(4,5,6,7),c(8,9,10,11),c(12,13,14),c(15,16,17,18)), "+",s6)))

pall.face <- p1.face * p2.face * p3.face * p4.face * p5.face * p6.face

# 各多面体の各辺の置換

p1.edge <- cycle(list(list(c(1,4),c(2,10),c(3,7),c(5,9),c(6,11),c(8,12))))

p5.edge <- cycle(list(lapply(list(c(1,8),c(2,4),c(3,11),c(5,7),c(6,12),c(9,10)),"+",s5)))

p2.edge <- cycle(list(lapply(list(c(1,4),c(2,15),c(3,8),c(5,11),c(6,14),c(7,16),c(9,18),c(10,12),c(13,17)),"+",s2)))

p3.edge <- cycle(list(lapply(list( c(1,4),c(2,15),c(3,8),c(5,11),c(6,12),c(7,16),c(9,18),c(10,13),c(14,17) ),"+",s3)))

p4.edge <- cycle(list(lapply(list( c(1,15),c(2,8),c(3,4),c(5,11),c(6,14),c(7,16),c(9,18),c(10,12),c(13,17) ),"+",s4)))

p6.edge <- cycle(list(lapply(list(c(1,16),c(2,4),c(3,8),c(5,15),c(6,12),c(7,9),c(10,14),c(11,17),c(13,18)),"+",s6)))

pall.edge <- p1.edge * p2.edge * p3.edge * p4.edge * p5.edge * p6.edge

# 面合わせの置換

p.p1p2 <- cycle(list(list(c(4,3+s2),c(5,2+s2),c(6,1+s2))))
p.p1p3 <- cycle(list(list(c(7,3+s3),c(8,2+s3),c(9,1+s3))))
p.p1p4 <- cycle(list(list(c(10,3+s4),c(11,2+s4),c(12,1+s4))))

p.p5p2 <- cycle(list(list(c(1+s5,12+s2),c(2+s5,14+s2),c(3+s5,13+s2))))
p.p5p3 <- cycle(list(list(c(4+s5,13+s3),c(5+s5,12+s3),c(6+s5,14+s3))))
p.p5p4 <- cycle(list(list(c(7+s5,14+s4),c(8+s5,13+s4),c(9+s5,12+s4))))

p.p2p3 <- cycle(list(list(c(8+s2,18+s3),c(9+s2,17+s3),c(10+s2,16+s3),c(11+s2,15+s3))))
p.p3p4 <- cycle(list(list(c(8+s3,11+s4),c(9+s3,10+s4),c(10+s3,9+s4),c(11+s3,8+s4))))
p.p4p2 <- cycle(list(list(c(4+s4,18+s2),c(5+s4,17+s2),c(6+s4,16+s2),c(7+s4,15+s2))))

p.p6p1 <- cycle(list(list(c(1+s6,1),c(2+s6,3),c(3+s6,2))))
p.p6p2 <- cycle(list(list(c(4+s6,7+s2),c(5+s6,6+s2),c(6+s6,5+s2),c(7+s6,4+s2))))
p.p6p3 <- cycle(list(list(c(8+s6,7+s3),c(9+s6,6+s3),c(10+s6,5+s3),c(11+s6,4+s3))))
p.p6p4 <- cycle(list(list(c(15+s6,15+s4),c(16+s6,18+s4),c(17+s6,17+s4),c(18+s6,16+s4))))
p.p6p5 <- cycle(list(list(c(12+s6,10+s5),c(13+s6,12+s5),c(14+s6,11+s5))))

pall.pp <- p.p1p2 * p.p1p3 * p.p1p4 * p.p5p2 * p.p5p3 * p.p5p4 * p.p2p3 * p.p3p4 * p.p4p2 * p.p6p1 * p.p6p2 * p.p6p3 * p.p6p4 * p.p6p5

# 3角メッシュの置換
pall.tri <- pall.face * pall.edge^(-1)

# 同一のエッジから出ているペアがp.edge
# 同一のエッジに入ってくるペアがp.opposite
pall.opposite <- pall.tri * pall.edge * pall.tri^(-1)

# 各多面体でのzig
pall.zig <- pall.face * pall.tri^2
# 各多面体でのzagはpall.edge
pall.zag <- pall.edge
# zig.zagの２歩
pall.zigzag <- pall.zig * pall.zag

# 多面体を渡り歩くzig.zag
# １歩ごとに向き合った面の対応辺に渡る。ただし向きを逆にする

pall.zig.pp <- pall.zig * pall.pp

# zagして向き合い面の対応辺に渡る
pall.zag.pp <- pall.opposite * pall.pp

# 面を渡り歩いてのzig.zag二歩

pall.zigzag.pp <- pall.zig.pp * pall.zag.pp

# ２歩ずつ記録

n.iter <- 20

path.hx <- matrix(1:78,nrow=1)

for(i in 1:n.iter){
	current <- path.hx[length(path.hx[,1]),]
	tmp1 <- as.matrix(pall.zig.pp)[,current]
	tmp2 <- as.matrix(pall.zigzag.pp)[,current]
	path.hx <- rbind(path.hx,tmp1,tmp2)
}

# 3D空間に埋め込んでみる
X <- matrix(c(1,0,3,cos(2*pi/3),sin(2*pi/3),3,cos(4*pi/3),sin(4*pi/3),3,1,0,0,cos(2*pi/3),sin(2*pi/3),0,cos(4*pi/3),sin(4*pi/3),0,0,0,2,0,0,1),byrow=TRUE,ncol=3)

el <- matrix(c(1,2,2,3,3,1,4,5,5,6,6,4,1,4,2,5,3,6,7,8,1,7,2,7,3,7,4,8,5,8,6,8),byrow=TRUE,ncol=2)
plot3d(X,box = F, axes = F)
segments3d(X[c(t(el)),])

# 箙有向辺の始点は２頂点の中点。その２頂点IDの情報
arrow.st1 <- matrix(c(1,2,1,3,2,3,1,2,2,7,7,1,2,3,3,7,2,7,1,3,1,7,3,7),byrow=TRUE,ncol=2)
arrow.st2 <- matrix(c(1,2,1,7,2,7,1,2,2,5,5,4,4,1,2,7,7,8,8,5,5,2,8,5,8,4,4,5,1,7,1,4,4,8,8,7),byrow=TRUE,ncol=2)
arrow.st3 <- matrix(c(2,3,2,7,7,3,2,3,3,6,6,5,5,2,3,7,7,8,8,6,6,3,5,6,6,8,8,5,2,7,2,5,5,8,8,7),byrow=TRUE,ncol=2)
arrow.st4 <- matrix(c(1,3,3,7,7,1,1,7,7,8,8,4,4,1,7,3,3,6,6,8,8,7,8,6,6,4,4,8,3,1,1,4,4,6,6,3),byrow=TRUE,ncol=2)
arrow.st5 <- matrix(c(4,8,8,5,5,4,5,8,8,6,6,5,6,8,8,4,4,6,6,4,4,5,5,6),byrow=TRUE,ncol=2)
arrow.st6 <- matrix(c(3,1,1,2,2,3,2,1,1,4,4,5,5,2,3,2,2,5,5,6,6,3,5,4,4,6,6,5,4,1,1,3,3,6,6,4),byrow=TRUE,ncol=2)
arrow.st.all <- rbind(arrow.st1,arrow.st2,arrow.st3,arrow.st4,arrow.st5,arrow.st6)

arrow.st.X <- (X[arrow.st.all[,1],] + X[arrow.st.all[,2],])/2

arrow.stend.id <- cbind(1:96,c(as.matrix(pall.face)))

segments3d(arrow.st.X[c(t(arrow.stend.id)),],col=2)

plot3d(X,box = F, axes = F)
segments3d(X[c(t(el)),])

segments3d(arrow.st.X[c(arrow.stend.id[path.hx[1,1],]),],col=2)
segments3d(arrow.st.X[c(arrow.stend.id[path.hx[2,1],]),],col=3)
segments3d(arrow.st.X[c(arrow.stend.id[path.hx[3,1],]),],col=4)
segments3d(arrow.st.X[c(arrow.stend.id[path.hx[4,1],]),],col=5)
segments3d(arrow.st.X[c(arrow.stend.id[path.hx[5,1],]),],col=6)
segments3d(arrow.st.X[c(arrow.stend.id[path.hx[6,1],]),],col=7)
segments3d(arrow.st.X[c(arrow.stend.id[path.hx[7,1],]),],col=8)

# zigzagサイクルの長さ
len.zigzag.cycle <- apply(path.hx,2,function(x){which(x==x[1])[2]-1})

tmp.list <- list()
for(i in 1:length(len.zigzag.cycle)){
	tmp.list[[i]] <- path.hx[1:len.zigzag.cycle[i],i]
	tmp.min <- which(tmp.list[[i]]==min(tmp.list[[i]]))
	tmp.len <- length(tmp.list[[i]])
	if(tmp.min==1){
		tmp.list[[i]] <- tmp.list[[i]][1:tmp.len]
	}else if(tmp.min<tmp.len){
		tmp.list[[i]] <- tmp.list[[i]][c(tmp.min:tmp.len,1:(tmp.min-1))]
	}else{
		tmp.list[[i]] <- tmp.list[[i]][c(tmp.min,1:(tmp.min-1))]
	}

}
init.val <- sapply(tmp.list,function(x){x[1]})
uniques <- which(!duplicated(init.val))

zigzag.pp.cycles <- list()
cnt <- 1
for(i in 1:length(uniques)){
	zigzag.pp.cycles[[cnt]] <- tmp.list[[uniques[i]]]
	cnt <- cnt + 1
}
p.zigzag.pp <- cycle(list(zigzag.pp.cycles))
