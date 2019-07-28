# 空間を小部屋に分けて拡散を離散的に扱う

# 「隣」の小部屋には砂粒が移動できて
# 「非隣」の小部屋には砂粒が移動できない

# 「隣」「非隣」をグラフで表してみる

install.packages("igraph")
library(igraph)

# ５部屋が横一列に並んでいる場合
edge.list <- rbind(c(1,2),c(2,3),c(3,4),c(4,5))

g <- graph.edgelist(edge.list,directed=FALSE)

# できた
plot(g)

# グラフオブジェクト g には、「隣・非隣情報」を行列(隣接行列)が付随している

g.mat <- get.adjacency(g)

g.mat

# 見にくいので、ちょっと工夫する
g.mat <- as.matrix(g.mat)
g.mat

# ５部屋の初期値を 1,2,4,8,16とする
x.init <- c(1,2,4,8,16)

# 以下の掛け算をしてみる
g.mat %*% x.init

# 現れた数字は隣の部屋の値の和になっている

# ある部屋の砂粒の10%が右隣りに、同じく10%が左隣に移動するとすれば
# g.matの成分値 1 を 0.1 にすればよさそうだ

p <- 0.1

g.mat.p <- p * g.mat
g.mat.p

# 移動量を再計算してみる

g.mat.p %*% x.init

# 確かに移動量が出た

# 移動して出て行った分はどうするか

# 各行の和を取れば、それが「(左右はいざ知らず)出て行った割合」がわかる
out.frac <- apply(g.mat.p,1,sum)
out.frac

# 出て行った分については、自身に関する変化だから、行列でいえば、対角成分

self.out <- diag(out.frac)
self.out

# 増減を併せれば

inANDout <- g.mat.p - self.out
inANDout

# これを初期値に掛けてみる

inANDout %*% x.init

# これは増減

# 増減する前の量は、行列を使えば、単位行列を掛けること
I <- diag(5)
I

I %*% x.init

# 微小時間後の小部屋の砂の量は、増減前の値に増減値を加えたもの
# I %*% x.init + inANDout %*% x.init = (I + inANDout) %*% x

(I + inANDout) %*% x.init

#############
# 準備OK!


# これを微小変化について繰り返せば、「拡散」の離散版

p <- 0.01

g.mat.p <- p * g.mat
g.mat.p

# 移動量を再計算してみる

g.mat.p %*% x.init

# 確かに移動量が出た

# 移動して出て行った分はどうするか

# 各行の和を取れば、それが「(左右はいざ知らず)出て行った割合」がわかる
out.frac <- apply(g.mat.p,1,sum)
out.frac

# 出て行った分については、自身に関する変化だから、行列でいえば、対角成分

self.out <- diag(out.frac)
self.out

# 増減を併せれば

inANDout <- g.mat.p - self.out
inANDout

# これを初期値に掛けてみる

inANDout %*% x.init

# これは増減

# 増減する前の量は、行列を使えば、単位行列を掛けること
I <- diag(5)
I

I %*% x.init


#
# 微小時間後の小部屋の砂の量は、増減前の値に増減値を加えたもの
# I %*% x.init + inANDout %*% x.init = (I + inANDout) %*% x

(I + inANDout) %*% x.init

# 繰り返し回数
n.step <- 2000

X <- matrix(0,n.step,5)
X[1,] <- x.init

for(i in 2:n.step){
	nowX <- X[i-1,]
	X[i,] <- (I + inANDout) %*% nowX
}

matplot(X,type="l")

###################
# 空間をグラフにできればOK
# 円周状に閉じたグラフはどう作る？

# edge.list <- rbind(c(1,2),c(2,3),c(3,4),c(4,5))
edge.list2 <- rbind(c(1,2),c(2,3),c(3,4),c(4,5),c(5,1))

g2 <- graph.edgelist(edge.list2,directed=FALSE)

# できた
plot(g2)

# 隣接行列を作る
g.mat2 <- get.adjacency(g2)

g.mat2

g.mat2 <- as.matrix(g.mat2)
g.mat2

# 微小時間化する
p <- 0.01

g.mat.p2 <- p * g.mat2
g.mat.p2

# 出て行く分
out.frac2 <- apply(g.mat.p2,1,sum)
out.frac2

self.out2 <- diag(out.frac2)
self.out2

# 増減を併せる

inANDout2 <- g.mat.p2 - self.out2
inANDout2

# 繰り返し回数
#n.step <- 4000

X2 <- matrix(0,n.step,5)
X2[1,] <- x.init

for(i in 2:n.step){
	nowX <- X2[i-1,]
	X2[i,] <- (I + inANDout2) %*% nowX
}

par(mfcol=c(1,2))
matplot(X,type="l")
matplot(X2,type="l")
par(mfcol=c(1,1))

############
# 小部屋を増やそう
g.mat
# じっと見る

n.cell <- 20

g.mat.n <- diag(n.cell)
par(mfcol=c(2,2))
image(g.mat.n)

g.mat.n <- g.mat.n[,c(2:n.cell,1)]
image(g.mat.n)
g.mat.n[1,n.cell] <- 0
image(g.mat.n)

g.mat.n <- g.mat.n + t(g.mat.n)
image(g.mat.n)
par(mfcol=c(1,1))

# もっと増やす
n.cell <- 100

g.mat.n <- diag(n.cell)
par(mfcol=c(2,2))
image(g.mat.n)

g.mat.n <- g.mat.n[,c(2:n.cell,1)]
image(g.mat.n)
g.mat.n[1,n.cell] <- 0
image(g.mat.n)

g.mat.n <- g.mat.n + t(g.mat.n)
image(g.mat.n)
par(mfcol=c(1,1))


# さて、実験
# 微小時間化する
p <- 0.01

g.mat.pn <- p * g.mat.n
g.mat.pn

# 出て行く分
out.fracn <- apply(g.mat.pn,1,sum)
out.fracn

self.outn <- diag(out.fracn)
self.outn

# 増減を併せる

inANDoutn <- g.mat.pn - self.outn
inANDoutn

# 繰り返し回数
n.step <- 100000

Xn <- matrix(0,n.step,n.cell) # 5 -> n.cell

x.init <- sin((1:n.cell)/4) + 1

plot(x.init,type="l")

Xn[1,] <- x.init

In <- diag(n.cell)

for(i in 2:n.step){
	nowX <- Xn[i-1,]
	Xn[i,] <- (In + inANDoutn) %*% nowX
}

matplot(Xn,type="l")

# 異なる初期値


Xn <- matrix(0,n.step,n.cell) # 5 -> n.cell

x.init <- runif(n.cell)^3

plot(x.init,type="l")

Xn[1,] <- x.init

In <- diag(n.cell)

for(i in 2:n.step){
	nowX <- Xn[i-1,]
	Xn[i,] <- (In + inANDoutn) %*% nowX
}

matplot(Xn,type="l")

