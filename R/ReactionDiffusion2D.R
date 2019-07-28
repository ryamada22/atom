# 格子グラフを作ろう

# ２次元座標を作る
x <- 1:50
y <- 1:50
xy <- expand.grid(x,y)
plot(xy)

# 頂点間距離を測る。ただし、マンハッタン距離
d <- as.matrix(dist(xy,method="manhattan"))

# マンハッタン距離=1のときにエッジを引く
d.01 <- d==1

library(igraph)

g <- graph.adjacency(d.01,mode = "undirected")

# レイアウト法としてgrid法を適用
plot(g,vertex.label="",vertex.size=5,layout=layout_on_grid(g,dim=2))

# 格子グラフの隣接行列は d.01

## X, Y
## 拡散係数
Dx <- 2.8
Dy <- 5

## 反応
# dx/dt = P x + Q x^3 + R*y
# dy/dt = S x  + U y

P <- 1
Q <- -1
R <- -1
S <- 1
U <- 1
# 小部屋数
n.cell <- length(xy[,1])
# dt
dt <- 0.01
# 繰り返し数
n.step <- 100

# 時刻別・小部屋別の量を記録する行列
X <- Y <- matrix(0,n.step,n.cell)

# 初期値
X[1,] <- rep(1,n.cell) + rnorm(n.cell,0,0.1)
Y[1,] <- rep(1,n.cell) + rnorm(n.cell,0,0.1)

# 拡散INOUT行列
In <- d.01
Out <- diag(apply(In,1,sum))

InOut <- In - Out

# 単位行列
I <- diag(n.cell)


for(i in 2:n.step){
	# まず拡散
	nowX <- X[i-1,]
	nowY <- Y[i-1,]
	X[i,] <- (I + dt * InOut * Dx) %*% nowX
	Y[i,] <- (I + dt * InOut * Dy) %*% nowY
	
	# ついで反応
	
	dx <- dt * (P*X[i,] + Q*X[i,]^3 + R*Y[i,])
	dy <- dt * (S*X[i,]*Y[i,] + U*Y[i,])
	
	X[i,] <- X[i,] + dx
	Y[i,] <- Y[i,] + dy
}

par(mfcol=c(1,2))
matplot(X,type="l")
matplot(Y,type="l")
par(mfcol=c(1,2))

for(i in 1:n.step){
	#if(i%%10==0){
		image(matrix(X[i,],length(x),length(y)),main="X")
	
		image(matrix(Y[i,],length(x),length(y)),main="Y")
	#}
	for(i in 1:100){
		eigen(matrix(rnorm(10^4),10^2,10^2))
	}
}

