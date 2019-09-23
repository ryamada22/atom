# ２つの対応する点座標セット行列の積である3x3行列のSVD分解から、二乗ノルム最小化回転行列を算出する

# 単純なデータを作る
library(GPArotation) # Random rotation matrixを作る
R <- Random.Start(3)

X <- matrix(rnorm(50*3),ncol=3)
X <- X/sqrt(apply(X^2,1,sum))
Y <- X + rnorm(length(X),0,0.001)
Y <- Y/sqrt(apply(Y^2,1,sum))
Y <- t(R %*% t(Y))

H <- t(X) %*% Y

#H <- t(X1.rep) %*% X2.rep
# SVD 分解
svd.out <- svd(H)
# 分解結果から、回転行列の推定
R. <- (svd.out$v) %*% t(svd.out$u)

# 正解行列Rと推定行列R.とを比較する
range(R - R.)

# 点集合座標行列を作る
n1 <- 50

tmp <- matrix(runif(n1*3),ncol=3)
tmp <- X1/sqrt(apply(X1^2,1,sum))
n.iter <- 3
X1 <- matrix(0,0,3)
for(i in 1:n.iter){
	tmp.R <- Random.Start(3)
	tmp2 <- tmp %*% tmp.R
	X1 <- rbind(X1,tmp2)
}

R <- Random.Start(3)


X2 <- t(R %*% t(X1)) + rnorm(length(X1),0,0.001)
X2 <- X2/sqrt(apply(X2^2,1,sum))

X2 <- X2[sample(1:length(X2[,1]),length(X2[,1])*0.9),]

n1 <- length(X1[,1])
n2 <- length(X2[,2])
# この点集合行列は、点の数も違えば、点の対応じょうほうもない
# すべての点ペアが「対応する」とみなして、点行列を大きくする

X1.rep <- X1[rep(1:n1,n2),]
X2.rep <- X2[rep(1:n2,each=n1),]

H <- t(X1.rep) %*% X2.rep
# SVD 分解
svd.out <- svd(H)
# 分解結果から、回転行列の推定
R. <- (svd.out$v) %*% t(svd.out$u)

# 正解行列Rと推定行列R.とを比較する
range(R - R.)

