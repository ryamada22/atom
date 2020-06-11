

library(rgl)

# 正四面体の４頂点座標を計算する関数
CategoryVector<-function (d = 3){
  df <- d - 1
  diagval <- 1:d
  diagval <- sqrt((d)/df) * sqrt((d - diagval)/(d - diagval + 1))
  others <- -diagval/(d - (1:d))
  m <- matrix(rep(others, d), nrow = d, byrow = TRUE)
  diag(m) <- diagval
  m[upper.tri(m)] <- 0
  as.matrix(m[, 1:df])
}

# 正四面体
Tetra <- CategoryVector(4)
plot3d(Tetra)
# 正四面体の４頂点座標
Va <- Tetra[1,]
Vb <- Tetra[2,]
Vc <- Tetra[3,] 
Vd <- Tetra[4,] 

theta <- 0.2
Va <- c(0,1/2,0)
Vb <- c(-sqrt(3)/2,0,0)
Vc <- c(0,-1/2,0)
Vd <- c(sqrt(3)/2 * cos(thea),0,sqrt(3)/2 * sin(theta))

# 正四面体の６辺ベクトルを適当な向きでとる
AB <- Vb - Va
AC <- Vc - Va
AD <- Vd - Va
BC <- Vc - Vb
CD <- Vd - Vc
BD <- Vd - Vb

sum(AB * CD) + sum(BC * (-AD)) - sum(AC * (-BD))



sum(AB * CD) + sum(BC * (-AD)) - sum(AC * (-BD))


# Cross積・ベクトル外積計算関数
# V x U = (v2u3-v3u2,v3u1-v1u3, v1u2-v2u1)

Here is a generalized cross product:
  
  xprod <- function(...) {
    args <- list(...)
    
    # Check for valid arguments
    
    if (length(args) == 0) {
      stop("No data supplied")
    }
    len <- unique(sapply(args, FUN=length))
    if (length(len) > 1) {
      stop("All vectors must be the same length")
    }
    if (len != length(args) + 1) {
      stop("Must supply N-1 vectors of length N")
    }
    
    # Compute generalized cross product by taking the determinant of sub-matricies
    
    m <- do.call(rbind, args)
    sapply(seq(len),
           FUN=function(i) {
             det(m[,-i,drop=FALSE]) * (-1)^(i+1)
           })
  }

# 箙変異に話を戻す
# 2つの三角形ABC と ACD とが辺ACを共有しているとする
# ACは四角形ABCDの対角線のようにも見える
# このACをDBに付け替えるのが箙の変異

# 箙の変異により、次のような変数変換が知られている
# X(PQ) を辺PQに付けた変数とする
# X(newAC = DB) -> (X(AB) * X(CD) + X(BC) * X(DA)) / X(AC)
# "->" の右辺の分母を、左辺に移行してみる（どんな代数なのかわからないが…）
# X(DB) * X(AC) -> (X(AB) * X(CD) + X(BC) * X(DA))
# "->" の両辺が、何らかの意味で等しいことがわかると、何かが見えてくると期待したい

# X(AB) * X(CD) を、辺ABベクトルと辺CDベクトルとのクロス積だとしてみよう
# そうすれば、"->"の左辺も右辺も３次元ベクトルである

XABxXCD <- xprod(AB,-CD) 
XBCxXDA <- xprod(BC,AD)



Uhen <- XABxXCD + XBCxXDA

Sahen.mat <- matrix(c(10^(-8),-BD[3],BD[2],BD[3],10^(-8),-BD[1],-BD[2],BD[1],10^(-8)),byrow=TRUE,3,3)

Sahen.mat.inv <- solve(Sahen.mat)
Sahen.mat.inv

Sahen.mat.inv %*% matrix(Uhen,ncol=1)



Sahen <- xprod(-BD,AC)

Uhen
Sahen

xprod((ABxCD + BCxDA),-AC)

BD
