n <- 20
S1 <- S2 <- diag(rep(1, n))
S1 <- S1[sample(n), ]
S2 <- S2[sample(n), ]

N <- 30

T1 <- t(S2) %*% t(S1) %*% S2
T1 <- S2
v <- matrix(1:n, ncol = 1)
v.hx <- matrix(0, N + 1, n)
v.hx[1, ] <- v
cnt <- 2
for (i in 1:(N/2)) {
  v.hx[cnt, ] <- S1 %*% matrix(v.hx[cnt - 1, ], ncol = 1)
  v.hx[cnt + 1, ] <- T1 %*% matrix(v.hx[cnt, ], ncol = 1)
  cnt <- cnt + 2
}
matplot(v.hx, type = "l")

eoutS <- eigen(S1)
eoutT <- eigen(T1)

eoutTS <- eigen(T1 %*% S1)

fS <- fractions(Arg(eoutS[[1]])/(2*pi))
fT <- fractions(Arg(eoutT[[1]])/(2*pi))
fTS <- fractions(Arg(eoutTS[[1]])/(2*pi))

strsplit(as.character(fS),"/")
