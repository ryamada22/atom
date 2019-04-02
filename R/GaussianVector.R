# Complex Gaussian random variable
my.cpx.rnorm <- function(n){
	x <- rnorm(n,0,1/sqrt(2))
	y <- rnorm(n,0,1/sqrt(2))
	return(x + 1i*y)
}

n <- 1000000
rcpx <- my.cpx.rnorm(n)
# ���ς�0
mean(rcpx) 
# E(Z * Conj(Z)) ��1
tmp <- rcpx * Conj(rcpx)
mean(tmp)

# E(Z^m * Conj(Z)^n)�́Am=n�̂Ƃ��Am!
# �����łȂ����0

m <- n <- 0:5
emn <- matrix(0,length(m),length(n))
for(i in 1:length(m)){
	for(j in 1:length(n)){
		tmp <- rcpx^m[i] * Conj(rcpx)^n[j]
		emn[i,j] <- mean(tmp)
	}
}
diag(emn)
diag(emn)/factorial(m) # ���̒l��1�ɂȂ�

# m = n �̂Ƃ��́Afactorial(m)�̒l�ɂȂ�A����ȊO�͐�Βl��0�ɂȂ�
matplot(Mod(emn)/factorial(m),type="l")

# m = n �̂Ƃ��́A�����ɍL���������
plot(Z * Conj(Z))
plot(Z * Z * Conj(Z) * Conj(Z))
# m != n�̂Ƃ��́A���f���ʂŁA���_���S�̗��_�����
plot(Z)
plot(Z * Z * Conj(Z))


image(emn)



# Gaussian Unitary Ensemble
# Complex Gaussian random variable
my.cpx.rnorm <- function(n){
	x <- rnorm(n,0,1/sqrt(2))
	y <- rnorm(n,0,1/sqrt(2))
	return(x + 1i*y)
}


my.makeMatFromVec <- function(v,d){
	ret <- diag(v[1:d])
	tmp <- v[(d+1):length(v)]
	z <- tmp[1:(length(tmp)/2)] + 1i * tmp[(length(tmp)/2+1):length(tmp)]
	z. <- Conj(z)
	ret[upper.tri(ret)] <- z
	ret <- ret + t(Conj(ret))
	diag(ret) <- v[1:d]
	return(ret)
}
# n : number of matrices to be generated
# d : dimension of matrices
my.GUE <- function(n,d){
	xd <- matrix(rnorm(n*d,0,1/sqrt(d)),nrow=n)
	xnond <- my.cpx.rnorm(n*d*(d-1)/2)/sqrt(d)
	xynond <- cbind(matrix(Re(xnond),nrow=n),matrix(Im(xnond),nrow=n))
	return(cbind(xd,xynond))
}
n <- 10000
d <- 3
rmat <- my.GUE(n,d)

# dxd�s�񂪕Ԃ�
my.makeMatFromVec(rmat[1,],d)

# �e�s��dxd�s��̌̐������Ԃ�
outmat <- t(apply(rmat,1,my.makeMatFromVec,d=d))

# N^2�����̕��ς�0
apply(rmat,2,mean)
# N^2�����̕��U��1/2�܂���1/(2d)
apply(rmat,2,var)

# N^2�����̕��U�����U�s���
# �Ίp������(1/N,1/N,...,1/2N,1/2N,...)�Ŕ�Ίp������0
covmat <- cov(rmat)
persp(covmat)

# ���_�I�ȁA����<BX,X>���K�肷��s��B�����
B <- diag(c(rep(d,d),rep(2*d,d*(d-1)/2*2)))

# ����B��p���āA�e�����_��d^2�x�N�g���̓��ς��v�Z����
BX <- B %*% t(rmat)

BX.X <- apply(t(BX) * rmat,1,sum)


# �����A�����_���s��̂��ꂼ��̓��̃g���[�X���v�Z����
Tr.Xsq <- rep(0,n)

for(i in 1:n){
	tmp <- matrix(outmat[i,],d,d)
	Tr.Xsq[i] <- sum(diag(tmp%*%tmp))
}
# ���҂�d��������֌W�ɂ���
plot(BX.X,Tr.Xsq)

# �����_���s��̐�������m�������P�ɕW��������

my.standard.rmat <-function(X,d){
	X[,1:d] <- X[,1:d] * sqrt(d)
	X[,(d+1):length(X[1,])] <- X[,(d+1):length(X[1,])] * sqrt(2*d)
	return(X)
}

rmat.st <- my.standard.rmat(rmat,d)

apply(Mod(rmat.st)^2,2,mean)




# ���K�����x�N�g���𔭐�������
# ����l�����_���s��W�𔭐�������
library(rWishart)
d <- 8 # �ϐ��̐�
W <- rWishart(1,d,diag(1,d))[,,1]
eigen(W) # �m�F
#W <- diag(1,d)

rmat <- rmvnorm(n=500000,sigma=W)

apply(rmat,2,mean)
apply(rmat,2,var)

# �������̕ϐ���I��
k <- 6
vs <- sample(1:d,k)

r.selected <- rmat[,vs]
EXselected <- mean(apply(r.selected,1,prod))

pairs <- my.all.pairs(k)

EXselected2 <- 0
for(i in 1:length(pairs[,1])){
	tmp0 <- 1
	for(j in 1:(length(pairs[1,])/2)){
		tmp <- r.selected[,pairs[i,(2*j-1):(2*j)]]
		tmp.j <- mean(apply(tmp,1,prod))
		tmp0 <- tmp0 * tmp.j
	}
	#tmp1 <- r.selected[,pairs[i,1:2]]
	#tmp2 <- r.selected[,pairs[i,3:4]]
	#m1 <- mean(apply(tmp1,1,prod))
	#m2 <- mean(apply(tmp2,1,prod))
	EXselected2 <- EXselected2 + tmp0
}

EXselected
EXselected2


my.pair.text <- function(x){
	tmp.mat <- matrix(x,byrow=TRUE,ncol=2)
	tmp2 <- apply(tmp.mat,1,paste0,collapse="",sep=".")
	return(tmp2)
}
# k: ����
my.all.pairs <- function(k){
	n <- k/2
	ret <- list()
	for(i in 1:(n)){
		ret[[i]] <- combn(k-(i-1)*2,2)		
	}
	ret[[1]] <- ret[[1]]
	lenret <- sapply(ret,ncol)
	tmplist <- list()
	for(i in 1:length(lenret)){
		tmplist[[i]] <- 1:lenret[i]
	}
	cmb <- expand.grid(tmplist)
	ans <- matrix(0,length(cmb[,1]),k)
	cnt <- 1
	for(i in 1:length(cmb[,1])){
		tmp <- c()
		v <- 1:k
		for(j in 1:length(cmb[1,])){
			tmp <- c(tmp,v[ret[[j]][,cmb[i,j]]])
			v <- v[-ret[[j]][,cmb[i,j]]]
		}
		ans[cnt,] <- tmp
		cnt <- cnt + 1
	}
	# �Q�d�ɐ����Ă��邩��E�E�E
	anstext <- t(apply(ans,1,my.pair.text))
	anstext.sort <- t(apply(anstext,1,sort))
	anstext.u <- unique(anstext.sort)
	ans.unique <- matrix(as.numeric(unlist(apply(anstext.u,1,strsplit,"\\."))),ncol=k,byrow=TRUE)
	return(ans.unique)
}

my.all.pairs(6)