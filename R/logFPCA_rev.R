my.logFPCA <- function(P,MIN.val=10^(-10)){
	P[which(P<=0)] <- MIN.val # caring zero prob.density
	Q <- log(P)
	H <- Q %*% t(Q)
	eigen.out <- eigen(H)
	e.val <- eigen.out[[1]]
	V <- eigen.out[[2]]
	# only the components with > 0 eigen values are required
	n.positives <- length(which(e.val > 0))
	e.val. <- e.val[1:n.positives]
	# only the positive eigen vectors are required for theta
	Theta <- V[,1:n.positives] %*% diag(sqrt(e.val.))

	F <- MASS::ginv(Theta) %*% Q
	return(list(Theta = Theta,F=F,e.value= e.val,V=V))
}

my.logLike.Ftheta <- function(sample,Fx,theta){
	psi <- log(sum(exp(theta %*% Fx)))
	
	LL <- sum((theta %*% Fx - psi) [sample])
	return(LL)
}
s <- sample(1:length(P[1,]),50,replace=TRUE,prob=P[1,])

thetas <- rnorm(length(out$F[,1]))

theta1 <- seq(from=-1000,to=100,length=500)

LLs <- rep(0,length(theta1))

for(i in 1:length(theta1)){
	LLs[i] <- my.logLike.Ftheta(s,out$F,c(theta1[i],thetas[2:length(thetas)]))
}

plot(theta1,LLs)

library(deef)


P <- Distset2D$P

# Let's make discrete samples
samples <- matrix(0,length(P[,1]),length(P[1,]))
N <- 30 * length(P[1,])

# zero-inflation
for(i in 1:length(P[,1])){
	for(j in 1:length(P[1,])){
		if(runif(1) < 0.2){
			samples[i,j] <- 0
		}else{
			samples[i,j] <- rpois(1,N*P[i,j])
		}
		
	}
}

P <- samples/apply(samples,1,sum)

P. <- P
P.[which(P.<= 0)] <- 10^(-10)
Q <- log(P.)

out <- my.logFPCA(P)

# original discrete values are nicely recovered
range(exp(out$Theta %*% out$F) - P)

# Smoothen F with k-window averaging

out.F <- out$F
k <- 40
F.sm <- matrix(0,length(out.F[,1]),length(out.F[1,])-k)
exp.F.sm <- F.sm
for(i in 1:(length(out.F[1,])-k)){
	F.sm[,i] <- apply(out.F[,i:(i+k)],1,mean)
	exp.F.sm[,i] <- apply(exp(out.F[,i:(i+k)]),1,mean)
}
for(i in 1:(k/2)){
	F.sm <- cbind(F.sm[,1],F.sm,F.sm[,length(F.sm[1,])])
	exp.F.sm <- cbind(exp.F.sm[,1],exp.F.sm,exp.F.sm[,length(exp.F.sm[1,])])
}

F.sm.per.exp <- log(exp.F.sm)

dim(F.sm)
dim(out.F)
P.sm.2 <- exp(out$Theta[,1:4] %*% F.sm[1:4,])
P.sm.2 <- P.sm.2/apply(P.sm.2,1,sum)
P.sm.3 <- exp(out$Theta %*% (F.sm.per.exp- apply(F.sm.per.exp,1,mean)))
P.sm.3 <- P.sm.3/apply(P.sm.3,1,sum)


Theta.sm <- Q %*% MASS::ginv(F.sm)

Q.sm <- Theta.sm %*% F.sm

P.sm <- exp(Q.sm)
P.sm <- P.sm/apply(P.sm,1,sum)

for(i in 1:length(P[,1])){
	matplot(cbind(P.sm[i,],Distset2D$P[i,]),type="l")
}

est.testings <- my.est.density(samples,F.sm)

k <- 30

matplot(cbind(est.testings[k,],Distset2D$P[k,]),type="l")

my.eta <- function(new.sample,F){
	1/(sum(new.sample)) * F %*% new.sample
}

my.theta <- function(new.sample,F,MIN.val=10^(-10)){
	#new.sample <- new.sample/sum(new.sample)
	new.sample <- new.sample/apply(new.sample,1,sum)
	new.sample[which(new.sample<=0)] <- MIN.val
	Q <- log(new.sample)
	F %*% t(Q)
}

my.est.density <- function(new.sample,F){
	theta <- my.theta(new.sample,F)
	if(! is.matrix(theta)){
		theta <- matrix(theta,nrow=1)
	}
	tmp <- t(theta) %*% F
	tmp - mean(tmp)
	est.new.sample <- exp(tmp)
	est.new.sample/apply(est.new.sample,1,sum)
}

trainings <- sample(1:length(P[,1]), 100)
testings <- (1:length(P[,1]))[-trainings]


out.trainings <- my.logFPCA(P[trainings,])

samples <- matrix(0,length(testings),length(P[1,]))
N <- 1000

for(i in 1:length(testings)){
	for(j in 1:length(new.sample)){
		samples[i,j] <- rpois(1,N*P[testings[i],j])
	}
}

est.testings <- my.est.density(samples,out$F)

for(i in 1:length(testings)){
	matplot(cbind(est.testings[i,],P[testings[i],]),type="l")
}

samples.trainings <- matrix(0,length(trainings),length(P[1,]))
N <- 1000

for(i in 1:length(trainings)){
	for(j in 1:length(new.sample)){
		samples.trainings[i,j] <- rpois(1,N*P[trainings[i],j])
	}
}

est.trainings <- my.est.density(samples.trainings,out$F)

for(i in 1:length(trainings)){
	matplot(cbind(est.trainings[i,],P[trainings[i],]),type="l")
}


#matplot(t(samples),type="l")


new.sample <- rep(0,length(P[target,]))
for(i in 1:length(new.sample)){
	new.sample[i] <- rpois(1,N*P[target,i])
}
plot(new.sample)

eta <- my.eta(new.sample,out$F)
theta <- my.theta(new.sample,out$F)



est.new.sample <- my.est.density(new.sample,out$F)

matplot(cbind(est.new.sample,P[target,]),type="l")

