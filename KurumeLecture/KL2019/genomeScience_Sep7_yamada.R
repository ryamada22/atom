n <- 100000
X <- sample(c(1,0),n,replace=TRUE)
Y <- sample(c(1,0),n,replace=TRUE)


tab <- table(X,Y)
chisq.test(tab,correct=FALSE)

fisher.test(tab)

# nを小さくする

n <- 30
X <- sample(c(1,0),n,replace=TRUE)
Y <- sample(c(1,0),n,replace=TRUE)


tab <- table(X,Y)
chisq.test(tab,correct=FALSE)

fisher.test(tab)

###
# 何度もやる

n.iter <- 1000 # 繰り返し回数
p.chisq <- p.fisher <- rep(0,n.iter)
for(i in 1:n.iter){
	n <- 30
	X <- sample(c(1,0),n,replace=TRUE)
	Y <- sample(c(1,0),n,replace=TRUE)


	tab <- table(X,Y)
	p.chisq[i] <- chisq.test(tab,correct=FALSE)$p.value

	p.fisher[i] <- fisher.test(tab)$p.value

}

hist(p.chisq)
hist(p.fisher)

plot(sort(p.chisq))
plot(sort(p.fisher))

