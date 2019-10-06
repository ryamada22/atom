library(onion)
my.runit.imq <- function(n=1){
	tmp <- matrix(rnorm(n*3),ncol=3)
	tmp <- Hi * tmp[,1] + Hj * tmp[,2] + Hk * tmp[,3]
	tmp <- tmp/sqrt(Norm(tmp))
	return(tmp)
}
q1 <- my.runit.imq()
q2 <- my.runit.imq()
q3 <- -(q1+q2)
q4 <- my.runit.imq()
q5 <- -q4+q3

tmp <- matrix(rnorm(n*3),ncol=3)
tmp <- tmp/sqrt(apply(tmp^2,1,sum))
v1 <- tmp[1,]
v2 <- tmp[2,]
v3 <- -(v1+v2)
v4 <- tmp[3,]
v5 <- -v4 + v3

v1 + v2 + v3
-v3 + v4 + v5

apply(tmp^2,1,sum)
apply(tmp,2,sum)

u1 <- exp(q1)
u2 <- exp(q2)
u3 <- exp(q3)
u4 <- exp(q4)
u5 <- exp(q5)

u1 * u2 * u3
1/u3 * u4 * u5

u1. <- u1/(u3+1)*u3
u2. <- u2 * (u3+1)
u3. <- 1/u3
u4. <- u4/(u3+1)*u3
u5. <- u5 * (u3+1)

u1. * 1/u3. * u5.
u2. * u4. * u3.
