# 三角メッシュを貼り合わせる
# t1,t2は三角面の３列行列
# f1,f2はt1,t2の行を指定する
# s1,s2はf1[t1,(s1,s1+1,s1+2)] と f2[s2,(s2,s2-1,s2-2)]とを貼り合わせる
my.paste.tri <- function(t1,t2,f1,f2,s1=1,s2=1){
	# t1の頂点IDをt2の頂点IDより後の値にする
	t2 <- t2 + max(t1)
	# t2の張り合わせ三角形の頂点ID
	v2s <- t2[f2,]
	
	# その上で、max.v2-c(0,1,2) を f1[t1,]に適切に変換する
	
	if(s1==1){
		v1s <- t1[f1,c(1,2,3)]
	}else if(s1==2){
		v1s <- t1[f1,c(2,3,1)]
	}else{
		v1s <- t1[f1,c(3,1,2)]
	}
	if(s2==1){
		v2s <- t2[f2,c(1,3,2)]
	}else if(s2==2){
		v2s <- t2[f2,c(2,1,3)]
	}else{
		v2s <- t2[f2,c(3,2,1)]
	}
	t2.new <- t2
	t2.new[which(t2 == v2s[1])] <- v1s[1]
	t2.new[which(t2 == v2s[2])] <- v1s[2]
	t2.new[which(t2 == v2s[3])] <- v1s[3]
	
	ret <- rbind(t1[-f1,],t2.new[-f2,])
	u <- unique(c(ret))
	r <- rank(u)
	ret2 <- ret
	for(i in 1:length(u)){
		ret2[which(ret==u[i])] <- r[i]
	}
	return(ret2)
}
tetra <- rbind(c(1,2,3),c(1,3,4),c(1,4,2),c(3,2,4))

out <- my.paste.tri(tetra,tetra,1,3,s1 = 1,s2 = 2)
