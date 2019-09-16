# 三角メッシュ・双対グラフとは無関係に
# 置換行列Wcと、それを辺ペア情報によって変換した置換行列Wdとを
# 作成する
# 三角メッシュの双対グラフの辺接続行列 We = Wc + Wdが満足するべき
# 条件の一部を満足するWc,Wdペアのランダム発生

d <- 8 　　　#dは偶数
Id <- diag(d)
s <- sample(d)
s
Wc <- Id[,s]
sp <- sample(d)
sp
sp.pair <- matrix(sp,ncol=2)
sp.pair
P <- matrix(0,d,d)
for(i in 1:length(sp.pair[,1])){
  P[sp.pair[i,1],sp.pair[i,2]] <- 1
  P[sp.pair[i,2],sp.pair[i,1]] <- 1
}
Wd <- P %*% t(Wc) %*% P


#ここから自作です．
# Wc,Wdは３正則グラフの辺を正逆２方向の２本の辺に分け
# その有向辺が接続する２つの有向辺の選択を意味する
# 今、３正則グラフを平面グラフと見れば
# Wc,Wdの選択は、辺をたどるパスに関して
# 左折・右折に相当する
# 以下では、Wc -> Wd -> Wc -> Wd (左折・右折を繰り返す)によって
# 移り変わる、辺idの系列を発生させている
# すべての有向辺を起点として、辺id列を発生させている
# パスの長さは有向辺の本数dに対して、 2d+1 としている
# その心は、(W = Wc * Wd)という左折・右折の２歩に相当するWが
# 置換行列であり、Wは１個以上の、相互にdisjointな順列に分解されることが
# わかっており、有向辺の本数(W,Wc,Wdの行数・列数)がその順列の周期の
# 最大数であることがわかっているので
# 巡回パスの検出には 2d + 1であり、十分である
# ちなみに、2d + 1は、１個の初期状態と、(左折・右折)の組をd回繰り返したときの
# 位置の推移の記録に相当する

# 1,2,...,dのそれぞれの有向辺を起点とする長さdのサイクルを納めるリスト
series <- list()

for(i in 1:d){ # d本の有向辺のそれぞれについて処理
  # i-th 有向辺が起点であることを示す長さdのベクトル v を作る
  v <- rep(0,d)
  v[i] <- 1
  # 左折・右折の２歩を d 回繰り返して、それを記録するためのベクトルとして
  # series[[i]] を作成する
  series[[i]] <- rep(0,2*d+1)　#同じ番号でも偶数番目と奇数番目に現れる状態を区別します．
  series[[i]][1] <- i
  # v.now は「歩いているエージェント」が立っている場所をベクトルとして表現したもの  
  v.now <- v # 初期状態
  for(j in 1:d){ # (左折・右折)の組のd回繰り返し
    v.now <- Wc %*% v.now # 左折
    series[[i]][2*j] <- which(v.now == 1)
    v.now <- Wd %*% v.now # 右折
    series[[i]][2*j+1] <- which(v.now == 1)
  }
}
series


cycles <- list()
cnt <- 1 # 新規サイクルを記録するためのカウンタ
found <- rep(0,2*d)
#2d個の部屋はそれぞれ，1(奇数番目),2(奇),...,d(奇),1(偶数番目),...,d(偶)に割り当てます．
for(i in 1:d){
  if(found[i] == 0){# すでに奇数歩目として通過している場合は処理しない
    tmp <- series[[i]]
    #tmpから奇数番目の成分を抽出したベクトルを作り，その番地を取得後，もとのtmpの番地に戻します．
    # odd.sames[1] は起点に相当
    # odd.sames[2] は初めて(左折・右折)後に自身に戻ってきたことに相当
    # ただし、歩数としては、odd.sames[2] * 2 -1歩目に戻ってくる
    # したがって、巡回は1,2,..,(odd.sames[2]*2-1-1)歩までの記録
    odd.sames <- which(tmp[seq(1,2*d+1,2)] == i)
    new.cycle <- tmp[1:(odd.sames[2]*2-2)]
    # 新規サイクルを記録してカウンタを更新
    cycles[[cnt]] <- new.cycle
    cnt <- cnt + 1
    #奇数番目，偶数番目それぞれ別個に記録します．
    # すべての有向辺が、奇数歩目にも偶数歩目にも通過したかどうかのフラグ
    # 長さ 2d のベクトル
    # 前半の長さdは奇数歩目についての情報
    # 後半の長さdは偶数歩目についての情報
    found[new.cycle[seq(1,length(new.cycle),2)]] <- 1
    found[new.cycle[seq(2,length(new.cycle),2)] + d] <- 1
  }
}
cycles
found

table(unlist(cycles))　#全ての辺を必ず2回ずつ通過している！？

# 処理を関数化すると、使いやすくなります。
# 上記のコードは３処理に相当するので、それぞれを関数にしてみます
# 処理を書いたら、基本的には、その処理を
####
## fx.name <- function(args){
##   
## }
# の {  } の間にコピーペーストすると、ほぼ完成しています
# 追加でやるべきことは、(args)のところに
# 処理のための「入力オブジェクト」を入れ
# { } の末尾に、return(hoge) と書いて
# 処理によって欲しかった結果を出力とすることでできます

# 第一処理。Wc,Wdのランダム作成を関数化してみます

my.randWcWd <- function(d){
	# d <- 8 　　　#dは引数になるのでコメントアウトします
	Id <- diag(d)
	s <- sample(d)
	s
	Wc <- Id[,s]
	sp <- sample(d)
	sp
	sp.pair <- matrix(sp,ncol=2)
	sp.pair
	P <- matrix(0,d,d)
	for(i in 1:length(sp.pair[,1])){
	  P[sp.pair[i,1],sp.pair[i,2]] <- 1
	  P[sp.pair[i,2],sp.pair[i,1]] <- 1
	}
	Wd <- P %*% t(Wc) %*% P
	
	# 出力します
	return(list(Wc=Wc,Wd=Wd))
}

# 関数を作ったら検算
d <- 100
WcWd <- my.randWcWd(d)
WcWd

# 同様に、Wc,Wdを与えて、seriesを返す関数を作ります
# 返り値オブジェクトがseriesであって、すでにそれをオリジナル処理コードで
# 書き出していたので、出力は修正不要です
# 引数(args)にWc,Wdを与えます
# べた書きコードでは上流で作成したオブジェクトを参照している場合があり
# それについて注意が必要です
# この場合は、d
my.seriesFromWcWd <- function(Wc,Wd){
	# 1,2,...,dのそれぞれの有向辺を起点とする長さdのサイクルを納めるリスト
	series <- list()
	# dが未設定なので書き足す
	d <- length(Wc[,1])
	
	for(i in 1:d){ # d本の有向辺のそれぞれについて処理
	  # i-th 有向辺が起点であることを示す長さdのベクトル v を作る
	  v <- rep(0,d)
	  v[i] <- 1
	  # 左折・右折の２歩を d 回繰り返して、それを記録するためのベクトルとして
	  # series[[i]] を作成する
	  series[[i]] <- rep(0,2*d+1)　#同じ番号でも偶数番目と奇数番目に現れる状態を区別します．
	  series[[i]][1] <- i
	  # v.now は「歩いているエージェント」が立っている場所をベクトルとして表現したもの  
	  v.now <- v # 初期状態
	  for(j in 1:d){ # (左折・右折)の組のd回繰り返し
	    v.now <- Wc %*% v.now # 左折
	    series[[i]][2*j] <- which(v.now == 1)
	    v.now <- Wd %*% v.now # 右折
	    series[[i]][2*j+1] <- which(v.now == 1)
	  }
	}
	series

}

# 検算
Wc <- WcWd$Wc
Wd <- WcWd$Wd
series.2 <- my.seriesFromWcWd(Wc,Wd)
series.2

# cyclesを返す処理も関数化しておきます

my.cyclesFromSeries <- function(series){
	cycles <- list()
	cnt <- 1 # 新規サイクルを記録するためのカウンタ
	found <- rep(0,2*d)
	#2d個の部屋はそれぞれ，1(奇数番目),2(奇),...,d(奇),1(偶数番目),...,d(偶)に割り当てます．
	d <- length(series)
	for(i in 1:d){
	  if(found[i] == 0){
	    tmp <- series[[i]]
	    #tmpから奇数番目の成分を抽出したベクトルを作り，その番地を取得後，もとのtmpの番地に戻します．
	    # odd.sames[1] は起点に相当
	    # odd.sames[2] は初めて(左折・右折)後に自身に戻ってきたことに相当
	    # ただし、歩数としては、odd.sames[2] * 2 -1歩目に戻ってくる
	    # したがって、巡回は1,2,..,(odd.sames[2]*2-1-1)歩までの記録
	    odd.sames <- which(tmp[seq(1,2*d+1,2)] == i)
	    new.cycle <- tmp[1:(odd.sames[2]*2-2)]
	    # 新規サイクルを記録してカウンタを更新
	    cycles[[cnt]] <- new.cycle
	    cnt <- cnt + 1
	    #奇数番目，偶数番目それぞれ別個に記録します．
	    # すべての有向辺が、奇数歩目にも偶数歩目にも通過したかどうかのフラグ
	    # 長さ 2d のベクトル
	    # 前半の長さdは奇数歩目についての情報
	    # 後半の長さdは偶数歩目についての情報
	    found[new.cycle[seq(1,length(new.cycle),2)]] <- 1
	    found[new.cycle[seq(2,length(new.cycle),2)] + d] <- 1
	  }
	}
	# 知りたい情報がcycles,found,table(unlist(cycles))と複数なので、それらをリストで返します
	# cycles
	# found

	# table(unlist(cycles))　#全ての辺を必ず2回ずつ通過している！？
	
	return(list(cycles=cycles,found=found,tab =table(unlist(cycles))))

}

# 検算

cycles.out <- my.cyclesFromSeries(series.2)

cycles.out