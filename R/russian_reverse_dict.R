---
title: "逆引きロシア単語"
author: "Ryo Yamada"
date: "2021/12/21"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## キリル文字

全33文字、大文字と小文字。

```{r}
uppers <- "АБВГДЕЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯЁ"
lowers <- "абвгдежзийклмнопрстуфхцчшщъыьэюяё"
uppers <- unlist(strsplit(uppers,""))
lowers <- unlist(strsplit(lowers,""))
```

## 品詞

```{r}
hinshi.jp <- c("名詞","代名詞","動詞","形容詞")
hinshi <- c("N","P","V","A") # Noun, Pronoun,Verb,Addjective
```

## 文法特性

性、単・複、六格

```{r}
sex <- c("m","f","n")
num <- c("s","p")
kaku <- c("Syu","Sei","Yo","Tai","Zou","Zenchi") # 主・生・与・対・造・前置
```

## 品詞の基本特性

```{r}
N.params <- list(sex,num,kaku)
P.params <- list(sex,num,kaku)
V.params <- list(sex,num,kaku)
A.params <- list(sex,num,kaku)
```


## Yaruku辞書からの読み込み

冒頭には検索した文字列が、多彩な変化体のどれに合致するかの情報が提示されるが、今回の作業では、その情報は不要。

情報行にキリル文字列がある行は、
文字列変化の情報行。
その行には、その文字列に限った情報が書かれている。

キリル文字列がない行は、キリル文字列に関する情報のうち、
複数のキリル文字列に関する情報を提供する。

後になって、同類の情報が提供された場合には、上書きされるべき情報である。

```{r}
library(Dict)
hinshis <- c("名詞","代名詞","動詞","形容詞")
sex <- c("男性","女性","中性")

# 情報開始行の番数を返す
my.line.number.info <- function(str){
  sep1 <- unlist(strsplit(str,"\n"))
  hinshi.line <- 0
  
  for(i in 1:length(sep1)){
    if(is.element(sep1[i],hinshis)){
      hinshi.line <- i
      hinshi <- sep1[i]
      break
    }
  }
  return(i)
}


my.register.word <- function(str){
  sep1 <- unlist(strsplit(str,"\n"))
  hinshi <- sex <- ""
  hinshi.line <- 0
  
  for(i in 1:length(sep1)){
    if(is.element(sep1[i],hinshis)){
      hinshi.line <- i
      hinshi <- sep1[i]
    }
    if(is.element(sep1[i],sex)){
      sex <- sep1[i]
    }
  }
  sep2 <- list()
  for(i in 1:length(sep1)){
    sep2[[i]] <- unlist(strsplit(sep1[i]," "))
  }
}
```

```{}
名詞
女性
単数
主格 мера
生格 меры
与格 мере
対格 меру
造格 мерой
前置 мере
複数
主格 меры
生格 мер
与格 мерам
対格 меры
造格 мерами
前置 мерах
```

```{r}
in1 <- c("名詞
女性
単数
主格 мера
生格 меры
与格 мере
対格 меру
造格 мерой
前置 мере
複数
主格 меры
生格 мер
与格 мерам
対格 меры
造格 мерами
前置 мерах")
in1

```
```{r}
sep1 <- unlist(strsplit(in1,"\n"))
sep2 <- list()
for(i in 1:length(sep1)){
  sep2[[i]] <- unlist(strsplit(sep1[i]," "))
}
sep1
sep2
```