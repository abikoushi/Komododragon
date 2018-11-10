# Komododragon

## Installation

Install the latest version of this package from Github by pasting in the following.

~~~R
devtools::install_github("abikoushi/Komododragon")
~~~

## Simulation
~~~R
###シミュレーションでデータを生成
simfunc2 <- function(m,beta,alpha,X,tau){
  n <- nrow(X)
  X2 <-cbind(1,X)
  y <- rweibull(n,shape=m,scale=exp(X2%*%beta))
  z <- rbinom(n,1,plogis(X2%*%alpha))
  d <- (y < tau) & (z==1)
  y <- ifelse(d,y,tau)
  d <- as.integer(d)
  data.frame(y,X,d,z)
}

# n<-1000
# set.seed(1);X<-cbind(x1=runif(n,-1,1),x2=rbinom(n,1,0.5),x3=runif(n,-1,1))
# dat <-simfunc2(m=2,beta=c(1, 1,-1,0.5),alpha=c(2,1,-1,0.5),X=X,tau=3)

#y:来店間隔
#x1:説明変数
#x2:説明変数
#d:1-完全データ, 0-打切りデータ
#z:1-残存, 0-離脱 (本来は未観測の正解データ. 推定には使用しない)

mean(1-dat$d) #打切りデータの割合
mean(1-dat$z[dat$d==0]) #残存率

###最尤法でパラメータを推定

est1 <-estdropWeibull(formula1=y~x1+x2+x3,formula2=~x1+x2+x3,event=d,data = dat)

#formula1:来店間隔の式. (列名+列名+...+列名)の形で説明変数を選ぶ
#formula2:残存率の式

print(est1$opt$convergence) #0ならば収束

print(round(est1$opt$par,2)) #パラメータの推定値

print(-2*est1$opt$value+2*length(est1$opt$par)) #AIC

CIcalc(est1$opt) #パラメータの信頼区間

#残存している確率の推定値の分布
hist(est1$surv.prob[dat$d==0,],xlab = "probability",main="",breaks = "FD") 

###残存・離脱の答え合わせ
cl <-as.integer(est1$surv.prob[,1]>0.5)
conf_mat <- matrix(NA,2,2)
conf_mat[1,1]<-sum(cl[dat$d==0]==1 & (dat$z[dat$d==0]==1))/sum(dat$z[dat$d==0]==1)
conf_mat[1,2]<-sum(cl[dat$d==0]==1 & (dat$z[dat$d==0]==0))/sum(dat$z[dat$d==0]==0)
conf_mat[2,1]<-sum(cl[dat$d==0]==0 & (dat$z[dat$d==0]==1))/sum(dat$z[dat$d==0]==1)
conf_mat[2,2]<-sum(cl[dat$d==0]==0 & (dat$z[dat$d==0]==0))/sum(dat$z[dat$d==0]==0)

colnames(conf_mat)<-c("残存（正解）","離脱（正解）")
rownames(conf_mat)<-c("残存（予測）","離脱（予測）")

print(round(conf_mat,2))
~~~
