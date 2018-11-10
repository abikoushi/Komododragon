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
#dropout_weibull <-stan_model("~/dropbox/コモドドラゴン/dropout_weibull.stan")
#dat4stan <- list(n=nrow(dat),q=ncol(X)+1,y=dat$y,d=dat$d,X=cbind(1,X))
#smp1 <-sampling(dropout_weibull,dat4stan,iter=1)
#traceplot(smp1,"alpha")
###最尤法でパラメータを推定
library("numDeriv")

est1 <-estdropWeibull(formula1=y~x1+x2+x3,formula2=~x1+x2+x3,event=d,data = dat,method = "BFGS")
est2 <-estdropWeibull2(formula1=y~x1+x2+x3,formula2=~x1+x2+x3,event=d,data = dat)

est1$opt$par
est2$opt$par
est1$opt$value
est2$opt$value
formula1=y~x1+x2
formula2=~x1+x2
d=dat$d
data=dat
data_f <-model.frame(formula1,data)
y <- model.response(data_f)
X1 <- model.matrix(formula1,data)
X2 <- model.matrix(formula2,data)
q1 <- ncol(X1)
q2 <- ncol(X2)
lldropWeibull <- function(par){
  m <-par[1]
  beta <-par[1:q1+1]
  alpha <- par[1:q2+q1+1]
  eta <- drop(exp(X1 %*% beta))
  logp <- drop(plogis(X2 %*% alpha,log=TRUE))
  log1mp <- drop(plogis(X2 %*% alpha,log=TRUE,lower.tail = FALSE))
  sum(d*(dweibull(y,shape=m,scale=eta,log = TRUE)+logp)+
        (1-d)*log(exp(log1mp)+exp(logp+pweibull(y,shape=m,scale=eta,lower.tail = FALSE,log.p = TRUE))))
}
grlldropWeibull <- function(par){
  m <-par[1]
  beta <-par[1:q1+1]
  alpha <- par[1:q2+q1+1]
  eta <- drop(exp(X1 %*% beta))
  logp <- drop(plogis(X2 %*% alpha,log=TRUE))
  log1mp <- drop(plogis(X2 %*% alpha,log=TRUE,lower.tail = FALSE))
  p <- drop(plogis(X2 %*% alpha))
  dm <- sum(
    d*(1/m+log(y)-drop(X1%*%beta)+ ((y/eta)^m)*(-(log(y)-drop(X1 %*% beta))))-
      (1-d)*(p*exp(-(y/eta)^m))*((y/eta)^m)*(log(y)-drop(X1%*%beta))/
      (p*exp(-(y/eta)^m)-p+1)
  )
  dbeta_mat <- d*(m*X1*((y/eta)^m-1))+
    (1-d)*(m*p*X1*y*exp(-(y/eta)^m-drop(X1%*%beta))*((y/eta)^(m-1)))/
    (p*exp(-(y/eta)^m)-p+1)
  B <- pweibull(y,shape=m,scale=eta,lower.tail = FALSE)
  dalpha_mat <- d*X2/(1+exp(drop(X2%*%alpha)))+
    (1-d)*(X2/(exp(drop(X2%*%alpha))+1)-X2/(B*exp(drop(X2%*%alpha))+1))
  c(dm,apply(dbeta_mat,2,sum),apply(dalpha_mat,2,sum))
}
grad(lldropWeibull,rep(1,7))
grlldropWeibull(rep(1,7))
curve(1/(1+exp(2*x)),-2,2)
curve(plogis(2*x,lower.tail = FALSE),add = TRUE,col="red")
#formula1:来店間隔の式. (列名+列名+...+列名)の形で説明変数を選ぶ
#formula2:残存率の式

est1$opt$convergence #0ならば収束

round(est1$opt$par,2) #パラメータの推定値

-2*est1$opt$value+2*length(est1$opt$par) #AIC

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

round(conf_mat,2)
~~~
