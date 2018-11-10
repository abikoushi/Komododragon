#' @export
estdropWeibull <- function(formula1,formula2,event,data=NULL,ini=NULL,maxit=5000,method="BFGS",reltol=1e-8){
  d_char <- as.character(substitute(event))
  d <- data[,d_char,drop=TRUE]
  data_f <-model.frame(formula1,data)
  y <- model.response(data_f)
  X1 <- model.matrix(formula1,data)
  X2 <- model.matrix(formula2,data)
  q1 <- ncol(X1)
  q2 <- ncol(X2)
  if(is.null(ini)){
    ini <- c(1,numeric(q1+q2))
  }
  names(ini) <-c("m",colnames(X1),colnames(X2))
  lldropWeibull <- function(par,y,d,X1,X2,q1,q2){
    m <-par[1]
    beta <-par[1:q1+1]
    alpha <- par[1:q2+q1+1]
    eta <- drop(exp(X1 %*% beta))
    logp <- drop(plogis(X2 %*% alpha,log.p=TRUE))
    log1mp <- drop(plogis(X2 %*% alpha,log.p=TRUE,lower.tail = FALSE))
    sum(d*(dweibull(y,shape=m,scale=eta,log = TRUE)+logp)+
          (1-d)*log(exp(log1mp)+exp(logp+pweibull(y,shape=m,scale=eta,lower.tail = FALSE,log.p = TRUE))))
  }
  grlldropWeibull <- function(par,y,d,X1,X2,q1,q2){
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
    dalpha_mat <-  d*X2/(1+exp(drop(X2%*%alpha)))+
      (1-d)*((B-1)*X2*exp(drop(X2%*%alpha)))/
      ((exp(drop(X2%*%alpha))+1)*(B*exp(drop(X2%*%alpha))+1))
    c(dm,apply(dbeta_mat,2,sum),apply(dalpha_mat,2,sum))
  }
  opt<-optim(ini,lldropWeibull,
             gr=grlldropWeibull,
             control = list(fnscale=-1,maxit=maxit,reltol=reltol),
             method = method,
             hessian = TRUE,
             y=y,d=d,X1=X1,X2=X2,q1=q1,q2=q2)
  m <-opt$par[1]
  beta <-opt$par[1:q1+1]
  alpha <- opt$par[1:q2+q1+1]
  eta <- drop(exp(X1 %*% beta))
  p <- drop(plogis(X2 %*% alpha))
  lp <- matrix(NA,length(y),2)
  lp[,1] <- ifelse(d==1,1,p*pweibull(y,shape=m,scale=eta,lower.tail = FALSE))
  lp[,2] <- ifelse(d==1,0,1-p)
  list(opt=opt,formula1=formula1,formula2=formula2,surv.prob=lp/apply(lp,1,sum))
}

#' @export
CIcalc <-function(opt,alpha=0.95){
  m1alpha <-1-alpha
  se <-sqrt(-diag(solve(opt$hessian)))
  upper <-qnorm(1-m1alpha/2,opt$par,se)
  lower <-qnorm(m1alpha/2,opt$par,se)
  data.frame(variavles=names(opt$par),estimate=opt$par,lower,upper)
}

