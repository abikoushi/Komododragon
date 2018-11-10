\name{estdropWeibull}
\alias{estdropWeibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estdropWeibull(formula1, formula2, event, data = NULL, ini = NULL, maxit = 5000, method = "BFGS", reltol = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{
%%     ~~Describe \code{formula1} here~~
}
  \item{formula2}{
%%     ~~Describe \code{formula2} here~~
}
  \item{event}{
%%     ~~Describe \code{event} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{ini}{
%%     ~~Describe \code{ini} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{reltol}{
%%     ~~Describe \code{reltol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula1, formula2, event, data = NULL, ini = NULL,
    maxit = 5000, method = "BFGS", reltol = 1e-08)
{
    d_char <- as.character(substitute(event))
    d <- data[, d_char, drop = TRUE]
    data_f <- model.frame(formula1, data)
    y <- model.response(data_f)
    X1 <- model.matrix(formula1, data)
    X2 <- model.matrix(formula2, data)
    q1 <- ncol(X1)
    q2 <- ncol(X2)
    if (is.null(ini)) {
        ini <- c(1, numeric(q1 + q2))
    }
    names(ini) <- c("m", colnames(X1), colnames(X2))
    lldropWeibull <- function(par, y, d, X1, X2, q1, q2) {
        m <- par[1]
        beta <- par[1:q1 + 1]
        alpha <- par[1:q2 + q1 + 1]
        eta <- drop(exp(X1 \%*\% beta))
        logp <- drop(plogis(X2 \%*\% alpha, log.p = TRUE))
        log1mp <- drop(plogis(X2 \%*\% alpha, log.p = TRUE, lower.tail = FALSE))
        sum(d * (dweibull(y, shape = m, scale = eta, log = TRUE) +
            logp) + (1 - d) * log(exp(log1mp) + exp(logp + pweibull(y,
            shape = m, scale = eta, lower.tail = FALSE, log.p = TRUE))))
    }
    grlldropWeibull <- function(par, y, d, X1, X2, q1, q2) {
        m <- par[1]
        beta <- par[1:q1 + 1]
        alpha <- par[1:q2 + q1 + 1]
        eta <- drop(exp(X1 \%*\% beta))
        logp <- drop(plogis(X2 \%*\% alpha, log = TRUE))
        log1mp <- drop(plogis(X2 \%*\% alpha, log = TRUE, lower.tail = FALSE))
        p <- drop(plogis(X2 \%*\% alpha))
        dm <- sum(d * (1/m + log(y) - drop(X1 \%*\% beta) + ((y/eta)^m) *
            (-(log(y) - drop(X1 \%*\% beta)))) - (1 - d) * (p *
            exp(-(y/eta)^m)) * ((y/eta)^m) * (log(y) - drop(X1 \%*\%
            beta))/(p * exp(-(y/eta)^m) - p + 1))
        dbeta_mat <- d * (m * X1 * ((y/eta)^m - 1)) + (1 - d) *
            (m * p * X1 * y * exp(-(y/eta)^m - drop(X1 \%*\% beta)) *
                ((y/eta)^(m - 1)))/(p * exp(-(y/eta)^m) - p +
            1)
        B <- pweibull(y, shape = m, scale = eta, lower.tail = FALSE)
        dalpha_mat <- d * X2/(1 + exp(drop(X2 \%*\% alpha))) +
            (1 - d) * ((B - 1) * X2 * exp(drop(X2 \%*\% alpha)))/((exp(drop(X2 \%*\%
                alpha)) + 1) * (B * exp(drop(X2 \%*\% alpha)) +
                1))
        c(dm, apply(dbeta_mat, 2, sum), apply(dalpha_mat, 2,
            sum))
    }
    opt <- optim(ini, lldropWeibull, gr = grlldropWeibull, control = list(fnscale = -1,
        maxit = maxit, reltol = reltol), method = method, hessian = TRUE,
        y = y, d = d, X1 = X1, X2 = X2, q1 = q1, q2 = q2)
    m <- opt$par[1]
    beta <- opt$par[1:q1 + 1]
    alpha <- opt$par[1:q2 + q1 + 1]
    eta <- drop(exp(X1 \%*\% beta))
    p <- drop(plogis(X2 \%*\% alpha))
    lp <- matrix(NA, length(y), 2)
    lp[, 1] <- ifelse(d == 1, 1, p * pweibull(y, shape = m, scale = eta,
        lower.tail = FALSE))
    lp[, 2] <- ifelse(d == 1, 0, 1 - p)
    list(opt = opt, formula1 = formula1, formula2 = formula2,
        surv.prob = lp/apply(lp, 1, sum))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
