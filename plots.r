# Plots of Black-Scholes Model

source("black-scholes.r")

heat.persp <- function(t,s,v,expand=0.45,...) { 
  nrz <- nrow(v) 
  ncz <- ncol(v)
  nb.col <- 256
  color <- heat.colors(nb.col) 
  facet <- - (v[-1, -1] + v[-1, -ncz] + v[-nrz, -1] + v[-nrz, -ncz])
  facetcol <- cut(facet, nb.col) 

  persp(t,s,v,col=color[facetcol],expand=expand,...) -> res

  return(res)
}
  
the.lines <- function(t,s,v,res) {
  hs <- head(s,1)
  ts <- tail(s,1)
  ms <- 100
  ht <- head(t,1)
  tt <- tail(t,1)
  hv <- max(v,na.rm=TRUE)
  lv <- min(v,na.rm=TRUE)

  lines(trans3d(c(ht,ht,ht,tt,tt,tt),
                c(ms,ms,ms,ms,ms,ms),
                c(lv,hv,hv,hv,hv,lv),
                res),
        lty='dashed')
   
  text(trans3d(c(tt),
               c(ms),
               c(lv),
               res),
       labels=c("S=K"),
       pos=1,
       offset=.5)
}

S <- 100
K <- 100
N <- 10
r <- 0.10
tau <- 1
sigma <- 0.15

s <- sort(c(K-EPS/10,K+EPS/10,seq(0.5*K,1.5*K,length=31)))
t <- sort(c(EPS/10,seq(0,tau,length=31)))

opar <- par(mar=rep(0,4))

doit <- function(f, filename, name) {
  fn <- paste("images/", filename, sep="")
  png(file=fn,width=480,height=480)
  v <- outer(t,s,Vectorize(function(t, s) f(s,K,t,r,sigma)))
  heat.persp(t,s,v,
             zlab="",
             xlab="TTM",
             ylab="S",
             theta=108,
             phi=40,
             axes=FALSE) -> res
  the.lines(t,s,v,res)
  title(name, line=0)
  dev.off()
  system(paste("convert -trim", fn, fn))
}

doit(bscall.value, "bscall-value.png", "Black-Scholes Call Value")
doit(bscall.delta, "bscall-delta.png", "Black-Scholes Call Delta")
doit(bscall.gamma, "bscall-gamma.png", "Black-Scholes Call Gamma")
doit(bscall.vega, "bscall-vega.png", "Black-Scholes Call Vega")
doit(bscall.rho, "bscall-rho.png", "Black-Scholes Call Rho")

par(opar)
