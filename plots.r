source("black-scholes.r")

if (!exists("persp.orig")) persp.orig <- persp

persp <- function(t,s,v, expand=0.45, ...) { 
  nrz <- nrow(v) 
  ncz <- ncol(v)
	nb.col <- 256
	color <- heat.colors(nb.col) 
	facet <- - (v[-1, -1] + v[-1, -ncz] + v[-nrz, -1] + v[-nrz, -ncz])
	facetcol <- cut(facet, nb.col) 
	opar <- par(mai=c(0.2,0.2,0.2,0.2))

	hs <- head(s,1)
	ts <- tail(s,1)
	ms <- 100
	ht <- head(t,1)
	tt <- tail(t,1)
	hv <- max(v,na.rm=TRUE)
  lv <- min(v, na.rm=TRUE)
        
  persp.orig(t, s, v, col=color[facetcol], expand=expand, axes=FALSE, ...) -> res
  
	lines(trans3d(c(ht,ht),c(hs,ts),c(lv,lv), res), lty='dotted')
	lines(trans3d(c(ht,ht,ht,tt,tt,tt),c(ms,ms,ms,ms,ms,ms),c(lv,hv,hv,hv,hv,lv), res), lty='dotted')
   
  text(trans3d(c(tt,ht),c(ts,ts),c(lv,lv), res), labels=c("0","T"), pos=4, offset=.5)
  text(trans3d(c(tt),c(ms),c(lv), res), labels=c("K"), pos=1, offset=.5)

	par(opar)
}

S <- 100
K <- 100
T <- 1
N <- 10
r <- 0.10
sigma <- 0.15

s <- sort(c(K-EPS/10, K+EPS/10,seq(0.5*K,1.5*K,length=51)))
t <- sort(c(EPS/10, seq(0,T, length=51)))

v <- outer(t,s,Vectorize(function(t, s) bs.value("ce",s,K,t,r,sigma)))
persp(t,s,v,zlab="Value",theta=120,phi=30,zlim = c(0,1))

v <- outer(t,s,Vectorize(function(t, s) bs.delta("ce",s,K,t,r,sigma)))
persp(t,s,v,zlab="Delta Call", zlim = c(-1,1),theta=120,phi=30)

v <- outer(t,s,Vectorize(function(t, s) bs.delta("pe",s,K,t,r,sigma)))
persp(t,s,v,zlab="Delta Put", zlim = c(-1,1),theta=120,phi=30)

v <- outer(t,s,Vectorize(function(t, s) bs.gamma("ce",s,K,t,r,sigma)))
persp(t,s,v,zlab="Gamma Call",theta=120,phi=30)

v <- outer(t,s,Vectorize(function(t, s) bs.theta("ce",s,K,t,r,sigma)))
persp(t,s,v,zlab="Theta Call",theta=120,phi=30)

v <- outer(t,s,Vectorize(function(t, s) bs.theta("pe",s,K,t,r,sigma)))
persp(t,s,v,zlab="Theta Put",theta=120,phi=30)

v <- outer(t,s,Vectorize(function(t, s) bs.vega("ce",s,K,t,r,sigma)))
persp(t,s,v,zlab="Vega Call",theta=120,phi=30)

v <- outer(t,s,Vectorize(function(t, s) bs.rho("ce",s,K,t,r,sigma)))
persp(t,s,v,zlab="Rho Call",theta=120,phi=30)


png(file="bs-delta-call.png",width=480,height=480)
persp(t,s,v,zlab="Delta Call", main="Black-Scholes European Call Delta", zlim = c(-1,1),theta=120,phi=30)
dev.off()