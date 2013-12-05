# Black-Scholes model: value function and some greeks
# Dale Roberts (c) 2010

# S = Asset price
# K = Strike price
# r = Risk-free rate
# tau = Time to maturity
# sigma = Volatility of asset price

EPS <- 1./365. 

bscall.value <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
  if (tau < EPS) { # Small time to expiry
    return(max(S-K,0))
  } else {
    return(S*pnorm(d1) - K*exp(-r*(tau))*pnorm(d2))
  }
}

bsput.value <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
   if (tau < EPS) { # Small time to expiry
    return(max(K-S,0))
  } else {
    return(K*exp(-r*(tau))*pnorm(-d2) - S*pnorm(-d1))
  }
}

bscall.delta <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  if (tau < EPS) { # by hand for small time to expiry
    if (K-S < 0) {
      return(1.0)
    } else {
      return(0.0)
    }
  } else {
    return(pnorm(d1))
  }
}

bsput.delta <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  if (tau < EPS) { # by hand for small time to expiry
    if (K-S > 0) {
      return(-1.0)
    } else {
      return(0.0)
    }
  } else {
    return(pnorm(d1)-1)
  }
}

bscall.gamma <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  return(dnorm(d1)/(S*sigma*sqrt(tau)))
}

bscall.theta <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
  return(-S*dnorm(d1)*sigma/(2*sqrt(tau)) - r*K*exp(-r*tau)*pnorm(d2))
}

bsput.theta <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
  return(-S*dnorm(d1)*sigma/(2*sqrt(tau)) + r*K*exp(-r*tau)*pnorm(-d2))
}

bscall.vega <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  return(S*sqrt(tau)*dnorm(d1))
}

bscall.rho <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
  return(K*tau*exp(-r*tau)*pnorm(d2))
}

bsput.rho <- function(S, K, tau, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
  return(-K*tau*exp(-r*tau)*pnorm(-d2))
}
