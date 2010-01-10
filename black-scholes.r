# Black-Scholes model: value function and some greeks
# Dale Roberts (c) 2010

#	S = Asset price
# K = Strike price
# T = Time to maturity
# r = Risk-free rate
# sigma = Volatility of asset price

EPS <- 0.01

bs.value <- function(type = c("pe","ce"), S, K, T, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  if (type == "ce") { # European call
    if (T < EPS) { # Small time to expiry
      return(max(S-K,0))
    } else {
      return(S*pnorm(d1) - K*exp(-r*(T))*pnorm(d2))
    }
  } else {
    if (T < EPS) { # Small time to expiry
      return(max(K-S,0))
    } else {
      return(K*exp(-r*(T))*pnorm(-d2) - S*pnorm(-d1))
    }
  }
}

bs.delta <- function(type = c("pe","ce"), S, K, T, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
  if (type == "ce") {
    # by hand for small T
    if (T < EPS) 
    {
      if (K-S < 0) {
        return(1.0)
      } else {
        return(0.0)
      }
    } 
    # formula
    else {
      return(pnorm(d1))
    }
  } else { # Put
    if (T < EPS) {
      if (K-S > 0) {
        return(-1.0)
      } else {
        return(0.0)
      }
    } 
    else {
      return(pnorm(d1)-1)
    }
  }
}

bs.gamma <- function(type = c("pe","ce"), S, K, T, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
  return(dnorm(d1)/(S*sigma*sqrt(T)))
}

bs.theta <- function(type = c("pe","ce"), S, K, T, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  if (type == "ce")
    return(-S*dnorm(d1)*sigma/(2*sqrt(T)) - r*K*exp(-r*T)*pnorm(d2))
  else
    return(-S*dnorm(d1)*sigma/(2*sqrt(T)) + r*K*exp(-r*T)*pnorm(-d2))
}

bs.vega <- function(type = c("pe","ce"), S, K, T, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
  return(S*sqrt(T)*dnorm(d1))
}

bs.rho <- function(type = c("pe","ce"), S, K, T, r, sigma) {
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  if (type == "ce")
    return(K*T*exp(-r*T)*pnorm(d2))
  else
    return(-K*T*exp(-r*T)*pnorm(-d2))
}
