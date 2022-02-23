###############################################
### Approximation of Posterior Distribution ###
### of A Count Parameter in GP & NB Model   ###
###############################################
library(VGAM)

### Approximation of Posterior Distribution

P_D <- function(observed_count, b, m, prior = 1, opt = 0 ){
  # observed_count >= 2
  # prior = 1 (GP)
  # prior = 2 (NB)
  # opt = 0 (exact/true post_dist.) 
  # opt = 1 (approximation)
  # opt = 2 (adjusted form of approximation)
  # 0<b<1
  # 0<m<min{4, 1/(1-b)^2}, under GP
  # 0<m<1, under NB
  if (prior == 1 ){
    x <- observed_count
    s <- 100000
    k <- x:s
    
    if (opt == 0){ # exact post_dist. under GP
       up <- k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m))
       down <- sum ( k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m)) )
      
       p <- c(rep(0,x-1), up/down)
       return(p)
    } else if (opt ==1 ) { # approximated post_dist., under GP
      K <- x+1
      theta <- b*sqrt(m)
      return (dgamma(1:s, shape=K, scale=1/theta) )
    } else if (opt ==2 ) { # adjsut form of approximation, under GP
      up <- k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m))
      down <- sum ( k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m)) )
      
      p <- c(rep(0,x-1), up/down)
      
      mu <- sum ( p*(1:s) ) 
      sigma.sq <- sum(p*(1:s)^2) - mu^2
      
      K <- mu^2 / sigma.sq
      theta <- sigma.sq / mu
      
      return (dgamma(1:s, shape=K, scale=theta))
    } else {
      print("Entered the wrong form of approximation under GP.")
    } 
  }else if ( prior == 2) {
    
    d_nb <- function(y,r,p){
      d_nb <- dnbinom(y, r, 1-p)
      return(d_nb)
    }
    
    y <- observed_count
    s <- 10000
    k <- y:s
    
    
    if (opt == 0){ # exact/true post_dist., under NB
      t <- b*m/(1-m)
      r <- k*t
      p <- 1-m
      
      up <- d_nb(y,r,p)
      down <- sum( d_nb(y,r,p) )
      return( c(rep(0,y-1), up/down) )
      
    } else if (opt == 1){ # Approximatation of post_dist., under NB
      t <- b*m/(1-m)
      A <- (1/{2*t}+1)*y - 1/{2*t}
      
      up <- (m^t)^{k-y}*(k+A-y)^{y}
      L <- lerch(m^t, -y, A, tolerance = 1.0e-10, iter = 100)
      return( c(rep(0,y-1), up/L) )
      
    } else {
      print("No adjusetd or other form of approximation under NB.")
    }
  } else {
    print("Entered the wrong prior distribution.")
  }
}




###############################
### KL Divergence, under GP ###
###############################

# prior = 1, under GP
KL_D_GP <- function(observed_count, b, m, opt=1){
  # opt = 1: exact Vs. approximation
  # opt = 2: exact Vs. adjusted form of approximation
  x <- observed_count
  if(opt==1){
    n <- sum(P_D(x,b,m,1,0) > 10^{-15})
    return(sum(P_D(x,b,m,1,0)[x:(n+x-1)]*log(P_D(x,b,m,1,0)[x:(n+x-1)]/P_D(x,b,m,1,1)[x:(n+x-1)])))
  } else if(opt==2){
    n <- sum(P_D(x,b,m,1,0) > 10^{-15})
    return(sum( P_D(x,b,m,1,0)[x:(n+x-1)]*log(P_D(x,b,m,1,0)[x:(n+x-1)]/P_D(x,b,m,1,2)[x:(n+x-1)]) ))
  }
}


b <- c(rep(0.1, 4), rep(0.2, 4), rep(0.5,4))
m0 <- c(0.1,0.2,0.5)
m <- c(m0, 1.2, m0, 1.5, m0, 3)
X <- c(2,2,3,3,10,10,20,20,50,50)

M_gp <- matrix(ncol = 10, nrow = 12)
for (i in 1:12) {
  for (j in 1:10) {
    if (j%% 2 == 0){
      opt. <- 2
    } else {
      opt. <- 1
    }
    M_gp[i,j] <- KL_D_GP(X[j], b[i], m[i], opt.)
  }
}

T1 <- round(M_gp, 4)
T1  # matrix for Talbe 1
  
##############################  
### KL Divergence, under NB###
##############################

# prior = 2, under NB
D_KL_NB <- function(observed_count, b, m){
  # exact Vs. approximation
  y <- observed_count
  n <- sum(P_D(y,b,m,2,0) > 10^{-15})
  
  D_KL <- sum( P_D(y,b,m, 2, 0)[y:(n+y-1)]*log(P_D(y,b,m, 2, 0)[y:(n+y-1)]/P_D(y,b,m, 2, 1)[y:(n+y-1)]))
  return(D_KL)
}

b <- c(rep(c(0.1,0.2,0.5), 4))
m <- c(rep(0.2, 3), rep(0.4, 3),rep(0.6, 3),rep(0.8, 3) )
Y <- c(2,3,10,20,50)

M_nb <- matrix(ncol = 5, nrow = 12)
for (i in 1:12) {
  for (j in 1:5) {
    M_nb[i,j] <- D_KL_NB(Y[j], b[i], m[i])
  }
}

T2 <- round(M_nb, 5)
T2 # matrix for Table 2
