#################################################################
#Title: Approximation of the Posterior Distribution 
#of a Count Parameter in the Generalized Poisson and 
#the Negative Binomial Models
#
#Author: Hongxiang Li, Tsung Fei Khang
#Date: 25 February 2022
#Version: 1.0
#Contact: chelsea.divo@hotmail.com (HXL); tfkhang@um.edu.my (TFK)
#################################################################

library(VGAM)

#Function for computing the pmf of the posterior distribution of k
#with GP or NB models, for k up to s

pospmf <- function(observed_count, b, m, model = 1, opt = 0, s=10^4){
  # observed_count = 0 or 1, exact 
  # observed_count >= 2, approximation
  # model = 1 (GP)
  # model = 2 (NB)
  # opt = 0 (exact/true post_dist.) 
  # opt = 1 (approximation)
  # opt = 2 (adjusted form of approximation, under GP)
  # 0<b<1
  # 0<m<min{4, 1/(1-b)^2}, under GP
  # 0<m<1, under NB
  # s is the largest count argument for the posterior pmf
  
  if (s < 10^4){
    s1 <- 10^4
  } 
  else {
    s1 <- s
  }
  
  
  if (observed_count < 0){
     print("Observed count should be non-negative.")
  }
  
  else if (model == 1){
    x <- observed_count
    k <- x:s1
    
    if (x == 0 || x == 1){
           if (opt != 0){
              print("This ia an exact posterior distribution.")
           }
        opt <- 0
    }

    if (b >= 1 || b <=0 || m >= 4 || m >= (1-b)^{-2}){
      print("Entered wrong parameter.")
    } 
    
    else if (opt == 0){ # exact posterior distribution under GP
           
                 up <- k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m))
                 down <- sum ( k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m)) )
                 
                 if(x==0){
                   return((1-exp(-b*sqrt(m)))*(exp(-b*k*sqrt(m)))[1:(s+1)])
                 }
                 else if (x==1){
                   return( k*(1-exp(-b*sqrt(m)))^2*(exp(-b*(k-1)*sqrt(m)))[1:s] )
                 }
                 else{
                   p <- c(rep(0,x-1), up/down)[1:s]
                   return(p)
                 }
    }
    
    else if (opt ==1) { # approximate posterior distribution under GP
      K <- x+1
      theta <- b*sqrt(m)
      return (dgamma(1:s, shape=K, scale=1/theta) )
    } 
    
    else if (opt == 2) { # adjusted form of approximation under GP
        up <- k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m))
        down <- sum ( k * (b*k * sqrt(m) + x*(1-sqrt(m)))^(x-1) * exp(-b*k*sqrt(m)) )
        p <- c(rep(0,x-1), up/down)[1:s1]
        
        mu <- sum ( p*(1:s) ) 
        sigma.sq <- sum(p*(1:s)^2) - mu^2
        
        K <- mu^2 / sigma.sq
        theta <- sigma.sq / mu
        
        return (dgamma(1:s, shape=K, scale=theta))
    }
    else {
      print("Entered the wrong form of distribution under GP.")
    } 
  } 
  
  else if (model == 2){
    
    d_nb <- function(y,r,p){
      d_nb <- dnbinom(y, r, 1-p)
      return(d_nb)
    }
    
    y <- observed_count
    k <- y:s1
    
    if (y == 0 || y == 1){
        if (opt != 0){
          print("This ia an exact posterior distribution.")
        }
        opt <- 0
    }
    
    if (b <= 0 || b >= 1 || m <= 0 || m >= 1){
      print("Entered wrong parameter.")
    }
    
    else if (opt == 0){ # exact/true post_dist., under NB
      r <- k*b*m/(1-m)
      p <- 1-m
      
      up <- d_nb(y,r,p)
      down <- sum( d_nb(y,r,p) )
      
      if (y==0){
        return( c(up/down)[1:(s+1)])
      }
      else if (y==1){
        return( c(up/down)[1:s])
      }
      else{
        return( c(rep(0,y-1), up/down)[1:s])
      }
      
    } 
    
    else if (opt == 1){ # approximate posterior distribution under NB
      t. <- b*m/(1-m)
      A <- (1/{2*t.}+1)*y - 1/{2*t.}
      
      up <- (m^t.)^{k-y}*(k+A-y)^{y}
      L <- lerch(m^t., -y, A, tolerance = 1.0e-10, iter = 100)
      
      return( c(rep(0,y-1), up/L)[1:s])
    } 
    else {
      print("No adjusted form of distribution under NB.")
    }
  } 
  else {
    print("GP and NB models only.")
  }
}


###############################
### KL Divergence - GP model
###############################

# model = 1 - GP
KL_D_GP <- function(x, b, m, opt=1){
  # opt = 1: exact vs. approximation
  # opt = 2: exact vs. adjusted form of approximation
  if(opt==1){
    n <- sum(pospmf(x,b,m,1,0) > 10^{-15})
    return(sum(pospmf(x,b,m,1,0)[x:(n+x-1)]*log(pospmf(x,b,m,1,0)[x:(n+x-1)]/pospmf(x,b,m,1,1)[x:(n+x-1)])))
  } else if(opt==2){
    n <- sum(pospmf(x,b,m,1,0) > 10^{-15})
    return(sum( pospmf(x,b,m,1,0)[x:(n+x-1)]*log(pospmf(x,b,m,1,0)[x:(n+x-1)]/pospmf(x,b,m,1,2)[x:(n+x-1)]) ))
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
T1  # Table 1




##############################  
### KL Divergence - NB model
##############################

# model = 2 - NB
D_KL_NB <- function(y, b, m){
  # exact Vs. approximation
  n <- sum(pospmf(y,b,m,2,0) > 10^{-15})
  D_KL <- sum( pospmf(y,b,m,2,0)[y:(n+y-1)]*log(pospmf(y,b,m,2,0)[y:(n+y-1)]/pospmf(y,b,m,2,1)[y:(n+y-1)]))
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
T2 # Table 2


#Usage example
#
#excellent fit
#prob_exact <- pospmf(observed_count=50, b=0.1, m=1.2, opt=0)
#prob_approx <- pospmf(observed_count=50, b=0.1, m=1.2, opt=1)
#prob_adj <- pospmf(observed_count=50, b=0.1, m=1.2, opt=2)
#
#h <- which(prob_exact == max(prob_exact))

#plot(prob_exact[1:(2*h)], type="l", xlab="k", ylab="Probability")
#points(prob_approx[1:(2*h)], type="l", col="red")
#points(prob_adj[1:(2*h)], type="l", col="blue")
#legend(2*h-sqrt(2*h*100),max(prob_exact), lty=rep(1,3), 
#col=c("black","red","blue"), c("Exact","Approximation","Adjusted"), cex=0.7)

#poorer fit
#prob_exact <- pospmf(observed_count=2, b=0.5, m=3, opt=0)
#prob_approx <- pospmf(observed_count=2, b=0.5, m=3, opt=1)
#prob_adj <- pospmf(observed_count=2, b=0.5, m=3, opt=2)

#h <- which(prob_exact == max(prob_exact))

#plot(prob_exact[1:(10*h)], type="l", xlab="k", ylab="Probability")
#points(prob_approx[1:(10*h)], type="l", col="red")
#points(prob_adj[1:(10*h)], type="l", col="blue")
#legend(10*h-sqrt(10*h*2.5),max(prob_exact), lty=rep(1,3), 
#col=c("black","red","blue"), c("Exact","Approximation","Adjusted"), cex=0.7)

#function for checking m,p values lie in valid parameter space
parcheck <- function(bpar, mpar){

	b <- seq(0, 0.5, 0.01)
	m <- 1/(1-b)^2

	plot(b, m, type = "l", xlim = c(0, 1), ylim = c(0,4), col="white")
	polygon(c(0,b,0.5), c(0,m,0), col="gray", border=NA)
	points(b, m, type = "l", lty=2)
	halfx <- seq(0.5,1,0.01)
	polygon(c(0.5,halfx,1), c(0, rep(4, length(halfx)), 0), col="gray", border=NA)
	
	lines(c(0,0),c(0,1), lty=2)
	lines(c(0,1),c(0,0), lty=2)
	lines(c(1,1),c(0,4), lty=2)
	lines(c(0.5,1),c(4,4), lty=2)

	points(bpar, mpar, pch=4)

	if(mpar >=4 | mpar <= 0 | bpar <= 0 | bpar > 1 | mpar > (1-bpar)^(-2) )
	print("Out of parameter space.")
	else print("Valid parameter values")
	}

#examples:
#parcheck(0.5, 3) #valid parameters
#parcheck(0.05, 2)    #invalid parameters

