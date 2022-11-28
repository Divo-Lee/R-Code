#################################################################
#Title: Some approximation results for Bayesian posteriors
#that involve the Hurwitz-Lerch zeta distribution
#
#
#Authors: Hongxiang Li, Tsung Fei Khang
#Date: 28 November 2022
#Version: 1.1
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
    else if(opt == 1) {
      # approximate posterior distribution under GP
      t. <- exp(-b*sqrt(m))
      g <- (1-sqrt(m))/(b*sqrt(m))
      up <- (t.)^(k-x)*(k + g*x)^x
      L <- lerch(t., -x, (g*x + x), tolerance = 1.0e-10, iter = 100)
      return(c(rep(0,x-1), up/L)[1:s])
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
KL_D_GP <- function(x, b, m){
    n <- sum(pospmf(x,b,m,1,0) > 10^{-15})
    return(sum(pospmf(x,b,m,1,0)[x:(n+x-1)]*log(pospmf(x,b,m,1,0)[x:(n+x-1)]/pospmf(x,b,m,1,1)[x:(n+x-1)])))
}


## approximated as extended Hurwitz-Lerch zeta under GP
b <- c(rep(0.2, 6), rep(0.4, 6), rep(0.6,6), rep(0.8,6))
m0 <- c(0.2, 0.4, 0.6, 0.8, 1.2)
m <- c(m0, 1.5, m0, 2, m0, 2.5, m0, 3)
X <- c(2,3,10,20,50)

M_gp <- matrix(ncol = 5, nrow = 24)
for (i in 1:24) {
  for (j in 1:5) {
    M_gp[i,j] <- KL_D_GP(X[j], b[i], m[i])
  }
}

T1 <- round(M_gp, 4)
T1 # Table 1
# maybe get NA value for some inputs, 
# because the VGAM lerch() function get NA for some large parameters and variables


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


b <- c(rep(c(0.2,0.4,0.6,0.8), 4))
m <- c(rep(0.2, 4), rep(0.4, 4),rep(0.6, 4),rep(0.8, 4))
Y <- c(2,3,10,20,50)

M_nb <- matrix(ncol = 5, nrow = 16)
for (i in 1:16) {
  for (j in 1:5) {
    M_nb[i,j] <- D_KL_NB(Y[j], b[i], m[i])
  }
}

T2 <- round(M_nb, 4)
T2 # Table 2



#######################
### Posterior Mean - GP
#######################

# model 1 - GP

mean_post_gp <- function(x, b, m, opt=0){
  #opt = 0, exact
  #opt = 1, approximation 1, with Lerch ()
  #opt = 2, approximation 2
  if (opt == 0){
    mu <- sum(pospmf(x, b, m, model = 1, opt =0)*seq(1,10000))
    return(mu) # exact
  }
  else if (opt == 1){
    theta <- exp(-b*sqrt(m))
    w <- 1 + (1 - sqrt(m))/(b*sqrt(m))
    L1 <- lerch(theta, -x-1, w*x, tolerance = 1.0e-10, iter = 100)
    L2 <- lerch(theta, -x, w*x, tolerance = 1.0e-10, iter = 100)
    mu <- L1/L2 - x*(1 - sqrt(m))/(b*sqrt(m))
    return(mu)  # approximation 1, Eq.(18)
  } else if (opt == 2){
    mu <- x/b + 1/(b*sqrt(m))
    return(mu) # approximation 2, Theorem 2
  }
}

b <- c(rep(0.2, 6), rep(0.4, 6), rep(0.6,6), rep(0.8,6))
m0 <- c(0.2, 0.4, 0.6, 0.8, 1.2)
m <- c(m0, 1.5, m0, 2, m0, 2.5, m0, 3)
x <- c(2,2,2,3,3,3,10,10,10,20,20,20,50,50,50)


## Comparison 
Tab3 <- matrix(ncol = 15, nrow = 24)
for (i in 1:24) {
  for (j in 1:15) {
    if (j %% 3 == 1){
      opt <- 0
      Tab3[i,j] <- mean_post_gp(x[j], b[i], m[i], opt)
    } else if (j %% 3 == 2) {
      opt <- 1
      Tab3[i,j] <- mean_post_gp(x[j], b[i], m[i], opt) - Tab3[i,(j-1)]
    } else if (j %% 3 == 0) {
      opt <- 2
      Tab3[i,j] <- mean_post_gp(x[j], b[i], m[i], opt) - Tab3[i,(j-2)]
    }
  }
}

Tab3 <- round(Tab3, 1)
Tab3 # Table 3



###########################
### Posterior Variance - GP
###########################

# model 1 - GP

var_post_gp <- function(x, b, m, obt = 0){
  #opt = 0, exact
  #opt = 1, approximation 1, with Lerch ()
  #opt = 2, approximation 2
  if (opt == 0){
    mu <- sum(pospmf(x, b, m, model = 1, opt =0)*seq(1,10000))
    mu2 <- sum(pospmf(x, b, m, model = 1, opt =0)*(seq(1,10000))^2)
    var <- mu2 - mu^2 
    return(var) # exact
  } else if (opt == 1){
    theta <- exp(-b*sqrt(m))
    w <- 1 + (1 - sqrt(m))/(b*sqrt(m))
    L1 <- lerch(theta, -x-1, w*x, tolerance = 1.0e-10, iter = 100)
    L2 <- lerch(theta, -x, w*x, tolerance = 1.0e-10, iter = 100)
    L3 <- lerch(theta, -x-2, w*x, tolerance = 1.0e-10, iter = 100)
    var <- L3/L2 -(L1/L2)^2 
    return(var) # approximation 1, Eq.(19)
  } else if (opt == 2){
    var_new <- (x+1)/(m*b^2) 
    return(var_new) # approximation 2, Theorem 2
  }
}

b <- c(rep(0.2, 6), rep(0.4, 6), rep(0.6,6), rep(0.8,6))
m0 <- c(0.2, 0.4, 0.6, 0.8, 1.2)
m <- c(m0, 1.5, m0, 2, m0, 2.5, m0, 3)
x <- c(2,2,2,3,3,3,10,10,10,20,20,20,50,50,50)

Tab4 <- matrix(ncol = 15, nrow = 24)
for (i in 1:24) {
  for (j in 1:15) {
    if (j %% 3 == 1){
      opt <- 0 
      Tab4[i,j] <- var_post_gp(x[j], b[i], m[i], opt)
    } else if (j %% 3 == 2) {
      opt <- 1 
      Tab4[i,j] <- var_post_gp(x[j], b[i], m[i], opt)
    } else if (j %% 3 == 0) {
      opt <- 2 
      Tab4[i,j] <- var_post_gp(x[j], b[i], m[i], opt)
    }
  }
} # variance

Tab4 <- sqrt(Tab4) # standard deviation

for (i in 1:24) {
  for (j in 1:15) {
    if (j %% 3 == 1){
      Tab4[i,j] <- Tab4[i, j] 
    } else if (j %% 3 == 2){
      Tab4[i,j] <- Tab4[i,j] - Tab4[i,(j-1)]
    } else if (j %% 3 == 0){
      Tab4[i,j] <- Tab4[i,j] - Tab4[i,(j-2)]
    }
  }
}

Tab4 <- round(Tab4, 1) 
Tab4 # Table 4



#######################
### Posterior Mean - NB
#######################

# model 2 - NB

mean_post_nb <- function(y, b, m, opt=0){
  #opt = 0, exact
  #opt = 1, approximation 1, with Lerch ()
  #opt = 2, approximation 2
  if (opt == 0){
    mu <- sum(pospmf(y, b, m, model = 2, opt =0)*seq(1,10000))
    return(mu) # exact
  }
  else if (opt == 1){
    tau <- b*m/(1-m)
    a <- (1/(2*tau)+1)*y - 1/(2*tau)
    L1 <- lerch(m^tau, -y-1, a, tolerance = 1.0e-10, iter = 100)
    L2 <- lerch(m^tau, -y, a, tolerance = 1.0e-10, iter = 100)
    mu <- L1/L2 - (y-1)/(2*tau)
    return(mu) # approximation 1, Eq.(24)
  }
  else if (opt == 2){
    tau <- b*m/(1-m)
    mu_new <- -(y+1)/(tau*log(m)) - (y-1)/(2*tau)
    return(mu_new) # approximation 2, Theorem 3
  }
}


m <- c(rep(0.2, 4), rep(0.4, 4),rep(0.6, 4),rep(0.8, 4))
b <- c(rep(c(0.2, 0.4, 0.6, 0.8), 4))
y <- c(2,2,2,3,3,3,10,10,10,20,20,20,50,50,50)

#### Comparison of approximated means
Tab5 <- matrix(ncol = 15, nrow = 16)

for (i in 1:16) {
  for (j in 1:15) {
    if (j %% 3 == 1){
      opt <- 0
      Tab5[i,j] <- mean_post_nb(y[j], b[i], m[i], opt)
    } else if (j %% 3 == 2) {
      opt <- 1
      Tab5[i,j] <- mean_post_nb(y[j], b[i], m[i], opt) - Tab5[i,(j-1)]
    } else if (j %% 3 == 0) {
      opt <- 2
      Tab5[i,j] <- mean_post_nb(y[j], b[i], m[i], opt) - Tab5[i,(j-2)]
    }
  }
}

Tab5 <- round(Tab5, 1)
Tab5 # Table 5



###########################
### Posterior Variance - NB
###########################

# model 2 - NB

var_post_nb <- function(y, b, m, opt = 0){
  #opt = 0, exact
  #opt = 1, approximation 1, wtih Lerch ()
  #opt = 2, approximation 2
  if (opt == 0){
    mu <- sum(pospmf(y, b, m, model = 2, opt =0)*seq(1,10000))
    mu2 <- sum(pospmf(y, b, m, model = 2, opt =0)*(seq(1,10000))^2)
    var <- mu2 - mu^2 
    return(var) # exact
  } else if (opt == 1){
    tau <- b*m/(1-m)
    a <- (1/(2*tau)+1)*y - 1/(2*tau)
    L1 <- lerch(m^tau, -y-1, a, tolerance = 1.0e-10, iter = 100)
    L2 <- lerch(m^tau, -y, a, tolerance = 1.0e-10, iter = 100)
    L3 <- lerch(m^tau, -y-2, a, tolerance = 1.0e-10, iter = 100)
    var <- L3/L2 -(L1/L2)^2 
    return(var) # approximation 1, Eq.(25)
  } else if (opt == 2){
    tau <- b*m/(1-m)
    var_new <- (y+1)/(tau*log(m))^2 
    return(var_new) # approximation 2, Theorem 3
  }
}

m <- c(rep(0.2, 4), rep(0.4, 4),rep(0.6, 4),rep(0.8, 4))
b <- c(rep(c(0.2, 0.4, 0.6, 0.8), 4))
y <- c(2,2,2,3,3,3,10,10,10,20,20,20,50,50,50)

Tab6 <- matrix(ncol = 15, nrow = 16)
for (i in 1:16) {
  for (j in 1:15) {
    if (j %% 3 == 1){
      opt <- 0
      Tab6[i,j] <- var_post_nb(y[j], b[i], m[i], opt)
    } else if (j %% 3 == 2) {
      opt <- 1
      Tab6[i,j] <- var_post_nb(y[j], b[i], m[i], opt)
    } else if (j %% 3 == 0) {
      opt <- 2
      Tab6[i,j] <- var_post_nb(y[j], b[i], m[i], opt)
    }
  }
} # variance

Tab6 <- sqrt(Tab6) # standard deviation

for (i in 1:16) {
  for (j in 1:15) {
    if (j %% 3 == 1){
      Tab6[i,j] <- Tab6[i, j] 
    } else if (j %% 3 == 2){
      Tab6[i,j] <- Tab6[i,j] - Tab6[i,(j-1)]
    } else if (j %% 3 ==0){
      Tab6[i,j] <- Tab6[i,j] - Tab6[i,(j-2)]
    }
  }
}

Tab6 <- round(Tab6, 1) 
Tab6 # Table 6








#Usage example
#
#excellent fit
#prob_exact <- pospmf(observed_count=50, b=0.1, m=1.2, opt=0)
#prob_approx <- pospmf(observed_count=50, b=0.1, m=1.2, opt=1)
#
#h <- which(prob_exact == max(prob_exact))

#plot(prob_exact[1:(2*h)], type="l", xlab="k", ylab="Probability")
#points(prob_approx[1:(2*h)], type="l", col="red")
#legend(2*h-sqrt(2*h*100),max(prob_exact), lty=rep(1,2), 
#col=c("black","red"), c("Exact","Approximation"), cex=0.7)

#poorer fit
#prob_exact <- pospmf(observed_count=5, b=0.5, m=2, opt=0)
#prob_approx <- pospmf(observed_count=5, b=0.5, m=2, opt=1)

#h <- which(prob_exact == max(prob_exact))

#plot(prob_exact[1:(10*h)], type="l", xlab="k", ylab="Probability")
#points(prob_approx[1:(10*h)], type="l", col="red")
#legend(10*h-sqrt(10*h*2.5),max(prob_exact), lty=rep(1,2), 
#col=c("black","red"), c("Exact","Approximation"), cex=0.4)

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
