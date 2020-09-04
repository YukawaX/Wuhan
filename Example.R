#-------load packages--------
rm(list = ls())
library(tidyverse)
library(lubridate)
library(pinyin)
library(smooth)
library(Mcomp)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(ggpubr)
library(reshape2)
#----------Data Smoothing-----------
## Method 1: Boundary Kernel + Local Constant
{
  get_Boundary_Weight <- function(t, N, h = 5.5){
    
    get_muk <- function(l,p){
      f <- function(u) 0.75*u^l*(1-u^2)
      muk = integrate(f,p,1)
      return(muk$value)
    }
    
    get_K <- function(t){
      if(0.75*(1-t^2)>0){
        return(0.75*(1-t^2))
      }else{return(0)}
    }
    
    get_B <- function(pt,t){
      K = get_K(t)
      if(pt <= -1){
        return(K)
      }else if(pt <= 0){
        muk0 = get_muk(0,pt)
        muk1 = get_muk(1,pt)
        muk2 = get_muk(2,pt)
        B = (muk2-muk1*t) / (muk0*muk2-muk1^2)*K
        return(B)
      }else{return(NA)}
    }
    
    Weight = c()
    if(t > (N+1)/2){
      pt = (t-N) / h
      for(i in 1:N){w = get_B(pt,(t-i)/h);Weight = c(Weight,w)}
    }else{
      pt = (1-t) / h
      for(i in 1:N){ w = get_B(pt,-(t-i)/h); Weight = c(Weight,w)}
    }
    
    return(Weight)
  }
  
  ## One-point Smooth
  Weighted.Ave <- function(x,h){
    N = length(x)
    y = c()
    for(i in 1:N){
      y[i] = weighted.mean(x,get_Boundary_Weight(i,N,h))
    }
    return(y)
  }
}


## Smooth function ; 
mav <- function(x,h){
  
  dx=diff(x)
  index = 1:length(dx)
  dy = Weighted.Ave(dx,h)
  # dy = Local.Linear(dx,h)
  y = cumsum(c(x[1],dy))
  y = as.double(y)
  
  return(y)
}


#-------Simulation function--------
#' simulate the trajectory under SEIdR model
#'
#' @param t_max range of the time interval 
#' @param y a vector of length 4 consists of initial Exposed, Infected, Death, Recovered
#' @param alpha time-varying diagnosis rate alpha (vector)
#' @param beta time-varying beta
#' @param M total population
#' @param r.beta beta.E / beta.I, set to be 5 by default
sim.vseidr = function(t_max,y, alpha, beta, gamma.d,gamma.r,M,r.beta = 5) {
  t = 0 : t_max
  S = rep(0, t_max+1)
  E = rep(0, t_max+1)
  I = rep(0, t_max+1)
  R.d = rep(0, t_max+1)
  R.r = rep(0, t_max+1)
  #initial value 
  E[1] = y[1]
  I[1] = y[2]
  R.d[1] = y[3]
  R.r[1] = y[4]
  S[1] = M - E[1] - I[1] - R.d[1]-R.r[1]
  
  for (i in 1:(t_max)) {
    dS = S[i] / M * (beta[i] * E[i] + beta[i] / r.beta * I[i]) 
    dE = alpha[i] * E[i]
    dR.d = gamma.d[i] * I[i]
    dR.r = gamma.r[i] * I[i]
    
    dS = min(rpois(1, dS),S[i])
    dE = min(rpois(1, dE),E[i] + dS)
    dR.d = min(rpois(1, dR.d), I[i] + dE)
    dR.r = min(rpois(1, dR.r), I[i] + dE-dR.d)
    
    S[i+1] = S[i] - dS  
    E[i+1] = E[i] + dS - dE 
    I[i+1] = I[i] + dE - dR.d- dR.r 
    R.d[i+1] = R.d[i] + dR.d
    R.r[i+1] = R.r[i] + dR.r
  }
  return(data.frame(t = t, S = S, E = E, I = I, R = R.d+R.r, N = I+R.d+R.r,R.d = R.d, R.r = R.r))
}


#----------Estimation----------
#' regE method for parameter estimation in SEIdR model 
#'
#' @param w bandwidth
#' @param D duration of infection, set to be 14 by default
vSEIDR_Coef <- function(A4, w,D, alpha, M,r = 5){
  gamma.temp = rep(NA,nrow(A4))
  gamma.d.temp = rep(NA,nrow(A4))
  gamma.r.temp = rep(NA,nrow(A4))
  beta.temp = rep(NA,nrow(A4))
  A4$alpha = alpha
  A4$E = c(diff(A4$N) / alpha[-length(alpha)], NA)
  A4$S = (M-A4$N-A4$E) / M
  
  A4$R.diff2 = c(diff(A4$R),NA)
  A4$R.d.diff2 = c(diff(A4$R.d),NA)
  A4$R.r.diff2 = c(diff(A4$R.r),NA)
  A4$y.temp=c(diff(A4$E),NA)+(alpha*A4$E)
  A4$x.temp=A4$I/r + A4$E
  A4$x.c=A4$x.temp; A4$x.c[nrow(A4)-1]=NA;
  
  for (j in 1:(nrow(A4)-2)){
    W = c(get_Boundary_Weight(j, (nrow(A4)-2), w/2),0,0)
    gamma.temp[j+1] = sum(A4$R.diff2 *A4$I*W,na.rm=T)/sum(W*(A4$I)^2,na.rm=T)
    gamma.d.temp[j+1] = sum(A4$R.d.diff2 *A4$I*W,na.rm=T)/sum(W*(A4$I)^2,na.rm=T)
    gamma.r.temp[j+1] = sum(A4$R.r.diff2 *A4$I*W,na.rm=T)/sum(W*(A4$I)^2,na.rm=T)
    beta.temp[j+1] = sum(A4$y.temp*A4$x.c*W,na.rm=T)/sum(W*(A4$x.c)^2,na.rm=T)
  }
  coef.one = data.frame(gamma = gamma.temp,
                        gamma.d = gamma.d.temp,
                        gamma.r = gamma.r.temp,
                        beta = beta.temp,
                        beta.hat = beta.temp / A4$S,
                        Rt = beta.temp * (1/alpha + D/r) / A4$S)
  
  coef.one = cbind(A4,coef.one)
  coef.one = coef.one %>% dplyr::select(
                                        'N','I','R','R.d','R.r','E','S',
                                        'alpha','beta','beta.hat',
                                        'gamma','gamma.d','gamma.r',
                                        'Rt')
  coef.one$Rt[which(coef.one$Rt<0)]=0
  coef.one$beta.hat[which(coef.one$beta.hat<0)] = 0.001
  coef.one$R.r[which(coef.one$R.r<0)] = 0
  coef.one$R.d[which(coef.one$R.d<0)] = 0
  return(coef.one)
}
#------- a example---------
t_max = 150
alpha = rep(0.25,t_max+1)
beta = rep(0.3,t_max+1)
M = 50000000
init = c(50,1,1,0)
gamma.d = rep(1/14,t_max+1)
gamma.r = rep(1/14,t_max+1)
simulated <- sim.vseidr(t_max = t_max,y = init,alpha = alpha,beta = beta,M = M,gamma.d = gamma.d, gamma.r = gamma.r)
result <- vSEIDR_Coef(simulated,w = 7, D = 14, alpha = alpha,M = M,r = 5) #alpha is assumed to be known

