#-------- library packages------
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
library(ggsci)
library(doMC)
registerDoMC(4)
library(xtable)
R0 = 5.7
#------ Smoothing -------
## Method Boundary Kernel + Local Constant
{
  get_Boundary_Weight <- function(t, N, h){
    
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


#--------------some functions-----------
##estimation under vSEIdRm model
vSEIDR_Coef_new <- function(A4, w,D, alpha, M,pt,r = 5){
  A4 = left_join(A4,pt,by = "date")
  A4$pt[is.na(A4$pt)] = 0
  pt = A4$pt
  pE = pt #proportion of E 
  gamma.temp = rep(NA,nrow(A4))
  gamma.d.temp = rep(NA,nrow(A4))
  gamma.r.temp = rep(NA,nrow(A4))
  beta.temp = rep(NA,nrow(A4))
  A = A4$daily_sum
  A4$alpha = alpha
  A4$E = c(diff(A4$N) / alpha[-length(alpha)], NA)
  A[is.na(A)] = 0
  M = M + cumsum(A)
  A4$S = (M-A4$N-A4$E) / M
  
  A4$R.diff2 = c(diff(A4$R),NA)
  A4$R.d.diff2 = c(diff(A4$R.d),NA)
  A4$R.r.diff2 = c(diff(A4$R.r),NA)
  A4$y.temp=c(diff(A4$E),NA)+(alpha*A4$E) - pE * A
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
                        Rt = beta.temp * (1/alpha + D/r) / A4$S,
                        A = A
  )
  
  coef.one = cbind(A4,coef.one)
  coef.one$date = as.character(coef.one$date)
  coef.one = coef.one %>% dplyr::select('province','date','infected','dead','recovered',
                                        'N','I','R','R.d','R.r','E','S',
                                        'alpha','beta','beta.hat',
                                        'gamma','gamma.d','gamma.r',
                                        'Rt','pt','A')
  coef.one$Rt[which(coef.one$Rt<0)]=0
  coef.one$beta.hat[which(coef.one$beta.hat<0)] = 0.001
  coef.one$R.r[which(coef.one$R.r<0)] = 0
  coef.one$R.d[which(coef.one$R.d<0)] = 0
  return(coef.one)
}
## Estimation under vSEIdR model
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
  coef.one$date = as.character(coef.one$date)
  coef.one = coef.one %>% dplyr::select('province','date','infected','dead','recovered',
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
## poisson simulation function
sim.vseidr = function(t_max, pt,y,A, alpha, beta, gamma.d,gamma.r,M,r.beta = 5) {
  pt = pt$pt
  t = 0 : t_max
  S = rep(0, t_max+1)
  E = rep(0, t_max+1)
  I = rep(0, t_max+1)
  R.d = rep(0, t_max+1)
  R.r = rep(0, t_max+1)
  pE = pt
  M = M + cumsum(A)
  #initial value ²»¶¯ 
  E[1] = y[1]
  I[1] = y[2]
  R.d[1] = y[3]
  R.r[1] = y[4]
  S[1] = M[1] - E[1] - I[1] - R.d[1]-R.r[1]
  
  for (i in 1:(t_max)) {
    dS = S[i] / M[i] * (beta[i] * E[i] + beta[i] / r.beta * I[i]) 
    dE = alpha[i] * E[i]
    dR.d = gamma.d[i] * I[i]
    dR.r = gamma.r[i] * I[i]
    
    dS = min(rpois(1, dS),S[i])
    dE = min(rpois(1, dE),E[i] + dS)
    dR.d = min(rpois(1, dR.d), I[i] + dE)
    dR.r = min(rpois(1, dR.r), I[i] + dE-dR.d)
    
    S[i+1] = S[i] - dS + (1 - pE[i]) * A[i] 
    E[i+1] = E[i] + dS - dE + pE[i] * A[i]
    I[i+1] = I[i] + dE - dR.d- dR.r 
    R.d[i+1] = R.d[i] + dR.d
    R.r[i+1] = R.r[i] + dR.r
  }
  return(data.frame(t = t, S = S, E = E, I = I, R = R.d+R.r, N = I+R.d+R.r,R.d = R.d, R.r = R.r))
}
## deterministic simulation function
sim.vseidr.det = function(t_max, pt,y,A, alpha, beta, gamma.d,gamma.r,M,r.beta = 5) {
  t = 0 : t_max
  S = rep(0, t_max+1)
  E = rep(0, t_max+1)
  I = rep(0, t_max+1)
  R.d = rep(0, t_max+1)
  R.r = rep(0, t_max+1)
  pE = pt
  M = M + cumsum(A)
  E[1] = y[1]
  I[1] = y[2]
  R.d[1] = y[3]
  R.r[1] = y[4]
  S[1] = M[1] - E[1] - I[1] - R.d[1]-R.r[1]
  
  for (i in 1:(t_max)) {
    dS = S[i] / M[i] * (beta[i] * E[i] + beta[i] / r.beta * I[i]) 
    dE = alpha[i] * E[i]
    dR.d = gamma.d[i] * I[i]
    dR.r = gamma.r[i] * I[i]
    S[i+1] = S[i] - dS + (1 - pE[i]) * A[i] 
    E[i+1] = E[i] + dS - dE + pE[i] * A[i]
    I[i+1] = I[i] + dE - dR.d- dR.r 
    R.d[i+1] = R.d[i] + dR.d
    R.r[i+1] = R.r[i] + dR.r
  }
  return(data.frame(t = t, S = S, E = E, I = I, R = R.d+R.r, N = I+R.d+R.r,R.d = R.d, R.r = R.r))
}
## Get data; start at `st`
vSEIDR_Data <- function(opendir, h, st,end){
  A = read.csv(opendir, head = TRUE, fileEncoding = 'UTF-8',stringsAsFactors = F)
  if(A$city[1] == "total"){
    A1 = A[A$infected != 0, ] %>% dplyr::select(province, date, infected, dead, recovered)
  }else{
    A1 = A[A$infected != 0, ] %>% dplyr::select(city, date, infected, dead, recovered)
    A1$date = c(sub(".csv", "", A1$date))
    colnames(A1)[1] <- "province"
  }
  if(unique(A1$province) %in% c("Korea","Italy","UK")){
    A1 = A1 %>% mutate(date = ymd(date))
  }
  else{
    A1 = A1 %>% mutate(date = mdy(date))
  }
  A1 = A1 %>% mutate(R.d = dead) %>% mutate(R.r = recovered)
  A1$N=A1$infected
  
  A1 = A1[which(A1$date>=st&A1$date<=end),]
  
  A1$N=mav(A1$N, h)
  A1$R.d = mav(A1$R.d,h); A1$R.r = mav(A1$R.r,h)
  
  A1$R = A1$R.d + A1$R.r
  A1$I = A1$N - A1$R
  A1$date=as.Date(A1$date)
  
  return(A1)
}
#------------- Minimum distance estimation of parameter pt------------
path0 = "/Users/gujia/Documents/National Bureau of Statistics/Data/nCoV_Pinyin0418/"
p.t = c(mdy("01102020"),mdy("01232020"))
pE <- function(par,path0,province,p.t,A,M){
  p.E = par[1]
  p.2 = par[2]
  alpha0 = par[3]
  alpha1 = par[4]
  opendir = paste0(path0,province,".csv")
  A4 = vSEIDR_Data(opendir = opendir, h = 5.5, st = mdy("01222020"),end = mdy("03152020"))
  R0 = 5.7 
  p.seq = p.E / 100 * exp(R0 / 14 * 1:14)
  A.seq = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% dplyr::select(daily_sum)
  E0 = sum(A.seq * p.seq)
  A4 = left_join(A4,A,by = c("date","province"))
  t = p.t[2] + (1:28)
  pt = data.frame(date = t,pt = rep(p.2,28))
  alpha = c(alpha0,cumsum(rep((alpha1 - alpha0)/14,14)) +alpha0)  #linearly increasing
  alpha = c(alpha,rep(alpha1,dim(A4)[1] - length(alpha)))
  
  Coef = vSEIDR_Coef_new(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,pt = pt,r = 5)
  Coef.sim <- Coef %>% filter(date<=ymd("20200219"),date >= ymd("20200123"))
  sim = sim.vseidr.det(t_max = dim(Coef.sim)[1],pt = Coef.sim$pt,y = c(E0,Coef.sim$I[1],Coef.sim$R.d[1],Coef.sim$R.r[1]),
                       A = c(0,Coef.sim$A[-1]),alpha = Coef.sim$alpha,beta = Coef.sim$beta.hat, gamma.d  = Coef.sim$gamma.d,gamma.r = Coef.sim$gamma.r,M = M)
  N.part = (Coef %>% filter(date<=ymd("20200220"),date >= ymd("20200124")))$N
  mean(sqrt((sim$N[-1] - N.part)^2))
}

fitting <- function(par,path0,province,p.t,A,M){
  p.E = par[1]
  p.2 = par[2]
  alpha0 = par[3]
  alpha1 = par[4]
  opendir = paste0(path0,province,".csv")
  A4 = vSEIDR_Data(opendir = opendir, h = 5.5, st = mdy("01222020"),end = mdy("03152020"))
  R0 = 5.7
  p.seq = p.E / 100 * exp(R0 / 14 * 1:14)
  A.seq = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% dplyr::select(daily_sum)
  E0 = sum(A.seq * p.seq)
  A4 = left_join(A4,A,by = c("date","province"))
  t = p.t[2] + (1:28)
  pt = data.frame(date = t,pt = rep(p.2,28))
  alpha = c(alpha0,cumsum(rep((alpha1 - alpha0)/14,14)) +alpha0) #linearly increasing
  alpha = c(alpha,rep(alpha1,dim(A4)[1] - length(alpha)))
  Coef = vSEIDR_Coef_new(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,pt = pt,r = 5)
  Coef.origin = vSEIDR_Coef(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,r = 5) # not consider migration, i.e. original SEIdR model
  Coef.sim <- Coef %>% filter(date<=ymd("20200219"),date >= ymd("20200123"))
  sim = sim.vseidr.det(t_max = dim(Coef.sim)[1],pt = Coef.sim$pt,y = c(E0,Coef.sim$I[1],Coef.sim$R.d[1],Coef.sim$R.r[1]),
                       A = c(0,Coef.sim$A[-1]),alpha = Coef.sim$alpha,beta = Coef.sim$beta.hat, gamma.d  = Coef.sim$gamma.d,gamma.r = Coef.sim$gamma.r,M = M)
  N.part = (Coef %>% filter(date<=ymd("20200220"),date >= ymd("20200124")))$N
  data1 = data.frame(date = ymd("20200124")+(0:27),N = sim$N[-1],province = province,type = "fitting")
  data2 = data.frame(date = ymd("20200124")+(0:27),N = N.part,province = province,type = "observed")
  list(fitting = rbind(data1,data2),Coef = Coef, Coef.origin = Coef.origin)
}

#------------- Estimation of pt------------
provinces = c("Anhui","Beijing","Chongqin","Fujian","Guangdong","Guangxi","Hebei","Heilongjiang","Henan","Hunan",
              "Jiangsu","Jiangxi","Shaanxi","Shandong","Shanghai","Sichuan","Zhejiang")
population.all = read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/country_population.csv")

summary.pt <- foreach(i = 1:17, .combine = rbind) %dopar% {
  ui.1 = matrix(c(0,exp(R0) / 100,0,0,0),ncol = 1)#constraint matrix
  ui.2 = matrix(c(1,-1,0,0,0),ncol = 1)
  ui.3 = matrix(c(0,0,1,0,-1),ncol = 1)
  ui.4 = matrix(c(0,0,0,-1,1),ncol = 1)
  ui = cbind(ui.1,ui.2,ui.3,ui.4)
  ci = c(0,0,1/56,-1/2,0)
  province.temp = provinces[i]
  population = population.all$population[population.all$city == province.temp] * 10000
  A = read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/year2020_Wuhan.csv")%>%mutate(province = as.character(province)) %>% filter(province == province.temp) %>% mutate(date = ymd(date))
  if(province.temp %in% c("Anhui","Chongqin","Fujian","Hunan","Jiangxi","Shaanxi")){
    op = constrOptim(theta = c(1/1000,1/1000,1/11.5,1/7),pE,grad = NULL,ui = ui,ci = ci,path0 = path0,province = province.temp,p.t = p.t, A = A,M = population)
    print(op$counts)
    print(op$value)
    print(i)
  }
  else{
    op = constrOptim(theta = c(1/50,1/50,1/11.5,1/7),pE,grad = NULL,ui = ui,ci = ci,path0 = path0,province = province.temp,p.t = p.t, A = A,M = population)
    print(op$counts)
    print(op$value)
    print(i)
  }
  
  N = (data$infected)[data$date == ymd("20200315")]
  Gof = op$value
  percentage = op$par[1]
  percentage.2 = op$par[2]
  alpha0 = op$par[3]
  alpha1 = op$par[4]
  R0 = 5.7
  p.seq = percentage / 100 * exp(R0 / 14 * 1:14)
  A.seq = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% dplyr::select(daily_sum)
  Import = sum(floor(A.seq * p.seq))
  A.seq.1 = A %>% filter(date > p.t[2],date <= p.t[2] + 28) %>% dplyr::select(daily_sum)
  p.seq.1 = rep(percentage.2,28)
  Import.afterlockdown = sum(floor(A.seq.1 * p.seq.1))
  
  summary.temp <- data.frame(province = province.temp,Gof = Gof,E = Import,N = N,proportion = Import / N,Import.cont = Import.afterlockdown,eta = percentage, pt.max = max(p.seq),pt.2 = percentage.2,
                             alpha0 = alpha0,alpha1 = alpha1)
  summary.temp
}
write.csv(summary.pt,"/Users/gujia/Documents/National Bureau of Statistics/Output/pt.csv")
  

eta = summary.pt$eta
pt.2 = summary.pt$pt.2
alpha0 = summary.pt$alpha0
alpha1 = summary.pt$alpha1
fitting.performance = data.frame(NULL)
Coef = data.frame(NULL)
Coef.origin = data.frame(NULL) #estimates without considering the effect of imported caseseat 
for(i in 1:17){
  par.temp <- c(eta[i],pt.2[i],alpha0[i],alpha1[i])
  province.temp <- as.character(summary.pt$province[i])
  population = population.all$population[population.all$city == province.temp] * 10000
  A = read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/year2020_Wuhan.csv")%>%mutate(province = as.character(province)) %>% filter(province == province.temp) %>% mutate(date = ymd(date))
  fitting.results = fitting(par.temp,path0 = path0,province = province.temp,p.t = p.t, A = A, M = population)
  fitting.performance = rbind(fitting.performance,fitting.results$fitting)
  Coef = rbind(Coef,fitting.results$Coef)
  Coef.origin = rbind(Coef.origin,fitting.results$Coef.origin)
}




##plot of fitting performance
fitting.plot <- ggplot(fitting.performance,aes(date,N,color = type))+facet_wrap(~province,scales = "free")+geom_line()+scale_color_npg()+
  theme_minimal()+theme(
    axis.title = element_text(size = 8), 
    strip.text = element_text(size = 10,face = 'bold'),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=8, angle=45,hjust = 1, face = 'bold'),
    axis.text.y = element_text(size=10, face = 'bold'),
    plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size=10, face = 'bold'),
    legend.key.width  = unit(.3,"inches"),
    legend.key.height = unit(.3,"inches"))
fitting.plot
ggsave("/Users/gujia/Documents/National Bureau of Statistics/Output/fitting.png")


## Comparative plots of estimation of Rt under two models
save.dir = "/Users/gujia/Documents/National Bureau of Statistics/Output/"
write.csv(Coef,paste0(save.dir,"Results0730.csv"))
Coef<-Coef %>% mutate(date = ymd(date)) %>% filter(date <= ymd("20200221"))
Coef.origin <- Coef.origin %>% mutate(date = ymd(date)) %>% filter(date <= ymd("20200221"))
Coef.1 <- Coef[,1:19]
Coef.1$type = "vSEIdRm"
Coef.2 <- Coef.origin
Coef.2$type = "vSEIdR"
Coef.combine = rbind(Coef.1,Coef.2)
Migration.Wuhan2020 <- Coef %>% mutate(import.cases = floor(pt * A)) %>% filter(date > ymd("20200123"),date <= ymd("20200221"))



ggplot(Coef.combine,aes(date,Rt,color = type))+facet_wrap(~province,scales = "free")+scale_color_npg()+geom_line()+theme_minimal()+
  geom_point(data = Migration.Wuhan2020,aes(date,import.cases, color = "Estimated Imported Cases"))+scale_y_continuous(sec.axis = sec_axis(~ . * 1))+
  theme(axis.title = element_text(size = 8),
        strip.text = element_text(size = 10,face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, angle=45,hjust = 1, face = 'bold'),
        axis.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=10, face = 'bold'),
        legend.key.width  = unit(.3,"inches"),
        legend.key.height = unit(.3,"inches"),
        legend.position = c(0.75,0.1))
ggsave(paste0(save.dir,"Rt_Compare_Import.png"))













#Spatial Temporal Model
covariates <- read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/Covariates.csv")
covariates$province = as.character(covariates$province)
covariates <- covariates %>% mutate(GDP.per = GDP * 1e4 / population)
covariates.part <- covariates %>% filter(province %in% provinces)
covariates.part[,c(4:7,9)] <- scale(covariates.part[,c(4:7,9)]) / sqrt(nrow(covariates.part) - 1)
population = read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/country_population.csv")
migration.2020 <- read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/year2020_Wuhan.csv")
#Infection data
Infect.data <- data.frame(NULL)
for(i in 1:length(covariates.part$province)){
  province.selected.i = covariates.part$province[i]
  path0 = "/Users/gujia/Documents/National Bureau of Statistics/Data/nCoV_Pinyin0418/"
  opendir = paste0(path0,province.selected.i,".csv")
  A4 = vSEIDR_Data(opendir = opendir, h = 5.5, st = mdy("01222020"),end = mdy("03152020"))
  Infect.data = rbind(Infect.data,A4)
}


Migration.temporal <- function(theta,covariates,population,migration,infect){
  eta <- covariates %>% mutate(eta = exp(flow * theta[1]  + GDP.per * theta[2] + population * theta[3] + distance * theta[4]) ) %>% select(province,eta)
  alpha0 = 1 / 56
  alpha1 = 1 / 3
  province.selected <- covariates$province
  covariates$uniparam <- theta[5:(4 + length(covariates$province))]
  obj <- foreach(i = 1:length(province.selected),.combine = sum) %dopar% {
    province.selected.i = province.selected[i]
    uniparam.i <- covariates$uniparam[covariates$province == province.selected.i]
    M = population$population[population$city == province.selected.i] * 10000
    A = migration %>% mutate(province = as.character(province)) %>% filter(province == province.selected.i) %>% mutate(date = ymd(date))
    A4 = infect %>% filter(province == province.selected.i)
    R0 = 5.7 
    p.E = eta[i,2] * exp(uniparam.i)
    p.seq = p.E / (1e4) * exp(R0 / 14 * 1:14)
    A.seq = A %>% filter(date >= ymd("20200110"),date <= ymd("20200123")) %>% dplyr::select(daily_sum)
    E0 = sum(A.seq * p.seq)
    A4 = left_join(A4,A,by = c("date","province"))
    t = ymd("20200123") + (1:28)
    pt = data.frame(date = t,pt = rep(0,28))
    alpha = c(alpha0,cumsum(rep((alpha1 - alpha0)/14,14)) +alpha0)  #linearly increasing
    alpha = c(alpha,rep(alpha1,dim(A4)[1] - length(alpha)))
    
    Coef = vSEIDR_Coef_new(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,pt = pt,r = 5)
    Coef.sim <- Coef %>% filter(date<=ymd("20200219"),date >= ymd("20200123"))
    sim = sim.vseidr.det(t_max = dim(Coef.sim)[1],pt = Coef.sim$pt,y = c(E0,Coef.sim$I[1],Coef.sim$R.d[1],Coef.sim$R.r[1]),
                         A = c(0,Coef.sim$A[-1]),alpha = Coef.sim$alpha,beta = Coef.sim$beta.hat, gamma.d  = Coef.sim$gamma.d,gamma.r = Coef.sim$gamma.r,M = M)
    N.part = (Coef %>% filter(date<=ymd("20200220"),date >= ymd("20200124")))$N
    objective <- objective + mean(sqrt((sim$N[-1] - N.part)^2))
    objective
  }
  print(obj / length(province.selected))
  obj 
}
op = optim(par = op$par,fn = Migration.temporal,covariates = covariates.part,population = population,migration = migration.2020,infect = Infect.data,control = list(maxit =5000))
write.csv(op$par,"/Users/gujia/Documents/National Bureau of Statistics/Output/Spatial/op.csv")
fitting.temporal <- function(theta,covariates,population,migration,infect){
  eta <- covariates %>% mutate(eta = exp(flow * theta[1]  + GDP.per * theta[2] + population * theta[3] + distance * theta[4]) )%>% select(province,eta)
  alpha0 = 1 / 56
  alpha1 = 1 / 3
  province.selected <- covariates$province
  covariates$uniparam <- theta[5:(4 + length(covariates$province))]
  dat.fit <- data.frame(NULL)
  E <- c()
  pE.max <- c()
  eta.seq <- c()
  D.seq <- c()
  Coef.all <- data.frame(NULL)
  obj <- for(i in 1:length(province.selected)){
    province.selected.i = province.selected[i]
    M = population$population[population$city == province.selected.i] * 10000
    A = migration %>% mutate(province = as.character(province)) %>% filter(province == province.selected.i) %>% mutate(date = ymd(date))
    A4 = infect %>% filter(province == province.selected.i)
    R0 = 5.7 
    uniparam.i <- covariates$uniparam[covariates$province == province.selected.i]
    p.E = eta[i,2] * exp(uniparam.i)
    eta.seq <- c(eta.seq,p.E / 100)
    p.seq = p.E / (1e4) * exp(R0 / 14 * 1:14)
    pE.max <- c(pE.max,max(p.seq))
    A.seq = A %>% filter(date >= ymd("20200110"),date <= ymd("20200123")) %>% dplyr::select(daily_sum)
    E0 = sum(A.seq * p.seq)
    E <- c(E,E0)
    A4 = left_join(A4,A,by = c("date","province"))
    t = ymd("20200123") + (1:28)
    pt = data.frame(date = t,pt = rep(0,28))
    alpha = c(alpha0,cumsum(rep((alpha1 - alpha0)/14,14)) +alpha0)  #linearly increasing
    alpha = c(alpha,rep(alpha1,dim(A4)[1] - length(alpha)))
    
    Coef = vSEIDR_Coef_new(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,pt = pt,r = 5)
    Coef.sim <- Coef %>% filter(date<=ymd("20200219"),date >= ymd("20200123"))
    Coef.all <- rbind(Coef.all,Coef)
    sim = sim.vseidr.det(t_max = dim(Coef.sim)[1],pt = Coef.sim$pt,y = c(E0,Coef.sim$I[1],Coef.sim$R.d[1],Coef.sim$R.r[1]),
                         A = c(0,Coef.sim$A[-1]),alpha = Coef.sim$alpha,beta = Coef.sim$beta.hat, gamma.d  = Coef.sim$gamma.d,gamma.r = Coef.sim$gamma.r,M = M)
    N.part = (Coef %>% filter(date<=ymd("20200220"),date >= ymd("20200124")))$N
    D.seq <- c(D.seq,mean(sqrt((sim$N[-1] - N.part)^2)))
    data1 = data.frame(date = ymd("20200124")+(0:27),N = sim$N[-1],province = province.selected.i,type = "fitting")
    data2 = data.frame(date = ymd("20200124")+(0:27),N = N.part,province = province.selected.i,type = "observed")
    dat <- rbind(data1,data2)
    dat.fit <- rbind(dat.fit,dat)
  }
  Import <- data.frame(province = province.selected,E = floor(E),N = floor(covariates$N), proportion = E / covariates$N ,eta = eta.seq, pE.max = pE.max,D = D.seq)
  fitting.plot <- ggplot(dat.fit,aes(date,N,color = type))+facet_wrap(~province,scales = "free")+geom_line()+scale_color_npg()+
    theme_minimal()+theme(
      axis.title = element_text(size = 8), 
      strip.text = element_text(size = 10,face = 'bold'),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=8, angle=45,hjust = 1, face = 'bold'),
      axis.text.y = element_text(size=10, face = 'bold'),
      plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size=10, face = 'bold'),
      legend.key.width  = unit(.3,"inches"),
      legend.key.height = unit(.3,"inches"),
      legend.position = c(0.7,0.05))
  Coef.all <- Coef.all %>% mutate(date = ymd(date)) 
  list(import = Import,fitting = fitting.plot,Coef = Coef.all)
}
temporal.result <- fitting.temporal(theta = op$par,covariates = covariates.part,population = population,migration = migration.2020,infect = Infect.data)
ggsave("/Users/gujia/Documents/National Bureau of Statistics/Output/Spatial/fitting_spatial.png")
##Check
##ggplot(temporal.result$Coef,aes(date,Rt)) + facet_wrap(~province,scales = "free") + geom_line()
write.csv(temporal.result$import,"/Users/gujia/Documents/National Bureau of Statistics/Output/pt_temporal.csv")
write.csv(temporal.result$Coef,"/Users/gujia/Documents/National Bureau of Statistics/Output/Results0902")






#----------- Scenario Analysis--------------
provinces = data.frame(provinces = c("Anhui","Beijing","Chongqin","Fujian","Guangdong","Guangxi","Hebei","Heilongjiang","Henan","Hunan",
                                     "Jiangsu","Jiangxi","Shaanxi","Shandong","Shanghai","Sichuan","Zhejiang")) %>% mutate(provinces = as.character(provinces))
population = read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/country_population.csv") %>% mutate(city = as.character(city))
colnames(population) = c("provinces","population")
provinces = left_join(provinces,population,by = "provinces")
start = read.xlsx("/Users/gujia/Documents/National Bureau of Statistics/Data/start.xlsx")
start <- start %>% mutate(st = ymd(st),end = ymd(end))  
migration = read.csv("/Users/gujia/Documents/National Bureau of Statistics/Data/year2019_Wuhan.csv")
opendir = "/Users/gujia/Documents/National Bureau of Statistics/Output/Results0730.csv"

scenario.domestic.true <- function(provinces,start,opendir,migration,end.date = ymd("2020-02-21"),w = 7){
  discrete.integral <-function(beta){
    size = c()
    for(i in 1:(length(beta) - 1 )){
      size = c(size,beta[i] + beta[i+1])
    }
    sum(size) / 2
  }
  set.seed(123456)
  data.summary  = data.frame(NULL)
  beta.summary = data.frame(NULL)
  results = read.csv(opendir)
  results$province = as.character(results$province)
  pt.full <-read.csv("/Users/gujia/Documents/National Bureau of Statistics/Output/pt.csv")
  import = c()
  beta.all = c()
  for(i in 1:dim(provinces)[1]){
    print(i)
    province.temp = provinces$provinces[i]
    population.temp = provinces$population[i]
    population.temp = provinces$population[i] * 10000
    results.temp <- results %>% mutate(date = ymd(date)) %>% filter(province ==province.temp,date >= ymd("2020-01-23"),date < end.date) 
    alpha.temp = results.temp$alpha
    migration.temp <- migration %>% mutate(date = ymd(date)) %>% filter(province == province.temp,date >= ymd(results.temp$date[1]),date <= end.date)
    beta.temp = data.frame(NULL)
    A = migration.temp$daily_sum
    A[1] = 0
    beta = results.temp$beta.hat
    gamma.r = results.temp$gamma.r
    gamma.d = results.temp$gamma.d
    pt = results.temp %>% dplyr::select(date,pt)
    T0 = dim(pt)[1]
    pt.max = pt.full[pt.full$province == province.temp,]$pt
    pt.new = pt$pt
    pt.new = c(rep(pt.max,14),pt.max - cumsum(rep(pt.max / 14,14)),rep(0,T0 - 28))
    import = c(import,sum(pt.new * A))
    pt$pt = pt.new
    E0 = pt.full$E[pt.full$province == province.temp]
    I0 = results.temp$I[1]
    R0.d = results.temp$R.d[1]
    R0.r = results.temp$R.r[1]
    if(R0.d<=0){
      R0.d = 1
    }
    if(R0.r <= 0){
      R0.r = 1
    }
    gamma.r[gamma.r<0] = 0.0001
    gamma.d[gamma.d<0] = 0.0001
    y = c(E0,I0,R0.d,R0.r)
    T.max = nrow(migration.temp)
    beta.temp <- foreach(icount(50),.combine=rbind) %dopar%{
      vseidr = sim.vseidr(T.max - 1,y = y,A = A,alpha = alpha.temp,beta = beta,gamma.d = gamma.d,gamma.r = gamma.r,pt = pt,M = population.temp)
      vseidr$date = migration.temp$date
      vseidr.combine = left_join(vseidr,migration.temp,by = "date")
      beta.est = vSEIDR_Coef_new(vseidr.combine,w = w, D = 14, alpha = alpha.temp,M = population.temp,pt = pt)$beta.hat
      beta.est
    }
    beta.mean = colMeans(beta.temp)
    beta.corrected = data.frame(NULL)
    N = data.frame(NULL)
    for(b in 1:50){
      corrected =  beta + (beta.temp[b,] - beta.mean)
      corrected[1:2] = beta[1:2]
      corrected[corrected<0] = 0.001
      beta.corrected = rbind(beta.corrected,corrected)
      vseidr = sim.vseidr.det(T.max - 1,y = y,A = A,alpha = alpha.temp,beta = as.vector(t(corrected)),gamma.d = gamma.d,gamma.r = gamma.r,pt = pt,M = population.temp)
      N = rbind(N,vseidr$N)
    }
    beta.all = c(beta.all,discrete.integral(beta))
    beta.corrected.summary = data.frame(date = migration.temp$date,beta = beta,L = apply(X = beta.corrected,MARGIN = 2,FUN = quantile,0.025,na.rm = T),U = apply(X = beta.corrected,MARGIN = 2,FUN = quantile,0.975,na.rm = T),province = province.temp)
    beta.summary = rbind(beta.summary,beta.corrected.summary)
    summary.temp = data.frame(date = migration.temp$date,N = colMeans(N,na.rm = T),U = apply(X = N,MARGIN = 2,FUN = quantile,0.975),L = apply(X = N,MARGIN = 2,FUN = quantile,0.025),type = "nolockdown")
    summary.real <- results.temp %>% filter(date <= end.date) %>%dplyr::select(date,N)%>% mutate(date = date,N = N,U = N,L = N,type = "lockdown")
    summary = rbind(summary.temp,summary.real)
    summary$province = province.temp
    data.summary = rbind(data.summary,summary)
  }
  summarise.table = data.final <- data.summary %>% group_by(province,type) %>% dplyr::summarise(N = max(N)) %>%
    group_by(province) %>% dplyr::summarise(Lockdown = min(N), Nolockdown = max(N)) %>% mutate(Increase = Nolockdown - Lockdown,percentage = (Nolockdown - Lockdown)/Lockdown * 100 )
  summarise.table[,2:4] = floor( summarise.table[,2:4])
  p = ggplot(data.summary,aes(date,N,color = type))+facet_wrap(~province,scales = "free")+geom_line()+scale_color_npg()+
    geom_ribbon(aes(ymin = L,ymax = U),alpha = 0.2)+theme_minimal()+theme(legend.position ="none",axis.title = element_text(size = 8), 
                                                                          strip.text = element_text(size = 10,face = 'bold'),
                                                                          axis.title.x = element_blank(),
                                                                          axis.title.y = element_blank(),
                                                                          axis.text.x = element_text(size=8, angle=45,hjust = 1, face = 'bold'),
                                                                          axis.text.y = element_text(size=10, face = 'bold'),
                                                                          plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                          legend.title = element_blank(),
                                                                          legend.text = element_text(size=5, face = 'bold'),
                                                                          legend.key.width  = unit(.3,"inches"),
                                                                          legend.key.height = unit(.3,"inches"))
  p.beta = ggplot(beta.summary,aes(date,beta))+facet_wrap(~province,scales = "free")+geom_line()+scale_color_lancet()+
    geom_ribbon(aes(ymin = L,ymax = U),alpha = 0.2)+theme_minimal()+theme(legend.position ="none",axis.title = element_text(size = 8), 
                                                                          strip.text = element_text(size = 10,face = 'bold'),
                                                                          axis.title.x = element_blank(),
                                                                          axis.title.y = element_blank(),
                                                                          axis.text.x = element_text(size=8, angle=45,hjust = 1, face = 'bold'),
                                                                          axis.text.y = element_text(size=10, face = 'bold'),
                                                                          plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                          legend.title = element_blank(),
                                                                          legend.text = element_text(size=5, face = 'bold'),
                                                                          legend.key.width  = unit(.3,"inches"),
                                                                          legend.key.height = unit(.3,"inches"))
  beta.table = data.frame(province = provinces$provinces,beta.integral = beta.all)
  return(list(p = p, import = data.frame(province = provinces$provinces,import = import),table = summarise.table,p.beta = p.beta, beta.table = beta.table))
}
results.total = scenario.domestic.true(provinces,start,opendir,migration,end.date = ymd("2020-02-21"),w = 7)
ggsave("/Users/gujia/Documents/National Bureau of Statistics/Output/Nt2019.png")
ggsave("/Users/gujia/Documents/National Bureau of Statistics/Output/beta.png")










