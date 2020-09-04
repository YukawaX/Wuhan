rm(list = ls())
library(tidyverse)
library(doMC)
registerDoMC(4)
library(lubridate)
library(pinyin)
library(smooth)
library(Mcomp)
library(ggplot2)
library(openxlsx)
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
sim.vseidr.det = function(t_max, pt,y,A, alpha, beta, gamma.d,gamma.r,M,r.beta = 5) {
  pt = pt$pt
  t = 0 : t_max
  S = rep(0, t_max+1)
  E = rep(0, t_max+1)
  I = rep(0, t_max+1)
  R.d = rep(0, t_max+1)
  R.r = rep(0, t_max+1)
  pE = pt
  M = M + cumsum(A)
  #initial value 不动 
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
    
    # dS = min(rpois(1, dS),S[i])
    # dE = min(rpois(1, dE),E[i] + dS)
    # dR.d = min(rpois(1, dR.d), I[i] + dE)
    # dR.r = min(rpois(1, dR.r), I[i] + dE-dR.d)
    
    S[i+1] = S[i] - dS + (1 - pE[i]) * A[i] 
    E[i+1] = E[i] + dS - dE + pE[i] * A[i]
    I[i+1] = I[i] + dE - dR.d- dR.r 
    R.d[i+1] = R.d[i] + dR.d
    R.r[i+1] = R.r[i] + dR.r
  }
  return(data.frame(t = t, S = S, E = E, I = I, R = R.d+R.r, N = I+R.d+R.r,R.d = R.d, R.r = R.r))
}
#check = sim.vseir.det(19,y = c(100,10,5,5),pt = rep(0.0001,19),A = rep(1000,19),alpha = rep(0.2,19),beta = rep(0.15,19),gamma.d = rep(1/14,19),gamma.r = rep(1/14,19),M = 50000000)
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
  #initial value 不动 
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
  A[1] = 0
  A4$alpha = alpha
  A4$E = c(diff(A4$N) / alpha[-length(alpha)], NA)
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
  coef.one = coef.one %>% dplyr::select('N','I','R','R.d','R.r','E','S',
                                        'alpha','beta','beta.hat',
                                        'gamma','gamma.d','gamma.r',
                                        'Rt','pt','A')
  coef.one$Rt[which(coef.one$Rt<0)]=0
  coef.one$beta.hat[which(coef.one$beta.hat<0)] = 0.001
  return(coef.one)
}

#Scenario analysis for lockdown
provinces = data.frame(provinces = c("Anhui","Beijing","Chongqin","Fujian","Guangdong","Guangxi","Hebei","Heilongjiang","Henan","Hunan",
                                     "Jiangsu","Jiangxi","Shaanxi","Shandong","Shanghai","Sichuan","Zhejiang")) %>% mutate(provinces = as.character(provinces))
population = read.csv("/Users/gujia/Documents/统计局/Data/country_population.csv") %>% mutate(city = as.character(city))
colnames(population) = c("provinces","population")
provinces = left_join(provinces,population,by = "provinces")
start = read.xlsx("/Users/gujia/Documents/统计局/Data/start.xlsx")
start <- start %>% mutate(st = ymd(st),end = ymd(end))  #tuning parameter alpha
migration = read.csv("/Users/gujia/Documents/统计局/Data/year2019_Wuhan.csv")
opendir = "/Users/gujia/Documents/统计局/Output/Results0730.csv"



#-------scenario true comparison with year19 data------
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
  pt.full <-read.csv("/Users/gujia/Documents/统计局/Output/pt.csv")
  import = c()
  beta.all = c()
  for(i in 1:dim(provinces)[1]){
    print(i)
    province.temp = provinces$provinces[i]
    population.temp = provinces$population[i]
    population.temp = provinces$population[i] * 10000
    results.temp <- results %>% mutate(date = ymd(date)) %>% filter(province ==province.temp,date >= ymd("2020-01-23"),date<=end.date) 
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
ggsave("/Users/gujia/Documents/统计局/Output/Nt2019.png")
ggsave("/Users/gujia/Documents/统计局/Output/beta.png")


import.full <- left_join(summary.pt,results.total$import)
import.full <- import.full %>% dplyr::select(province,Import.cont,import) %>% mutate(import = as.integer(import))
import.full <- import.full %>% mutate(ratio = import / Import.cont)

results.total$table <- results.total$table %>% mutate(Lockdown = as.integer(Lockdown),Nolockdown = as.integer(Nolockdown), Increase = as.integer(Increase))
a = left_join(summary.pt,results.total$import)
b = left_join(a,results.total$table)
bb <- b %>% mutate(Increase.import = as.integer(import - Import.cont),Increase.confirmed = Increase,ratio = Increase.confirmed / Increase.import) %>% dplyr::select(province,Increase.import,Increase.confirmed,ratio)
bbb <- left_join(bb,results.total$beta.table)
xtable(bbb)

table.3 <- results.total$table %>% dplyr::select(province,Lockdown,Nolockdown) %>% mutate(ratio = Nolockdown / Lockdown)
xtable(table.3)
