rm(list = ls())
source("/Users/gujia/Documents/统计局/Data/Data_prepare.R")
# w=11; r=5; D=14; alpha=rep(0.125,nrow(A4)); M=47100000
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
sim.vseidr.det = function(t_max, pt,y,A, alpha, beta, gamma.d,gamma.r,M,r.beta = 5) {
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

#--------Estimation of pE----------
#pE stay as a constant from Jan 10th to Jan 23th , then gradually drops to zero (till Jan 30th)
# alpha:diagnosis rate, set to be zero before Jan 20th, since there're almost no confirmed 
#A = read.csv("/Users/gujia/Documents/统计局/Data/year2020.csv")%>%mutate(province = as.character(province)) %>% filter(province == "Guangdong") %>% mutate(date = ymd(date))
path0 = "/Users/gujia/Documents/统计局/Data/nCoV_Pinyin0418/"
p.t = c(mdy("01102020"),mdy("01232020"),mdy("01162020"))
pE <- function(p.E,path0,province,p.t,A,M){
  opendir = paste0(path0,province,".csv")
  A4 = vSEIDR_Data(opendir = opendir, h = 5.5, st = mdy("01202020"),end = mdy("03152020"))
  lambda = 1/14 * log(200)
  p.seq = p.E / 100 * exp(lambda * 0:13)
  A.seq = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% dplyr::select(daily_sum)
  E0 = sum(A.seq * p.seq)
  # n1 = A %>% filter(date >= p.t[1],date <= p.t[3]) %>% dplyr::select(daily_sum) %>% sum()
  # n2 = A %>% filter(date > p.t[3],date < p.t[2]) %>% dplyr::select(daily_sum) %>% sum()
  # E0 = p.E/2*n1 + p.E*n2
  # n = A %>% filter(date >= p.t[1],date < p.t[2]) %>% dplyr::select(daily_sum) %>% sum()
  # E0 = n * p.E
  A4 = left_join(A4,A,by = c("date","province"))
  t = p.t[2] + (1:21)
  lambda2 = 1/7 * log(200)
  pt = data.frame(date = t,pt = p.E / 2 - p.E / 400 * exp(lambda2 * (1:7)))
  alpha = c(rep(1/14,mdy("01232020") - A4$date[1]+1),1/14 + cumsum(rep(3/14/14,14)))
  alpha = c(alpha,rep(1/3.5,dim(A4)[1] - length(alpha)))
  Coef = vSEIDR_Coef_new(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,pt = pt,r = 5)
  Coef.sim <- Coef %>% filter(date<=ymd("20200220"),date >= ymd("20200123"))
  sim = sim.vseidr.det(t_max = dim(Coef.sim)[1],pt = Coef.sim$pt,y = c(E0,Coef.sim$I[1],Coef.sim$R.d[1],Coef.sim$R.r[1]),
                 A = c(0,Coef.sim$A[-1]),alpha = Coef.sim$alpha,beta = Coef.sim$beta.hat, gamma.d  = Coef.sim$gamma.d,gamma.r = Coef.sim$gamma.r,M = M)
  N.part = (Coef %>% filter(date<=ymd("20200221"),date >= ymd("20200124")))$N
  sqrt(mean(abs(sim$N[-1] - N.part)^2))
}
fitting <- function(p.E,path0,province,p.t,A,M){
  opendir = paste0(path0,province,".csv")
  A4 = vSEIDR_Data(opendir = opendir, h = 5.5, st = mdy("01202020"),end = mdy("03152020"))
  lambda = 1/14 * log(200)
  p.seq = p.E / 100 * exp(lambda * 0:13)
  A.seq = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% dplyr::select(daily_sum)
  E0 = sum(A.seq * p.seq)
  # n1 = A %>% filter(date >= p.t[1],date <= p.t[3]) %>% dplyr::select(daily_sum) %>% sum()
  # n2 = A %>% filter(date > p.t[3],date < p.t[2]) %>% dplyr::select(daily_sum) %>% sum()
  # E0 = p.E/2*n1 + p.E*n2
  # n = A %>% filter(date >= p.t[1],date < p.t[2]) %>% dplyr::select(daily_sum) %>% sum()
  # E0 = n * p.E
  A4 = left_join(A4,A,by = c("date","province"))
  t = p.t[2] + (1:21)
  lambda2 = 1/7 * log(200)
  pt = data.frame(date = t,pt = p.E / 2 - p.E / 400 * exp(lambda2 * (1:7)))
  alpha = c(rep(1/14,mdy("01232020") - A4$date[1]+1),1/14 + cumsum(rep(3/14/14,14)))
  alpha = c(alpha,rep(1/3.5,dim(A4)[1] - length(alpha)))
  Coef = vSEIDR_Coef_new(A4 = A4,w = 7,D = 14,alpha = alpha,M = M,pt = pt,r = 5)
  Coef.sim <- Coef %>% filter(date<=ymd("20200220"),date >= ymd("20200123"))
  sim = sim.vseidr.det(t_max = dim(Coef.sim)[1],pt = Coef.sim$pt,y = c(E0,Coef.sim$I[1],Coef.sim$R.d[1],Coef.sim$R.r[1]),
                       A = Coef.sim$A,alpha = Coef.sim$alpha,beta = Coef.sim$beta.hat, gamma.d  = Coef.sim$gamma.d,gamma.r = Coef.sim$gamma.r,M = M)
  N.part = (Coef %>% filter(date<=ymd("20200221"),date >= ymd("20200124")))$N
  data1 = data.frame(date = ymd("20200124")+(0:28),N = sim$N[-1],province = province,type = "fitting")
  data2 = data.frame(date = ymd("20200124")+(0:28),N = N.part,province = province,type = "observed")
  list(fitting = rbind(data1,data2),Coef = Coef)
}
provinces = c("Anhui","Beijing","Chongqin","Fujian","Guangdong","Guangxi","Hebei","Heilongjiang","Henan","Hunan",
              "Jiangsu","Jiangxi","Shaanxi","Shandong","Shanghai","Sichuan","Zhejiang")
population.all = read.csv("/Users/gujia/Documents/统计局/Data/country_population.csv")
pt.all = c()
E = c()
N = c()
Gof = c()
Import.afterlockdown = c()
fitting.performance = data.frame(NULL)
Coef = data.frame(NULL)
for(i in 1:length(provinces)){
  print(i)
  province.temp = provinces[i]
  opendir = paste0(path0,province.temp,".csv")
  data = read.csv(opendir) %>% mutate(date = mdy(as.character(date)))
  N = c(N,(data$infected)[data$date == ymd("20200315")])
  population = population.all$population[population.all$city == province.temp] * 10000
  A = read.csv("/Users/gujia/Documents/统计局/Data/year2020.csv")%>%mutate(province = as.character(province)) %>% filter(province == province.temp) %>% mutate(date = ymd(date))
  if(province.temp %in% c("Jiangxi")){
    percentage = optimize(pE,interval = c(1/7500,1/1000),path0 = path0,province = province.temp,p.t = p.t, A = A,M = population)$minimum
  }
  else{
  if(province.temp %in% c("Hunan","Anhui")){
    percentage = optimize(pE,interval = c(1/25000,1/1000),path0 = path0,province = province.temp,p.t = p.t, A = A,M = population)$minimum
  }
  else{
    percentage = optimize(pE,interval = c(1/100000,0.1),path0 = path0,province = province.temp,p.t = p.t, A = A,M = population)$minimum
  }
  }
  Gof = c(Gof,pE(percentage,path0 = path0,province = province.temp,p.t = p.t, A = A, M = population))
  fitting.results = fitting(percentage,path0 = path0,province = province.temp,p.t = p.t, A = A, M = population)
  fitting.performance = rbind(fitting.performance,fitting.results$fitting)
  Coef = rbind(Coef,fitting.results$Coef)
  # n = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% select(daily_sum) %>% sum()
  # n = n * percentage
  # n1 = A %>% filter(date >= p.t[1],date <= p.t[3]) %>% select(daily_sum) %>% sum()
  # n2 = A %>% filter(date > p.t[3],date < p.t[2]) %>% select(daily_sum) %>% sum()
  # n = percentage/2*n1 + percentage*n2
  lambda = 1/14 * log(200)
  p.seq = percentage / 100 * exp(lambda * 0:13)
  A.seq = A %>% filter(date >= p.t[1],date <= p.t[2]) %>% dplyr::select(daily_sum)
  n = sum(A.seq * p.seq)
  E = c(E,n)
  A.seq.1 = A %>% filter(date > p.t[2],date <= p.t[2] + 7) %>% dplyr::select(daily_sum)
  lambda2 = 1/7 * log(200)
  p.seq.1 = percentage / 2 - percentage / 400 * exp(lambda2 * (1:7))
  Import.afterlockdown = c(Import.afterlockdown,sum(A.seq.1 * p.seq.1))
  pt.all = c(pt.all,percentage)
}
summary.pt = data.frame(province =provinces,E = E,N = N,pt = pt.all)
summary.pt <- summary.pt %>% mutate(E = as.integer(E),proportion = E/N,fitting.performance = Gof,Import.cont = as.integer(Import.afterlockdown))
summary.pt
xtable(summary.pt[,c(1,6,2,7,3)])

fitting.plot <-ggplot(fitting.performance,aes(date,N,color = type))+facet_wrap(~province,scales = "free")+geom_line()+scale_color_npg()+
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
ggsave("/Users/gujia/Documents/统计局/Output/fitting.png")
write.csv(summary.pt,"/Users/gujia/Documents/统计局/Output/pt.csv")

save.dir = "/Users/gujia/Documents/统计局/Output/"
write.csv(Coef,paste0(save.dir,"Results0602.csv"))
Coef<-Coef %>% mutate(date = ymd(date)) %>% filter(date >= ymd("20200124"),date <= ymd("20200221"))
ggplot(Coef,aes(date,Rt))+facet_wrap(~province,scales = "free")+scale_color_npg()+geom_line()+theme_minimal()+
  theme(legend.position ="none",axis.title = element_text(size = 8),
        strip.text = element_text(size = 10,face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, angle=45,hjust = 1, face = 'bold'),
        axis.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=5, face = 'bold'),
        legend.key.width  = unit(.3,"inches"),
        legend.key.height = unit(.3,"inches"))
ggsave(paste0(save.dir,"Rt.png"))
