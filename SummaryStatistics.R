rm(list = ls())
setwd("/Users/gujia/Documents/统计局/Data/")
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
library(corrplot)
#------summary statistics and plots---------
year2020 <- read.csv("year2020_Wuhan.csv")
year2020 <- year2020 %>% mutate(province = as.character(province)) %>% filter(province != "Hubei")
#mean size of migration population between Jan 10th and Jan 23rd, which is two weeks 
year2020.normal <- year2020 %>% mutate(date = ymd(date)) %>% filter(date <= ymd("20200123"),date >= ymd("20200110")) %>%
  group_by(province) %>%  dplyr::summarise(migration.normal = mean(daily_sum))
#mean size of migration after lockdown
year2020.lockdown <- year2020 %>% mutate(date = ymd(date)) %>% filter(date >= ymd("20200124"),date <= ymd("20200206")) %>%
  group_by(province) %>% dplyr::summarise(migration.lockdown = mean(daily_sum))
#change in migration (percentage)
year2020.combine <- left_join(year2020.normal,year2020.lockdown,by = "province")
year2020.summary <- year2020.combine %>% mutate(percentage = (migration.normal - migration.lockdown) / migration.normal,reduction = migration.normal - migration.lockdown,type = "Migration")
setwd("/Users/gujia/Documents/统计局/Output/")
ggplot(year2020.summary,aes(x = type,y = percentage))+geom_boxplot()+theme_bw()+
  geom_text_repel(aes(x = type,label = province))+theme(legend.position = "none",axis.title = element_text(size = 13), 
                                                        strip.text = element_text(size = 20,face = 'bold'),
                                                        axis.title.x = element_blank(),
                                                        axis.title.y = element_text(size=20),
                                                        axis.text.x = element_text(size=15,hjust = 0.5, face = 'bold'),
                                                        axis.text.y = element_text(size=20, face = 'bold'),
                                                        plot.title = element_text(size=20, face = 'bold', hjust = 0.5),
                                                        legend.text = element_text(size=20, face = 'bold'),
                                                        legend.background = element_blank(),
                                                        legend.title = element_text(size=20, face = 'bold'),
                                                        legend.key.width  = unit(.3,"inches"),
                                                        legend.key.height = unit(.3,"inches"))
ggsave("percentage_change.png")

ggplot(year2020.summary,aes(x = type,y = log(reduction)))+geom_boxplot()+theme_bw()+
  geom_text_repel(aes(x = type,label = province))+theme(legend.position = "none",axis.title = element_text(size = 13), 
                                                        strip.text = element_text(size = 20,face = 'bold'),
                                                        axis.title.x = element_blank(),
                                                        axis.title.y = element_text(size=20),
                                                        axis.text.x = element_text(size=15,hjust = 0.5, face = 'bold'),
                                                        axis.text.y = element_text(size=20, face = 'bold'),
                                                        plot.title = element_text(size=20, face = 'bold', hjust = 0.5),
                                                        legend.text = element_text(size=20, face = 'bold'),
                                                        legend.background = element_blank(),
                                                        legend.title = element_text(size=20, face = 'bold'),
                                                        legend.key.width  = unit(.3,"inches"),
                                                        legend.key.height = unit(.3,"inches"))



ggsave("reduction.png")
#-------correlation between the migration and Confirmed Cases --------
setwd("/Users/gujia/Documents/统计局/Data/nCoV_Pinyin0418/")
province = year2020.summary$province
province.remove = province 
cor.matrix <- matrix(ncol = 14, nrow = 28)
for(i in 1:length(province.remove)){
  print(province.remove[i])
  infected = read.csv(paste0(province.remove[i],".csv")) %>% mutate(date = mdy(as.character(date))) %>%
    filter(date >= mdy("01242020"),date <= mdy("02202020")) %>% dplyr::select(infected)
  flow <- (year2020 %>% mutate(date = ymd(date)) %>% filter(province == province.remove[i],date >= ymd("20200110"),date <= ymd("20200123")) %>%
    dplyr::select(daily_sum))[,1] %>% rev %>% cumsum
  if(province.remove[i] == "Xizang"){
    infected = rep(1,28)
  }
  Infect.data <- data.frame(infected = infected)
  if(i == 1){
    Infect.province <- Infect.data
    Flow.province <- data.frame(flow = flow)
  }
  else{
    Infect.province <- cbind(Infect.province, Infect.data)
    Flow.province <- cbind(Flow.province,flow)
  }
}
Infect.province <- as.matrix(Infect.province)
Flow.province <- as.matrix(Flow.province)
###Calculate the correlation
for(i in 1:28){
  for(j in 1:14){
    cor.matrix[i,j] = cor(Infect.province[i,],Flow.province[j,])
  }
}
corrplot(t(cor.matrix),method = c("number"))
#diagnostic
temp <- data.frame(province = province.remove,flow = Flow.province[14,],Infect.province[28,])

cor.dat <- data.frame(cor.matrix)
colnames(cor.dat) <- 0:13
cor.dat$Time <- ymd("20200124") + 0:27
cor.dat.part <- cor.dat[,c(1:7,15)]
cor.dat.melt <- melt(cor.dat.part,id.vars = "Time",variable.name = "Days",value.name = "Correlation")
Sys.setlocale(category = "LC_ALL",locale = "en_US.UTF-8")
ggplot(cor.dat.melt,aes(Time,Correlation,color = Days))+scale_color_npg()+geom_line()+scale_x_date(date_labels = "%b-%d")+theme_minimal()+theme(axis.title = element_text(size = 8),
                                                                                                            strip.text = element_text(size = 15,face = 'bold'),
                                                                                                            axis.title.x = element_text(size = 15),
                                                                                                            axis.title.y = element_text(size = 15),
                                                                                                            axis.text.x = element_text(size=15,hjust = 1, face = 'bold'),
                                                                                                            axis.text.y = element_text(size=15, face = 'bold'),
                                                                                                            plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                                                            legend.title = element_blank(),
                                                                                                            legend.text = element_text(size=15, face = 'bold'),
                                                                                                            legend.key.width  = unit(.3,"inches"),
                                                                                                            legend.key.height = unit(.3,"inches"),
                                                                                                            legend.position = c(0.75,0.3))



ggsave("/Users/gujia/Documents/统计局/Output/Cor_Migration_Cases_line.png")
Infected.firstweek = c()
Infected.total = c()
setwd("/Users/gujia/Documents/统计局/Data/nCoV_Pinyin0418/")
for(i in 1:length(province)){
  print(province[i])
  temp = read.csv(paste0(province[i],".csv")) %>% mutate(date = mdy(as.character(date))) %>%
    filter(date >= mdy("01312020"),date <= mdy("03152020")) %>% dplyr::select(infected)
  print(min(temp))
  print(max(temp))
  Infected.firstweek = c(Infected.firstweek,min(temp))
  Infected.total = c(Infected.total,max(temp))
}
Infected = data.frame(province = province,firstweek = Infected.firstweek,total = Infected.total) %>% mutate(province = as.character(province))
year2020.summary <- left_join(year2020.summary, Infected ,by = "province") %>% arrange(desc(migration.normal))

#-----first week-----
cor.test(log(year2020.summary$migration.normal),year2020.summary$firstweek,alternative = "greater") 
year2020.summary.part = year2020.summary[!(year2020.summary$province %in% c("Zhejiang","Guangdong")),]
cor(log(year2020.summary.part$migration.normal),year2020.summary.part$firstweek) #0.89
ggplot(year2020.summary,aes(log(migration.normal),firstweek))+geom_point()+
  geom_text_repel(aes(x = log(migration.normal),y = firstweek,label = province))+theme_bw()+labs(x = "log(population outflow)")+theme(axis.title = element_text(size = 8),
                                                                                                   strip.text = element_text(size = 15,face = 'bold'),
                                                                                                   axis.title.x = element_text(size = 15),
                                                                                                   axis.title.y = element_text(size = 15),
                                                                                                   axis.text.x = element_text(size=15,hjust = 1, face = 'bold'),
                                                                                                   axis.text.y = element_text(size=15, face = 'bold'),
                                                                                                   plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                                                   legend.title = element_blank(),
                                                                                                   legend.text = element_text(size=15, face = 'bold'),
                                                                                                   legend.key.width  = unit(.3,"inches"),
                                                                                                   legend.key.height = unit(.3,"inches"),
                                                                                                   legend.position = c(0.75,0.3))



path0 = "/Users/gujia/Documents/统计局/Output/"
ggsave(paste0(path0,"firstweek_Wuhan.png"))

#-----total------
cor.test(log(year2020.summary$migration.normal),year2020.summary$total,alternative = "greater") #0.74
ggplot(year2020.summary,aes(log(migration.normal),total))+geom_point()+
  geom_text_repel(aes(x = log(migration.normal),y = total,label = province))+theme_bw()+labs(x = "log(population outflow)")+theme(axis.title = element_text(size = 8),
                                                                                              strip.text = element_text(size = 15,face = 'bold'),
                                                                                              axis.title.x = element_text(size = 15),
                                                                                              axis.title.y = element_text(size = 15),
                                                                                              axis.text.x = element_text(size=15,hjust = 1, face = 'bold'),
                                                                                              axis.text.y = element_text(size=15, face = 'bold'),
                                                                                              plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                                              legend.title = element_blank(),
                                                                                              legend.text = element_text(size=15, face = 'bold'),
                                                                                              legend.key.width  = unit(.3,"inches"),
                                                                                              legend.key.height = unit(.3,"inches"),
                                                                                              legend.position = c(0.75,0.3))



ggsave(paste0("/Users/gujia/Documents/统计局/Output/","total_Wuhan.png"))
year2020.summary$migration.lockdown[27] = 1
cor(log(year2020.summary$migration.lockdown),year2020.summary$total) #0.77

 ggplot(year2020.summary,aes(log(migration.lockdown),total))+geom_point()+
  geom_text_repel(aes(x = log(migration.lockdown),y = total,label = province))+theme_bw()+labs(x = "log(population outflow)")+theme(axis.title = element_text(size = 8),
                                                                                                                                    strip.text = element_text(size = 15,face = 'bold'),
                                                                                                                                    axis.title.x = element_text(size = 15),
                                                                                                                                    axis.title.y = element_text(size = 15),
                                                                                                                                    axis.text.x = element_text(size=15,hjust = 1, face = 'bold'),
                                                                                                                                    axis.text.y = element_text(size=15, face = 'bold'),
                                                                                                                                    plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                                                                                    legend.title = element_blank(),
                                                                                                                                    legend.text = element_text(size=15, face = 'bold'),
                                                                                                                                    legend.key.width  = unit(.3,"inches"),
                                                                                                                                    legend.key.height = unit(.3,"inches"),
                                                                                                                                    legend.position = c(0.75,0.3))
 
 
 
ggsave(paste0(path0,"total_lockdown.png"))
year2020.summary$percentage[27] = 1
cor(year2020.summary$percentage,year2020.summary$total) #-0.62
ggplot(year2020.summary,aes(percentage,total))+geom_point()+
  geom_text_repel(aes(x = percentage,y = total,label = province))+theme_bw()
ggsave(paste0("/Users/gujia/Documents/统计局/Output/","total_percentage.png"))

#------butterfly------
library(Cairo)
#利用cairo_pdf输出中文字体
library(mapdata)
library(maptools)
library(ggplot2)
library(plyr)
library(scales)
library(colorspace)
library("ggplot2")
library("dplyr")
library("grid")
library("showtext")
library("Cairo")
year2020.butterfly <- year2020.summary %>% mutate(id = 1:30)
normal.thousand = c()
lockdown.thousand = c()
for(i in 1:dim(year2020.butterfly)[1]){
  temp.normal = year2020.butterfly[i,]$migration.normal 
  temp.lockdown = year2020.butterfly[i,]$migration.lockdown
  if(temp.normal > 100){
    normal.thousand = c(normal.thousand,paste0(as.character(round(temp.normal / 1000,1)),"k"))
  }
  else{
    normal.thousand = c(normal.thousand,"< 0.1k")
  }
  if(temp.lockdown > 100){
    lockdown.thousand = c(lockdown.thousand,paste0(as.character(round(temp.lockdown / 1000,1)),"k"))
  }
  else{
    lockdown.thousand = c(lockdown.thousand,"< 0.1k")
  }
  
}
year2020.butterfly$normal.thousand = normal.thousand
year2020.butterfly$lockdown.thousand = lockdown.thousand
p1<-ggplot(year2020.butterfly)+ geom_hline(yintercept=-mean(year2020.butterfly$migration.normal),linetype=2,size=.25,colour="grey")+
  geom_bar(aes(x=id,y=-migration.normal),stat="identity",fill="#C44E4C",colour=NA)+
  geom_text(aes(x=id,y=12000,label=province),vjust=.5,size = 4.5)+
  scale_x_reverse()+
  geom_text(aes(x=id,y=-(migration.normal + 10000),label=normal.thousand),size=4.5)+
  coord_flip()+
  theme_void()+labs(title = "(a) before lockdown")+theme(plot.title = element_text(hjust = 0.5,size = 25))
p2<-ggplot(year2020.butterfly)+
  geom_hline(yintercept=mean(year2020.butterfly$migration.lockdown),linetype=2,size=.25,colour="grey")+
  geom_bar(aes(x=id,y=migration.lockdown),stat="identity",fill="#E2BB1E",colour=NA)+
  scale_x_reverse()+
  geom_text(aes(x=id,y=migration.lockdown+3000,label=lockdown.thousand),size=4.5)+
  coord_flip()+
  theme_void()+labs(title = "(b) after lockdown")+theme(plot.title = element_text(hjust = 0.5,size = 25))
setwd("/Users/gujia/Documents/统计局/Output/")
p3 = ggarrange(p1,p2,ncol = 2,widths = c(2,1))
ggsave(filename = "butterfly.png",plot = p3,width = 15,height = 10,dpi = 300)
# CairoJPEG(file="butterfly.jpeg",width=1200,height=696,dpi = 300)
# showtext.begin()
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(7,11)))
# vplayout<-function(x,y){viewport(layout.pos.row =x,layout.pos.col=y)}
# print(p2,vp=vplayout(2:7,9:11))
# print(p1,vp=vplayout(2:7,1:8))
# grid.text(label="(b) After lockdown",x=.80,y=.88,gp=gpar(col="black",fontsize=25,fontfamily="sans",draw=TRUE),just="left")
# grid.text(label="(a) Before lockdown",x=.20,y=.88,gp=gpar(col="black",fontsize=25,fontfamily="sans",draw=TRUE),just="left")
# grid.text(label="Population flow from Hubei province before and after Jan 23rd",x=.50,y=.95,gp=gpar(col="black",fontsize=30,fontfamily="sans",draw=TRUE,just="centre"))
# showtext.end()
# dev.off()
#--------maps----------
library("maptools")
china_map = readShapePoly("/Users/gujia/Documents/统计局/Data/china-province-border-data/bou2_4p.shp") 
ggplot(china_map,aes(x=long,y=lat,group=group)) +
  geom_polygon(fill="white",colour="grey") +
  coord_map("polyconic")
x <- china_map@data
xs <- data.frame(x,id=seq(0:924)-1) 
china_map1 <- fortify(china_map) 
library(plyr)
china_map_data <- join(china_map1, xs, type = "full")
setwd("/Users/gujia/Documents/统计局/Output/")
#write.csv(china_map_data,"China.csv")
china_map_data = read.csv("China.csv",encoding = "UTF-8") %>% mutate(province = as.character(NAME))
china_data <- join(china_map_data, year2020.summary, type="full")

midpos <- function(x) mean(range(x,na.rm=TRUE))#取形状内的平均坐标
centres <- ddply(china_data,.(province),colwise(midpos,.(long,lat)))
p.map <- ggplot(china_data,aes(x = long,y = lat))+
  geom_polygon(aes(group = group,fill = percentage),colour="grey40") +
  scale_fill_gradient(low="red",high="white") +  
  coord_map("polyconic")+
  geom_text(aes(label=province),data=centres,size = 6)+
  theme(               #清除不需要的元素
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size=20, face = 'bold'),
    legend.background = element_blank(),
    legend.title = element_text(size=20, face = 'bold'),
    legend.key.width  = unit(.3,"inches"),
    legend.key.height = unit(.3,"inches")
  )
ggsave("percentage_map.png",plot = p.map,width = 20,height = 15,dpi = 300)
# p.normal = ggplot(china_data,aes(x = long,y = lat))+
#   geom_polygon(aes(group = group,fill = migration.normal),colour="grey40") +
#   scale_fill_gradient(low="white",high="steelblue",trans = 'log10',breaks = c(1e2,1e3,1e4,1e5)) +  
#   coord_map("polyconic")+
#   geom_text(aes(label=province),data=centres,size = 6)+
#   guides(fill=guide_legend(title="Before lockdown"))+
#   theme(               #清除不需要的元素
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     legend.text = element_text(size=20, face = 'bold'),
#     legend.background = element_blank(),
#     legend.title = element_text(size=20, face = 'bold'),
#     legend.key.width  = unit(.3,"inches"),
#     legend.key.height = unit(.3,"inches")
#   )
# 
# p.lockdown = ggplot(china_data,aes(x = long,y = lat))+
#   geom_polygon(aes(group = group,fill = migration.lockdown),colour="grey40") +
#   scale_fill_gradient(low="white",high="steelblue",trans = 'log10',breaks = c(1e2,1e3,1e4,1e5)) +  
#   coord_map("polyconic")+
#   geom_text(aes(label=province),data=centres,size = 6)+
#   guides(fill=guide_legend(title="After lockdown"))+
#   theme(               #清除不需要的元素
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     legend.text = element_text(size=20, face = 'bold'),
#     legend.background = element_blank(),
#     legend.title = element_text(size=20, face = 'bold'),
#     legend.key.width  = unit(.3,"inches"),
#     legend.key.height = unit(.3,"inches")
#   )

china_data.melt = melt(china_data,id.vars = c("group","long","lat","province"),measure.vars = c("migration.normal","migration.lockdown"))
p.melt = ggplot(china_data.melt,aes(x = long,y = lat))+facet_grid(~variable)+
  geom_polygon(aes(group = group,fill = value),colour="grey40") +
  scale_fill_gradient(low="white",high="steelblue",trans = 'log10',breaks = c(1e2,1e3,1e4,1e5)) +  
  coord_map("polyconic")+
  #geom_text(aes(label=province),data=centres,size = 6)+
  guides(fill=guide_legend(title="population size"))+
  theme(               #清除不需要的元素
    strip.text = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size=20, face = 'bold'),
    legend.background = element_blank(),
    legend.title = element_text(size=20, face = 'bold'),
    legend.key.width  = unit(.3,"inches"),
    legend.key.height = unit(.3,"inches")
  )
ggsave("change_map_Wuhan.png",plot = p.melt,width = 20,height = 15,dpi = 300)
