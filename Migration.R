library(tidyverse)
library(lubridate)
library(pinyin)
library(smooth)
library(Mcomp)
library(ggplot2)
library(openxlsx)
opendata = "/Users/gujia/Documents/统计局/Data/nCoV_Pinyin0418/"
import_data = "export_to_province_new.xlsx"
data = read.xlsx(paste0(opendata,import_data),sheet = 3)
data <- data %>% mutate(date = ymd(date))

##setting 1 shutdown
provinces = data.frame(provinces = c("Zhejiang","Guangdong","Jiangsu","Anhui",
                                     "Jiangxi","Beijing","Heilongjiang","Shanghai","Henan","Shaanxi","Chongqin","Hunan"))
data1 = data %>% mutate(province = as.character(province))
data1 = data1[,-2]
data.shut = data.frame(NULL)
for(i in 1:dim(provinces)[1]){
  data.temp = data1 %>% filter(province == provinces$provinces[i],date <= ymd("2020-02-29"))
  #last = data.temp[dim(data.temp)[1],]
  #time.seq = last$date + 1:14
  #migration.seq = last$daily_sum - cumsum(rep(last$daily_sum/14,14))
  #add = data.frame(province = provinces$provinces[i],date = time.seq,daily_sum = migration.seq)
  #data.temp = rbind(data.temp,add)
  data.shut = rbind(data.shut,data.temp)
}
data.shut$type = "shutdown"




data1 = data %>% mutate(province = as.character(province))
data1 = data1[,-2]
data.noshut = data.frame(NULL)
for(i in 1:dim(provinces)[1]){
  data.temp = data1 %>% filter(province == provinces$provinces[i])
  data.temp = data.temp %>%filter(date <= ymd("2020-01-23"))
  start.time = data.temp %>%filter(date %in% (ymd("2020-01-18")+ 1:5))
  ave = mean(start.time$daily_sum)
  #constant+gradual decline
  final = data1 %>% filter(province == provinces$provinces[i],date == ymd("2020-01-23")+28)
  dif = ave - final$daily_sum
  migration.seq = c(rep(ave,14),ave - cumsum(rep(dif/14,14)))
  time.seq = ymd("2020-01-23")+1:28
  add = data.frame(province = provinces$provinces[i],date = time.seq,daily_sum = migration.seq)
  data.temp = rbind(data.temp,add)
  data.noshut = rbind(data.noshut,data.temp)
}
data.noshut$type = "noshutdown"

data.combine = rbind(data.shut,data.noshut)
ggplot(data.combine,aes(date,daily_sum,color = type))+facet_wrap(~province,scales = "free")+theme_bw()+geom_line()+theme(legend.position ="none",axis.title = element_text(size = 8), 
                                 strip.text = element_text(size = 10,face = 'bold'),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_text(size=10, angle=45,hjust = 1, face = 'bold'),
                                 axis.text.y = element_text(size=10, face = 'bold'),
                                 plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                 legend.title = element_blank(),
                                 legend.text = element_text(size=5, face = 'bold'),
                                 legend.key.width  = unit(.3,"inches"),
                                 legend.key.height = unit(.3,"inches"))+
  geom_vline(xintercept = ymd("2020-01-23"),linetype = 2,color = "orange")+scale_x_date(breaks = ymd("2020-01-23") + (-3:5)*7,date_labels=format("%m/%d"))


Output = "/Users/gujia/Documents/统计局/Output/"
ggsave(paste0(Output,"Scenario_Migration.png"))
write.csv(data.combine,paste0(Output,"Migration.csv"))


data = read.csv("/Users/gujia/Documents/统计局/Data/nCoV_Pinyin0418/Migration_full.csv")
data = data[data$type != "noshutdown",]
data$type = as.character(data$type)
data$type[data$type == "shutdown"] = "lockdown" 
data <- data %>% mutate(date = ymd(date))
ggplot(data,aes(date,daily_sum,color = type))+facet_wrap(~province,scales = "free")+theme_minimal()+scale_color_npg()+geom_line()+theme(legend.position ="none",axis.title = element_text(size = 8), 
                                                                                                                         strip.text = element_text(size = 10,face = 'bold'),
                                                                                                                         axis.title.x = element_blank(),
                                                                                                                         axis.title.y = element_blank(),
                                                                                                                         axis.text.x = element_text(size=10, angle=45,hjust = 1, face = 'bold'),
                                                                                                                         axis.text.y = element_text(size=10, face = 'bold'),
                                                                                                                         plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                                                                         legend.title = element_blank(),
                                                                                                                         legend.text = element_text(size=5, face = 'bold'),
                                                                                                                         legend.key.width  = unit(.3,"inches"),
                                                                                                                         legend.key.height = unit(.3,"inches"))+
  geom_vline(xintercept = ymd("2020-01-23"),linetype = 2,color = "orange")+scale_x_date(breaks = ymd("2020-01-23") + (-3:5)*7,date_labels=format("%m/%d"))
ggsave("/Users/gujia/Documents/统计局/Output/2019.png")

year2019<-read.csv("/Users/gujia/Documents/统计局/Data/year2019_Wuhan.csv")
year2020<-read.csv("/Users/gujia/Documents/统计局/Data/year2020_Wuhan.csv")

year2020 <- year2020 %>% mutate(province = as.character(province),date = ymd(date))
year2020$type = "2020"

year2019 <- year2019 %>% mutate(province = as.character(province),date = ymd(date))
year2019$type = "2019"

two.year <- rbind(year2019,year2020)

two.year <- two.year %>% filter(province %in% provinces, date >= ymd("20200120"),date <= ymd("20200221"))

ggplot(two.year,aes(date,daily_sum,color = type))+facet_wrap(~province,scales = "free")+theme_minimal()+scale_color_npg()+geom_line()+theme(legend.position ="none",axis.title = element_text(size = 8), 
                                                                                                                                        strip.text = element_text(size = 10,face = 'bold'),
                                                                                                                                        axis.title.x = element_blank(),
                                                                                                                                        axis.title.y = element_blank(),
                                                                                                                                        axis.text.x = element_text(size=10, angle=45,hjust = 1, face = 'bold'),
                                                                                                                                        axis.text.y = element_text(size=10, face = 'bold'),
                                                                                                                                        plot.title = element_text(size=10, face = 'bold', hjust = 0.5),
                                                                                                                                        legend.title = element_blank(),
                                                                                                                                        legend.text = element_text(size=5, face = 'bold'),
                                                                                                                                        legend.key.width  = unit(.3,"inches"),
                                                                                                                                        legend.key.height = unit(.3,"inches"))+
  geom_vline(xintercept = ymd("2020-01-23"),linetype = 2,color = "orange")+scale_x_date(breaks = ymd("2020-01-23") + (-3:5)*7,date_labels=format("%m/%d"))
ggsave("/Users/gujia/Documents/统计局/Output/2019_Wuhan.png")
