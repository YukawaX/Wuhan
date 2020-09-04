migration <- read.csv("/Users/gujia/Documents/统计局/Data/year2020_Wuhan.csv")
provinces = c("Anhui","Beijing","Chongqin","Fujian","Guangdong","Guangxi","Hebei","Heilongjiang","Henan","Hunan",
              "Jiangsu","Jiangxi","Shaanxi","Shandong","Shanghai","Sichuan","Zhejiang")
dat.combined <- data.frame(NULL)
for(i in 1:length(provinces)){
  Confirmed.temp <- read.csv(paste0("/Users/gujia/Documents/统计局/Data/nCoV_Pinyin0418/",provinces[i],".csv")) %>% dplyr::select(province,date,infected) %>%
    mutate(date = mdy(date),infected = c(diff(infected),NA))
  Migration.temp <- migration %>% filter(province == provinces[i]) %>% dplyr::select(province,date,daily_sum) %>% mutate(date = ymd(date))
  dat.temp <- left_join(Confirmed.temp,Migration.temp)
  dat.combined <- rbind(dat.combined,dat.temp)
}
dat.combined$daily_sum[is.na(dat.combined$daily_sum)] = 0
library(reshape2)
dat.long <- melt(dat.combined,id.vars = c("province","date"),variable.name = "type",value.name = "Number")
dat.long.ind <- dat.long %>% filter(date >= ymd("20200125"),date <= ymd("20200221"))
ggplot(dat.long.ind,aes(date,Number,color = type))+facet_wrap(~province,scales = "free")+geom_line()+scale_color_npg()+
  theme_minimal()+theme(legend.position = "none",
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
ggsave("./统计局/Output/Migration_Confirmed.png")
