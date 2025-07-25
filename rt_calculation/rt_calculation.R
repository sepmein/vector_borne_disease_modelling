library(EpiEstim)
library(tidyverse)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)

#设置路径
path <- "F:\\AI学习"

d1<-read.csv("F:\\AI学习\\daily_cases_foshan.csv",header=TRUE)

d1$date <- as.Date(d1$date)

names(d1) <- c("dates", "I")

# Rt estimation -----------------------------------------------------------
#set time window
t_start <- seq(4, length(d1$I) - 2)
t_end <- t_start + 2

## setting config
conf <- make_config(
  list(
    mean_si = 3.5,
    std_si = 2.6, 
    t_start = t_start,
    t_end = t_end
  )
)

## setting method
meth <- "parametric_si"

## estimating Rt
DataRt <- estimate_R(
  d1,
  method = meth,
  config = conf
)

#result
Rt_result <- data.frame(date=d1$dates[DataRt$R$t_end],
                        meanR=c(round(DataRt$R$`Mean(R)`,4)),
                        lbd=c(DataRt$R$`Quantile.0.025(R)`),
                        ubd=c(DataRt$R$`Quantile.0.975(R)`)) %>% 
  mutate(CI = paste0(round(meanR,4)," (",round(lbd,4),", ",round(ubd,4),")"))

#plot
ggplot(data = Rt_result,aes(x=date,y=meanR))+
  geom_ribbon(aes(ymin=lbd,ymax=ubd),fill= "#AD002AFF",alpha=0.2)+
  geom_line(size=1,colour= "#AD002AFF")+
  geom_hline(yintercept = 1,size=1,lty=2)+
  scale_x_date(date_breaks = "1 days",date_labels = "%m/%d",expand = c(0,0), limits=as.Date(c("2025-07-13","2025-07-22")))+
  labs(x="Date",y="Estimated Rt")+
  theme(legend.title = element_blank(),legend.position = c(0.85,0.95),
        legend.background = element_blank(),legend.key.size = unit(15,"pt"),
        legend.key = element_blank(),legend.text=element_text(size=15,hjust = 0),
        axis.text = element_text(size = 20,colour = "black"),
        axis.title = element_text(size=20,color = "black"),
        axis.line = element_line(colour = "black",size = 1),
        axis.ticks = element_line(color = "black",size = 1),
        panel.background = element_blank(),panel.grid = element_line(colour = "grey"),
        panel.border = element_rect(fill = NA,size = 1,colour = "black"),
        plot.margin=unit(rep(2,4),'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size=20,colour = "black")
  )

ggsave(filename = paste0("F:\\AI学习\\rt_报告数.png"),width = 16,height = 9,dpi = 300)
write.xlsx(Rt_result,paste0(path,"rt_result.xlsx"))
