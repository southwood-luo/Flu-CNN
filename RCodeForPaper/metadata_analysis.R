library(ggplot2)
library(ggpubr)

## read metadata----
getwd()
df = read.csv('data/complete_acc/complete_meta.csv',row.names = 1)
df = read.csv('data/complete_acc/metadata.dat',sep = '\t',row.names = 1)
df[1,]
colnames(df) = c('')
table(df$Host)

# line chart----
df1 = df[df$year %in% c(1933:2022),]
events_by_year <- table(df1$year) 
data = data.frame(Year = names(events_by_year), Quantity = as.numeric(events_by_year))
str(data)
ggplot(data, aes(x = as.numeric(Year), y = Quantity)) + 
  geom_line(color='#0072B2', size=1.25) +
  scale_x_continuous(breaks=seq(1900,2025,5)) +
  scale_y_continuous(breaks=seq(0, 11000, 1000)) +
  xlab("Year") + ylab("Quantity") +
  theme(panel.background = element_rect(fill='white', color='grey50')) +
  theme(panel.grid.major = element_line(color='grey80', linetype='dashed')) +
  theme(axis.line = element_line(color='black')) +
  theme(axis.text = element_text(size = 13)) +
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold",color='#666666'))
ggsave(filename = 'results/complete_meta/line_chart.png',width = 12,height = 5)
dev.off()

## subtype ----
df = read.csv('data/complete_test_acc/complete_meta.csv',row.names = 1)
df$host_group = str_to_title(df$host_group)
subtype_table = data.frame(table(df$subtype))
colnames(subtype_table) = c('Subtype','Quantity')
subtype_table = subtype_table[order(subtype_table$Quantity),]
subtype_big = as.character(subtype_table$Subtype[subtype_table$Quantity > 500])
df$subtype2 = ifelse(df$subtype %in% subtype_big,df$subtype,'Others')
df$subtype2 = factor(df$subtype2,levels=c('Others',subtype_big))
str(df)
colnames(df)[2] = 'Host'
windowsFonts(Times_New_Roman=windowsFont("Times New Roman"))
pdf("results/complete_meta/subtype_meta1.pdf",width = 6,height = 10)
ggplot(df,aes(x=subtype2))+
  geom_bar(aes(fill=Host),position = "stack",width = .9,color = "black")+
  scale_x_discrete("",position = "top")+ 
  scale_fill_manual(values=c("#6495ED",'#B22222'))+
  scale_y_reverse("Total strain number",expand = c(0,0),position = "right",limits=c(35000,0))+
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = size_use),
    # legend.text = element_text(size = size_use),legend.title = element_text(size = size_use),
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+#theme_bw()+
  coord_flip()+theme(plot.margin=unit(rep(1,4),'lines'))
dev.off()
ggsave("results/complete_meta/subtype_meta1.png",width = 8,height = 12,units = "cm")

pdf("results/complete_test/complete/complete_meta/subtype_meta2.pdf",width = 6,height = 10)
size_use = 22
ggplot(df,aes(x=subtype2))+
  geom_bar(aes(fill=Host),position = "fill",width = .9,color = "black")+
  scale_x_discrete("")+
  scale_y_continuous("Total strain proportion",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual(values=c("#6495ED",'#B22222'))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = size_use),
    legend.text = element_text(size = size_use),legend.title = element_text(size = size_use),
    axis.title.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  coord_flip()
dev.off()
ggsave("results/complete_test/complete/complete_meta/subtype_meta2.png",width = 8,height = 12,units = "cm")


## location_pie_chart----
library(ggplot2);library(scatterpie)
worldMap6 <- fortify(map_data("world"))
table(worldMap6$region)
worldMap6$fill<-ifelse(worldMap6$region=="Cambodia"| worldMap6$region=="Mexico"| worldMap6$region=="Egypt"| 
                         worldMap6$region=="Slovakia"| worldMap6$region=="Thailand"| worldMap6$region=="Vietnam"| 
                         worldMap6$region=="Somalia"| worldMap6$region=="Korea"| worldMap6$region=="Greece"| 
                         worldMap6$region=="Singapore"| worldMap6$region=="Ghana"| worldMap6$region=="Saudi Arabia"| 
                         worldMap6$region=="Switzerland"| worldMap6$region=="Colombia"| worldMap6$region=="Brazil"| 
                         worldMap6$region=="Turkey"| worldMap6$region=="UK"| worldMap6$region=="Pakistan"| 
                         worldMap6$region=="Israel"| worldMap6$region=="Iraq"| worldMap6$region=="Bangladesh"| 
                         worldMap6$region=="Senegal"| worldMap6$region=="Georgia"| worldMap6$region=="Finland"| 
                         worldMap6$region=="Kyrgyzstan"| worldMap6$region=="Belarus"| worldMap6$region=="Burundi"| 
                         worldMap6$region=="Nigeria"| worldMap6$region=="Canada"| worldMap6$region=="Guinea-Bissau"| 
                         worldMap6$region=="Bolivia"| worldMap6$region=="Democratic Republic of the Congo"| 
                         worldMap6$region=="Oman"| worldMap6$region=="Afghanistan"| worldMap6$region=="Qatar"| 
                         worldMap6$region=="Croatia"| worldMap6$region=="Portugal"| worldMap6$region=="Romania"| 
                         worldMap6$region=="Poland"| worldMap6$region=="Czech Republic"| worldMap6$region=="Germany"| 
                         worldMap6$region=="Chile"| worldMap6$region=="Netherlands"| worldMap6$region=="Kenya"| 
                         worldMap6$region=="Ireland"| worldMap6$region=="Laos"| worldMap6$region=="Kuwait"| 
                         worldMap6$region=="Australia"| worldMap6$region=="Sri Lanka"| worldMap6$region=="New Zealand"| 
                         worldMap6$region=="Serbia"| worldMap6$region=="Mali"| worldMap6$region=="China"| 
                         worldMap6$region=="Malaysia"| worldMap6$region=="Republic of the Congo"| worldMap6$region=="Slovenia"| 
                         worldMap6$region=="Guatemala"| worldMap6$region=="Madagascar"| worldMap6$region=="Argentina"| 
                         worldMap6$region=="Rwanda"| worldMap6$region=="Tanzania"| worldMap6$region=="Belgium"| 
                         worldMap6$region=="Haiti"| worldMap6$region=="USA"| worldMap6$region=="Hungary"| 
                         worldMap6$region=="India"| worldMap6$region=="Lebanon"| worldMap6$region=="Gambia"| 
                         worldMap6$region=="Russia"| worldMap6$region=="Uruguay"| worldMap6$region=="South Korea"| 
                         worldMap6$region=="South Africa"| worldMap6$region=="Armenia"| worldMap6$region=="Jordan"| 
                         worldMap6$region=="Sierra Leone"| worldMap6$region=="Cuba"| worldMap6$region=="Nepal"| 
                         worldMap6$region=="Malawi"| worldMap6$region=="Iran"| worldMap6$region=="Norway"| 
                         worldMap6$region=="Myanmar"| worldMap6$region=="Peru"| worldMap6$region=="Kazakhstan"| 
                         worldMap6$region=="United Arab Emirates"| worldMap6$region=="Spain"| worldMap6$region=="Indonesia"| 
                         worldMap6$region=="France"| worldMap6$region=="Denmark"| worldMap6$region=="Austria"| 
                         worldMap6$region=="Japan"| worldMap6$region=="Algeria"| worldMap6$region=="Italy"| 
                         worldMap6$region=="Mozambique"| worldMap6$region=="Azerbaijan"| worldMap6$region=="Sweden"| 
                         worldMap6$region=="Ecuador","A","B")
mp0<-ggplot(data=worldMap6,aes(x=long,y=lat,group=group,fill=fill))+geom_polygon(colour='grey55')
mp3<-mp0+ylim(-60,90)
mp3
data33 <- read.csv('data/data_result.csv')
data2<- data.frame(data33)
# data3<-t(data2[5])
data3<-t(data2[7])
data4<-data3[1:96]
data2$radius<-data4
mp3+geom_scatterpie(aes(x=long, y=lat, group=region, r=radius/200),geom_line(lwd = 0.000000001), 
                    data=data2,cols=names(data2[4:6]), alpha=.8)+
  scale_fill_manual(values=c('grey70','white',"#00BA38","#619CFF","#F8766D"),
                    labels=c("Source country","Not found","Single replicon", 
                             "Double replicons","More than two replicons"),name="") + 
  coord_equal()+geom_scatterpie_legend(data2$radius/200, x=210, y=-30,n=5,labeller = function(x) {200*x})+
  theme(legend.text=element_text(size=16),legend.position = c(0.915,0.6))+
  theme(axis.title = element_blank(),axis.text.x =element_blank(),axis.ticks.x = element_blank(),
        axis.text.y =element_blank(),axis.ticks.y = element_blank())+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="white"))
