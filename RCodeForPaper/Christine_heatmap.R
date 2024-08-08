zoonotic_meta = read.csv('data/zoonotic/zoonotic_strain2.csv')
zoonotic_meta1 = zoonotic_meta[zoonotic_meta$Group == 'Avian-isolated suspected zoonotic',]
zoonotic_meta1 = zoonotic_meta1[order(zoonotic_meta1$Subtype,zoonotic_meta1$Country),]
zoonotic_meta1$genomeID = as.character(zoonotic_meta1$genomeID)
zoonotic_meta2 = zoonotic_meta[zoonotic_meta$Group == 'Human-isolated confirmed zoonotic',]
zoonotic_meta2 = zoonotic_meta2[order(zoonotic_meta2$Subtype,zoonotic_meta2$Country),]
zoonotic_meta2$genomeID = as.character(zoonotic_meta2$genomeID)
str(zoonotic_meta1)

library(ggplot2)
library(ComplexHeatmap)
library(circlize)

datacore1 = read.csv('data/zoonotic/zoonotic_strain_3parts_smallsub2.csv',row.names = 1)
# sort(datacore1$info)
datacore1 = datacore1[order(datacore1$info),]
data1 = datacore1[datacore1$info == 'Avian-isolated suspected zoonotic',]
data2 = datacore1[datacore1$info == 'Human-isolated confirmed zoonotic',]
data = datacore1[,1:8]

data[is.na(data)]=0
my36colors <- c('#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#B22222',"#6495ED",
                         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                         '#968175')
mycolor = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#7f7f7f')
str(datacore1)
table(datacore1$info)
data[data < 0.5]=0
data[data >= 0.5]=1

#avian_zoonotic
data1[data1 < 0.5]=0;data1[data1 >= 0.5]=1
data1 = data1[zoonotic_meta1$genomeID,]
png('results/zoonotic_heatmap/smallsub_zoonotic_avian2.png',width = 80,height = 3000)
Heatmap(as.matrix(data1[,1:8]),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = T,show_column_names = F,col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()

#human_zoonotic
data2[data2 < 0.5]=0;data2[data2 >= 0.5]=1
data2 = data2[zoonotic_meta2$genomeID,]
png('results/zoonotic_heatmap/smallsub_zoonotic_human.png',width = 80,height = 1000)
Heatmap(as.matrix(data2[,1:8]),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()


png('results/zoonotic_heatmap/smallsub_zoonotic.png',width = 300,height = 2000)
Heatmap(as.matrix(data),cluster_columns = FALSE,name = "Probability",cluster_rows = F,
        show_row_names = F,left_annotation = rowAnnotation(df = data.frame(Type=datacore1$info),width=unit(15,"cm"),
                                                           col = list(Type = c("Human-isolated confirmed zoonotic" = "#D6E7A3", 
                                                                               "Avian-isolated suspected zoonotic" = "#E63863"))), 
        # top_annotation = HeatmapAnnotation(Family = anno_barplot(orthonum9)),
        col = colorRamp2(c(0,0.5,1),c("white","#57C3F3",'#E95C59')),show_heatmap_legend = T,
        na_col = "white")
dev.off()
ggsave(filename = 'results/zoonotic_heatmap/smallsub_zoonotic.png')


##typical
unzoonotic = read.csv('data/zoonotic/zoonotic_strain_3parts_smallsub2_not.csv',row.names = 1)
unzoonotic = na.omit(unzoonotic)
unzoonotic1 = unzoonotic[unzoonotic$info == 'Avian',]
unzoonotic2 = unzoonotic[unzoonotic$info == 'Human',]

unzoonotic_meta1 = zoonotic_meta[zoonotic_meta$Group == 'Avian',]
unzoonotic_meta1 = unzoonotic_meta1[order(unzoonotic_meta1$Subtype,unzoonotic_meta1$Country),]
unzoonotic_meta1$genomeID = as.character(unzoonotic_meta1$genomeID)
unzoonotic_meta2 = zoonotic_meta[zoonotic_meta$Group == 'Human',]
unzoonotic_meta2 = unzoonotic_meta2[order(unzoonotic_meta2$Subtype,unzoonotic_meta2$Country),]
unzoonotic_meta2$genomeID = as.character(unzoonotic_meta2$genomeID)
str(unzoonotic_meta1)

#avian_unzoonotic
unzoonotic1[unzoonotic1 < 0.5]=0;unzoonotic1[unzoonotic1 >= 0.5]=1
unzoonotic1 = unzoonotic1[unzoonotic_meta1$genomeID,]
rownames(unzoonotic_meta1) = unzoonotic_meta1$genomeID
unzoonotic_meta1 = unzoonotic_meta1[rownames(unzoonotic1),]
unzoonotic1 = na.omit(unzoonotic1)
png('results/zoonotic_heatmap/smallsub_unzoonotic_avian.png',width = 80,height = 1000)
Heatmap(as.matrix(unzoonotic1[,1:8]),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()

#human_unzoonotic
unzoonotic2[unzoonotic2 < 0.5]=0;unzoonotic2[unzoonotic2 >= 0.5]=1
unzoonotic2 = unzoonotic2[unzoonotic_meta2$genomeID,]
rownames(unzoonotic_meta2) = unzoonotic_meta2$genomeID
unzoonotic_meta2 = unzoonotic_meta2[rownames(unzoonotic2),]
unzoonotic2 = na.omit(unzoonotic2)
png('results/zoonotic_heatmap/smallsub_unzoonotic_human.png',width = 80,height = 1000)
Heatmap(as.matrix(unzoonotic2[,1:8]),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()




##我的数据，所有----
my_all = read.csv('data/zoonotic/all_acc.csv',row.names = 1)
meta = read.csv('data/zoonotic/genomeIDToinfo.csv',row.names = 1)


#禽
avian_genomeID = meta[meta$host_group == 'avian',]
avian_genomeID = avian_genomeID[order(avian_genomeID$subtype,avian_genomeID$country),]
avian_genomeID = na.omit(avian_genomeID)
my_avian = my_all[rownames(avian_genomeID),]
my_avian = na.omit(my_avian)
avian_genomeID = avian_genomeID[rownames(my_avian),]
my_avian_data = my_avian[,4:11]
my_avian_data[my_avian_data < 0.5]=0;my_avian_data[my_avian_data >= 0.5]=1

str(avian_genomeID)
table(avian_genomeID$year)
avian_genomeID$country2 = ifelse(avian_genomeID$country %in% c('China','USA'),avian_genomeID$country,'others')
avian_genomeID$subtype2 = ifelse(avian_genomeID$subtype %in% c('H7N9','H9N2','H5N1','H3N2','H1N1','H2N2'),avian_genomeID$subtype,'others')
avian_genomeID$year2 = ifelse(avian_genomeID$year %in% c(as.character(1902:1990)),'-1990',
                              ifelse(avian_genomeID$year %in% c(as.character(1991:2000)),'1991-2000',
                                     ifelse(avian_genomeID$year %in% c(as.character(2001:2010)),'2001-2010',
                                            ifelse(avian_genomeID$year %in% c(as.character(2011:2015)),'2011-2015',
                                                   ifelse(avian_genomeID$year %in% c(as.character(2016:2020)),'2016-2020',
                                                          ifelse(avian_genomeID$year %in% c(as.character(2021:2022)),'2021-2022','others'))))))
avian_genomeID$year2 = ifelse(avian_genomeID$year %in% c(as.character(2018:2020)),avian_genomeID$year,
                             'others')

png('results/zoonotic_heatmap/smallsub_all_avian4.png',width = 3000,height = 200)
Heatmap(as.matrix(t(my_avian_data[10000:23665,])),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,column_names_rot = 90,
        top_annotation = columnAnnotation(df = data.frame(subtype=avian_genomeID$subtype2[10000:23665],
                                                          country=avian_genomeID$country2[10000:23665],
                                                          year=avian_genomeID$year2[10000:23665]),height=unit(50,"cm"),
                                          col = list(subtype = c("H7N9" = "#6495ED","H5N1" = "#E59CC4", "H3N2" = "#968175",
                                                              "H1N1" = "#625D9E", "H9N2" = "#E63863","others" = '#CCC9E6'),
                                                     country = c("China" = "#F3B1A0","USA" = "#E39A35","others" = '#CCC9E6'))), 
        # top_annotation = HeatmapAnnotation(Family = anno_barplot(orthonum9)),
        na_col = "white",
        col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()


#人
human_genomeID = meta[meta$host_group == 'human',]
human_genomeID = human_genomeID[order(human_genomeID$subtype,human_genomeID$country),]
human_genomeID = na.omit(human_genomeID)
my_human = my_all[rownames(human_genomeID),]
my_human = na.omit(my_human)
human_genomeID = human_genomeID[rownames(my_human),]
my_human_data = my_human[,4:11]
my_human_data[my_human_data < 0.5]=0;my_human_data[my_human_data >= 0.5]=1

str(human_genomeID)
table(human_genomeID$year2)
human_genomeID$country2 = ifelse(human_genomeID$country %in% c('China','USA'),human_genomeID$country,'Others')
human_genomeID$country2 = ifelse(human_genomeID$country %in% c('China','USA','Zambia','HongKong','Thailand','Indonesia',
                                                               'Vietnam','Egypt','Cambodia'),human_genomeID$country,'Others')

human_genomeID$subtype2 = ifelse(human_genomeID$subtype %in% c('H7N9','H9N2','H5N1','H3N2','H1N1','H2N2','H5N6'),human_genomeID$subtype,'Others')
human_genomeID$year2 = ifelse(human_genomeID$year %in% c(as.character(1905:1990)),'-1990',
                              ifelse(human_genomeID$year %in% c(as.character(1991:2000)),'1991-2000',
                                     ifelse(human_genomeID$year %in% c(as.character(2001:2010)),'2001-2010',
                                            ifelse(human_genomeID$year %in% c(as.character(2011:2020)),'2011-2020',
                                                    ifelse(human_genomeID$year %in% c(as.character(2021:2022)),'2021-2022','Others')))))
human_genomeID$year2 = ifelse(human_genomeID$year %in% c(as.character(2018:2020)),human_genomeID$year,
                              'others')

png('results/zoonotic_heatmap/smallsub_all_human3.png',width = 3000,height = 250)
Heatmap(as.matrix(t(my_human_data[,])),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,column_names_rot = 90,
        top_annotation = columnAnnotation(df = data.frame(subtype=human_genomeID$subtype2[],
                                                          country=human_genomeID$country2[],
                                                          year=human_genomeID$year2[]),height=unit(50,"cm"),
                                          col = list(subtype = c("H7N9" = "#6495ED","H5N1" = "#E59CC4", "H3N2" = "#968175","H2N2" = "#D6E7A3",
                                                                 "H1N1" = "#625D9E", "H9N2" = "#E63863","others" = '#CCC9E6'),
                                                     country = c("China" = "#B22222","USA" = "#6495ED","others" = '#CCC9E6'),
                                                     year = c('-1990' = '#476D87','1991-2000' = '#AB3282','2001-2010' = '#8C549C',
                                                              '2011-2020' = '#E4C755','2021-2022' = '#AA9A59',"others" = '#CCC9E6'))), 
        # top_annotation = HeatmapAnnotation(Family = anno_barplot(orthonum9)),
        na_col = "white",
        col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)

dev.off()
png('results/zoonotic_heatmap/smallsub_all_human.png',width = 3000,height = 200)
Heatmap(as.matrix(t(my_human_data[,])),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,column_names_rot = 90,
        # top_annotation = columnAnnotation(df = data.frame(subtype=human_genomeID$subtype2[],
        #                                                   country=human_genomeID$country2[],
        #                                                   year=human_genomeID$year2[]),height=unit(50,"cm"),
        #                                   col = list(subtype = c("H7N9" = "#6495ED","H5N1" = "#E59CC4", "H3N2" = "#968175","H2N2" = "#D6E7A3",
        #                                                          "H1N1" = "#625D9E", "H9N2" = "#E63863","others" = '#CCC9E6'),
        #                                              country = c("China" = "#B22222","USA" = "#6495ED","others" = '#CCC9E6'),
        #                                              year = c('-1990' = '#476D87','1991-2000' = '#AB3282','2001-2010' = '#8C549C',
        #                                                       '2011-2020' = '#E4C755','2021-2022' = '#AA9A59',"others" = '#CCC9E6'))), 
        # top_annotation = HeatmapAnnotation(Family = anno_barplot(orthonum9)),
        na_col = "white",
        col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()
png('results/zoonotic_heatmap/smallsub_all_human5.png',width = 1000,height = 50)
Heatmap(as.matrix(t(my_human_data[55338:56373,])),cluster_columns = FALSE,cluster_rows = F,
        show_row_names = F,show_column_names = F,column_names_rot = 90,
        # heatmap_legend_param = list(legend_direction = "left", title_position = "topcenter"),
        # top_annotation = columnAnnotation(df = data.frame(Subtype=human_genomeID$subtype2[55320:56373],
        #                                                   Country=human_genomeID$country2[55320:56373],
        #                                                   Year=human_genomeID$year2[55320:56373]),height=unit(20,"cm"),
        #                                   # annotation_legend_side = "left",
        #                                   col = list(Subtype = c("H7N9" = "#6495ED","H5N1" = "#E59CC4", "H3N2" = "#968175","H2N2" = "#D6E7A3",
        #                                                          "H1N1" = "#625D9E", "H9N2" = "#E63863",'H5N6' = '#91D0BE',"Others" = '#CCC9E6'),
        #                                              Country = c("China" = "#B22222","USA" = "#6495ED",'Zambia' = '#476D87',
        #                                                          'HongKong' = '#F3B1A0','Thailand' = '#D6E7A3','Indonesia' = '#E59CC4',
        #                                                          'Vietnam' = '#5F3D69','Egypt' = '#58A4C3','Cambodia' = '#E4C755',
        #                                                          "Others" = '#CCC9E6'),
        #                                              Year = c('-1990' = '#476D87','1991-2000' = '#57C3F3','2001-2010' = '#8C549C',
        #                                                       '2011-2020' = '#F7F398','2021-2022' = '#CCC9E6',"Others" = '#CCC9E6'))
        #                                   ),
        na_col = "white",
        col = colorRamp2(c(0,1),c("#6495ED",'#B22222')),
        show_heatmap_legend = F)
dev.off()
table(human_genomeID$subtype[55320:56373])
table(human_genomeID$country[55320:56373])
my36colors <- c('#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#B22222',"#6495ED",
                         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                         '#968175')


lgd = Legend(at = c(1, 2), title = "Legend")
draw(lgd, x = 0.9, y = 0.5)



