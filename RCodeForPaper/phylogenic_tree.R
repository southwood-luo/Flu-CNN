library(ggtreeExtra)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(treeio)
library(tidytree)
library(dplyr)
library(ggstar)
library(TDbook)
library(tidyr)
library(stringr)

rm(list = ls())
seg = 'PA'
# read tree and metadata----
tr = read.tree(paste0('data/treeflie/H5N1-H7N9-H9N2_',seg,'_ALL_mafft.fasta.treefile'))
labels <- tr$tip.label
label_parts = as.vector(label_parts)
metada = read.table('data/treeflie/meta.dat',sep = '\t',quote = '')
colnames(metada) = c('genomeID','strain_name','country','location','subtype','year','host','host_group')
metada$id = labels

# subtype
table(metada$subtype)
table(metada$host_group)
colors1 <- c('#ff7f0e', '#2ca02c', '#9467bd')
metada <- metada %>%
  mutate(subtype__colour = ifelse(subtype == "H5N1", colors1[1], 
                                  ifelse(subtype == "H7N9", colors1[2],colors1[3])))


#country
length(table(metada$country))
n_colors <- 100 
colors_vec <- sample(c("#000000", colors()), n_colors, replace = TRUE)
metada$country__colour <- colors_vec[(match(metada$country, sort(unique(metada$country)))) %% n_colors + 1]

metada$year <- ifelse(is.na(metada$year), 2020,metada$year)
metada$year_group <- ifelse(metada$year %in% c(1979:2000), 1,
                            ifelse(metada$year %in% c(2001:2005), 2,
                                   ifelse(metada$year %in% c(2006:2010), 3,
                                          ifelse(metada$year %in% c(2011:2015), 4,
                                                 ifelse(metada$year %in% c(2016:2020), 5,
                                                        6)))))
metadata_df = metada
year_groups <- unique(metadata_df$year_group)
colors_year <- c("#9467bd", "#FFA500", "#8c564b", "#2ca02c", "#4169E1", "#FF1493")
for (i in seq_along(year_groups)) {
  y <- year_groups[i]
  year_temp = metadata_df[metadata_df$year_group == y,]
  year_group_temp <- as.numeric(names(table(year_temp$year)))
  year_group_size <- length(year_group_temp)
  color_temp = colors_year[i]
  rgb_color <- col2rgb(color_temp, alpha = TRUE)
  alpha_per <- floor(255/year_group_size)
  for (j in seq_along(year_group_temp)) {
    hex_color <- rgb(rgb_color[1], rgb_color[2], rgb_color[3],
                     alpha=alpha_per*j, maxColorValue=255)
    year = year_group_temp[j]
    metadata_df$year__colour = ifelse(metadata_df$year == year, hex_color,metadata_df$year__colour)
  }
}
metada$year__colour = metadata_df$year__colour

#host
table(metada$host_group)
colors3 <- c("#A50F15", "#08519C")
metada <- metada %>%
  mutate(host_group__colour = ifelse(host_group == "Human", colors3[1], colors3[2]))

metadata <- metada %>%
  select(c("id", "subtype", "subtype__colour", 
           "country", "country__colour",
           "year", "year__colour", 
           "host", "host__colour",
           "host_group", "host_group__colour"))

countrycolors <- metada %>%
  select(c("country", "country__colour")) %>%
  distinct()

subtypecolors <- metada %>%
  select(c("subtype", "subtype__colour")) %>%
  distinct()

yearcolors <- metada %>%
  select(c("year", "year__colour")) %>%
  distinct()
yearcolors <- yearcolors[order(yearcolors$year, decreasing=TRUE),]

hostcolors <- metada %>%
  select(c("host", "host__colour")) %>%
  distinct()

hostgroupcolors <- metada %>%
  select(c("host_group", "host_group__colour")) %>%
  distinct()

metadata$country <- factor(metadata$country, levels=countrycolors$country)
metadata$year <- factor(metadata$year, levels=yearcolors$year)
metadata$host_group <- factor(metadata$host_group, levels=hostgroupcolors$host_group)
metadata$subtype <- factor(metadata$subtype, levels=subtypecolors$subtype)

pp <- ggtree(tr, layout="fan", open.angle=15, size=0.01)

pp <- pp %<+% metadata

pp1 <-pp +
  geom_tippoint(
    mapping=aes(colour=subtype),
    size=0.3,
    stroke=0,
    alpha=1
  ) +
  scale_colour_manual(
    name="Subtype",
    values=subtypecolors$subtype__colour,
    guide=guide_legend(keywidth=0.3,
                       keyheight=0.3,
                       ncol=2,
                       override.aes=list(size=2,alpha=1),
                       order=1)
  ) +
  theme(
    legend.title=element_text(size=5),
    legend.text=element_text(size=4),
    legend.spacing.y = unit(0.02, "cm")
  )

pp2 <-pp1 +
  geom_fruit(
    geom=geom_star,
    mapping=aes(fill=host_group),
    starshape=13,
    color=NA,
    size=2.5,
    starstroke=0,
    offset=0.1,
  ) +
  scale_fill_manual(
    name="Host",
    values=hostgroupcolors$host_group__colour,
    guide=guide_legend(
      keywidth=0.3,
      keyheight=0.3,
      order=3
    ),
    na.translate=FALSE
  )

pp3 <-pp2 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=year),
    width=0.015,
    offset=0.1
  ) +
  scale_fill_manual(
    name="Year",
    values=yearcolors$year__colour,
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
  ) +
  theme(
    legend.title=element_text(size=6), 
    legend.text=element_text(size=4.5),
    legend.spacing.y = unit(0.02, "cm")
  )

ggsave(filename = paste0('results/',seg,'_branch.pdf'), plot = pp3,width = 7,height = 7)
