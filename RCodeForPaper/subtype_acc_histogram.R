subtypes = c('H1N1','H1N2','H2N1','H2N2','H3N2','H3N8','H5N1','H7N2','H7N9','H9N2','H10N3','H10N8')
subtype_acc_segs_plot = subtype_acc_segs[subtype_acc_segs$subtype %in% subtypes,]
subtype_acc_segs_plot$Segment = factor(subtype_acc_segs_plot$Segment,levels = c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS','Mix'))
subtype_acc_segs_plot$subtype = factor(subtype_acc_segs_plot$subtype,levels = c('H1N1','H1N2','H2N1','H2N2','H3N2','H3N8','H5N1','H7N2','H7N9','H9N2','H10N3','H10N8'))
# windowsFonts(Times_New_Roman=windowsFont("Times New Roman"))
ggplot(subtype_acc_segs_plot,aes(Segment,Accuracy,fill=Method))+
  geom_bar(stat="identity",position="dodge",alpha=0.8) +  theme_light() + 
  scale_fill_manual(values=color_set)+theme(plot.title = element_text(face='bold.italic',size = 20),
                                            axis.title.x = element_text(face = 'bold',size =30),
                                            axis.title.y = element_text(face = 'bold',size=30),
                                            legend.title = element_text(face = 'bold',size=20),
                                            legend.text = element_text(face = 'bold',size=20),
                                            strip.text = element_text(face = 'bold.italic',size=20),
                                            axis.line = element_line(color = "black",
                                                                     size = 1, linetype = "solid"),
                                            axis.text.x = element_text(size = 13,
                                                                       face = "bold"),
                                            axis.text.y = element_text(size = 13, 
                                                                       face = "bold"),
                                            legend.position = "right",
                                            legend.background = element_rect(colour = 'grey90',linewidth = 2, size=.5),
                                            axis.ticks = element_line(size = 1),
                                            axis.ticks.length = unit(5, "pt"),
  ) + facet_wrap(~subtype)
ggsave(filename = 'results/subtype_seg_bar.pdf',width = 16,height = 8.5,limitsize = FALSE)
