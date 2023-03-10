#Modified from https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R
# ExtractCore graphing functions ------------------------------------------------------------------------------------
# Modofied form Stopnisek and Shade 2019 - Current Opinion in Microbiology
# https://www.sciencedirect.com/science/article/pii/S1369527419300426?via%3Dihub

#Graphs the increase in the Bray-Curtis with the exclusion of OTU

#
#BC_ranked = ExtractCoreFlex generate data format 
#max = Max number of OTUs graphed
#per_cut = Percentage Bray-Curtis threshold used (whole number)
#x_breaks = Breaks for the x-axis
#ann_offset = x offset for the annotation on the graph

PlotBCincreaseFlex <- function(BC_ranked, max,per_cut,x_breaks,ann_offset){
  BC_ranked <- as.data.frame(BC_ranked[[2]])
  dim(BC_ranked) %>% print()
  BC_ranked$rank <- factor(BC_ranked$rank, levels=BC_ranked$rank)
  #BC_ranked <- BC_ranked[complete.cases(BC_ranked), ]
  plot <- ggplot(BC_ranked[1:max,], 
                 aes(x=rank[1:max], y=proportionBC)) +
    geom_point(pch=21, fill='white', alpha=0.25, size=1.5) +
    theme(strip.background = element_blank(),axis.text.x = element_text(size=6, angle=90)) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 16, face = "bold")) +
    
    theme(axis.text.y = element_text(angle = 0, size = 14, hjust = 0.5, vjust = 0.5)) +
    scale_x_discrete(limits=BC_ranked$rank[1:max], # making the x axis more readable
                     breaks=BC_ranked$rank[1:max][seq(1,length(BC_ranked$rank[1:max]),by=x_breaks)]) +
    geom_vline(xintercept=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.0+(per_cut*0.01))])), 
               lty=4, col="red", cex=0.5) +
    labs(x="Ranked OTUs",y="Bray-Curtis similarity") +
    annotate(geom="text", 
             x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.0+(per_cut*0.01))]))+ann_offset, 
             y=0.2, 
             label=paste("Last ",per_cut,"% increase\n(",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=
                                                                                           1+(per_cut*0.01))]))," OTUs)", 
                         sep=""), color="red", size=8, )+theme_cowplot(font_size = 26)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.15)) 
  return(plot)
}

#Graphing the richness of OTU table at each Bray-Curtis level

PlotBCThreshold_Rich<- function(BC_ranked,max_thresh,fin_thresh){
  BC_ranked_1 <- as.data.frame(BC_ranked[[2]])
  BC_ranked_1$rank <- factor(BC_ranked_1$rank, levels=BC_ranked_1$rank)
  thresh_rich<-data.frame(x_thresh=c(2:max_thresh),
                          y_rich=c(rep(NA)))
  
  for (i in 2:max_thresh) {
    thresh_rich$y_rich[i-1]<-last(as.numeric(BC_ranked_1$rank[(BC_ranked_1$IncreaseBC>=1.0+(i*0.01))]))/
      nrow(BC_ranked[[4]])
  }
  
  plot_rich=ggplot(thresh_rich,aes(x=x_thresh,y=y_rich))+geom_point(size=3)+geom_line()+
    geom_vline(xintercept=fin_thresh,color="red")+
    geom_text(x=fin_thresh+4,y=max(thresh_rich$y_rich)-quantile(thresh_rich$y_rich, probs =0.25),
              label=paste("Last ",fin_thresh,"% increase\n(",
                          last(as.numeric(BC_ranked_1$rank[(BC_ranked_1$IncreaseBC>=1+
                                                              (fin_thresh*0.01))]))," OTUs)",sep = ""),
              color="red",size=8)+
    ylab("Proportion of richness")+xlab("BC Threshold")+theme_cowplot(font_size = 25)
  return(plot_rich)
}

#Graphing the number of reads in table at each Bray-Curtis level

PlotBCThreshold_Abun<- function(BC_ranked,max_thresh,fin_thresh){
  BC_ranked_1 <- as.data.frame(BC_ranked[[2]])
  BC_ranked_1$rank <- factor(BC_ranked_1$rank, levels=BC_ranked_1$rank)
  BC_ranked_2 <- as.data.frame(BC_ranked[[4]])
  BC_ranked_3 <- as.data.frame(BC_ranked[[3]])
  thresh_abun<-data.frame(x_thresh=c(2:max_thresh),
                          y_abun=c(rep(NA)))
  
  for (i in 2:max_thresh) {
    thresh_abun$y_abun[i-1]<-sum(BC_ranked_2[BC_ranked_2$otu%in%BC_ranked_3$
                                               otu[1:last(as.numeric(BC_ranked_1$rank[(BC_ranked_1$IncreaseBC>=1.0+
                                                                                         (i*0.01))]))],"otu_rel"])
  }
  
  plot_abun=ggplot(thresh_abun,aes(x=x_thresh,y=y_abun))+geom_point(size=3)+geom_line()+
    geom_vline(xintercept=fin_thresh,color="red")+
    ylab("Proportion of reads")+xlab("BC Threshold")+theme_cowplot(font_size = 25)
  return(plot_abun)
}
