# PlotNeutral function ----------------------------------------------------------------

#Example: PlotNeutral(core_otus)

# This function will plot the result of a neutral model fit, and will color
# the core otus that are driven by the host (bove the model) or dispersal limited
# (below the model). The neutral OTUs are plotted as grey circles.

PlotNeutral <- function(obs){
  obs1 <- 
    as.data.frame(obs[2])
  obs1 <- obs1[!is.na(obs1$p), ]
  obs2 <- 
    as.data.frame(obs[1])
  plotn <- ggplot(obs1, aes(x=log10(otu_rel), y=otu_occ)) +
    # geom_point(data=subset(obs1, fill%in%"no" | fit_class%in%"As predicted"),
    #            aes(x=log10(otu_rel), y=otu_occ), pch=21, fill="white", alpha=0.25, size=1.5) +
    geom_point(data=subset(obs1, fill%in%"no" | fit_class%in%"As predicted"),
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill="white", alpha=0.25, size=3) +
    geom_point(data=subset(obs1, fill%in%"core" & fit_class%in%"As predicted"),
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill="gray60", color="gray60", size=3) +
    geom_point(data=subset(obs1, fill%in%"core" & fit_class%in%"Above prediction"), 
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='#CB54D6', color= '#CB54D6',size=3) +
    geom_point(data=subset(obs1, fill%in%"core" & fit_class%in%"Below prediction"), 
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='#09104d', color= '#09104d', size=3) +
    geom_line(color='red', data=obs1, size=0.7, aes(y=obs1$freq.pred, x=log10(obs1$p)), alpha=0.25) +
    geom_line(color='black', lty='twodash', size=0.7, data=obs1, aes(y=obs1$pred.upr, x=log10(obs1$p)), alpha=0.25)+
    geom_line(color='black', lty='twodash', size=0.7, data=obs1, aes(y=obs1$pred.lwr, x=log10(obs1$p)), alpha=0.25)+
    labs(x="Log10(mean abundance)", y="Occupancy") + #mean relative abundance
    theme(axis.title = element_text(angle = 0, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position = "right") +
    annotate("text", -Inf, Inf, label = paste("italic(R)^2 ==", round(obs2$Rsqr, 3)),
             parse = TRUE, size = 8, hjust = -0.2, vjust = 1.2) +
    annotate("text", -Inf, Inf, label = paste("italic(m) ==", round(obs2$m, 3)),
             parse = TRUE, size = 8, hjust = -0.2, vjust = 3.2)+theme_cowplot(font_size = 26)
  return(plotn)
}
