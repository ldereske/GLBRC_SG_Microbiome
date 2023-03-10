
## ---------------------------
##
## Script name: Community dispersion of the switchgrass microbiomes cross 
## the Lux Arbor Growing Season
##
## Purpose of script: Graphing and characterization of the dispersion in 
## bacterial and fungal communities 
## of roots and soils of switchgrass (Panicum virgatum)
##
## Author: Lukas Bell-Dereske
##
## Email: lukas.dereske@gmail.com
##
## ---------------------------


library(here)
#here::i_am("R_code/OTU_bacterial_fungal_community_analyses_202110011.R")
library(phyloseq)
library(plyr); library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(ggsignif)
library(vegan)
library(car)

Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "English")



#Load community data and make the phyloseq object

#Bacteria
GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_OTU_tab=
  read.table(here::here("Publish_data","OTU_GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.txt"),
             sep = "\t", header = T, row.names = 1)


GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_tax_tab=
  as.matrix(read.csv(here::here("Publish_data","Tax_table_GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.csv"),
                     header = T, row.names = 1))
head(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_tax_tab)

GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_metadata=
  read.csv(here::here("Publish_data","metadata_GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.csv"),
           header = T, row.names = 1)

GLBRC018_OTU_bact_MMPRNT_mock_bact_rar=
  phyloseq(otu_table(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_OTU_tab,taxa_are_rows = T),
           tax_table(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_tax_tab),
           sample_data(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_metadata))
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)
#1109
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)
#31839
rm(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_OTU_tab)
rm(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_tax_tab)
rm(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar_metadata)


GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_OTU_tab=
  read.table(here::here("Publish_data","OTU_GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.txt"),
             sep = "\t", header = T, row.names = 1)

GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_OTU_tab[1:10,1:10]


GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_tax_tab=
  as.matrix(read.csv(here::here("Publish_data","Tax_table_GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.csv"),
                     header = T, row.names = 1))
head(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_tax_tab)

GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_metadata=
  read.csv(here::here("Publish_data","metadata_GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.csv"),
           header = T, row.names = 1)

GLBRC018_OTU_fung_MMPRNT_mock_fung_rar=
  phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_OTU_tab,taxa_are_rows = T),
           tax_table(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_tax_tab),
           sample_data(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_metadata))

nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#1086
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#3090
rm(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_OTU_tab)
rm(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_tax_tab)
rm(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar_metadata)



#####Lux Arbor Community Composition #####

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,plotType=="G5"&siteID=="LUX")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
#496
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5) > 0, GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
ntaxa(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
#23991
rm(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,plotType=="G5"&siteID=="LUX")
nsamples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
#490
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5) > 0, GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
ntaxa(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
#2723
rm(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)

#Amp129 is an outlier so I am going to remove it 

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,SampleID!="Amp129")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
#495
rm(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)





#Lux Arbor Betadisp####


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,method = "bray")
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,method = "bray")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betamod=betadisper(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis, with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,interaction(Root_soil,collectionDate,FertStatus)))
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis=as.data.frame(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betamod$distances)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,by="row.names")
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis)[colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis)==
                                                        "GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betamod$distances"]="betadisp"


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betamod=betadisper(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis, with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map,interaction(Root_soil,collectionDate,FertStatus)))
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis=as.data.frame(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betamod$distances)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map,by="row.names")
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis)[colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis)==
                                                        "GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betamod$distances"]="betadisp"




#Stats

#Roots
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis_root=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Root")

root_betadisp_LUX_bact_mod=lmer((betadisp)^-1~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis_root)
plot(root_betadisp_LUX_bact_mod)
hist(resid(root_betadisp_LUX_bact_mod))
qqPlot(resid(root_betadisp_LUX_bact_mod))
shapiro.test(resid(root_betadisp_LUX_bact_mod))
#W = 0.99228, p-value = 0.6366
simulateResiduals(fittedModel = root_betadisp_LUX_bact_mod, plot = T)

anova(root_betadisp_LUX_bact_mod)
#collectionDate            1.78757 0.35751     5 119.12  2.7225 0.02292 *

#emmeans(root_betadisp_LUX_bact_mod,pairwise~FertStatus|collectionDate)

emmeans(root_betadisp_LUX_bact_mod,pairwise~collectionDate)

#Soils
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_soil=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Soil")

soil_betadisp_LUX_bact_mod=lmer((betadisp)^-1~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_soil)
plot(soil_betadisp_LUX_bact_mod)
hist(resid(soil_betadisp_LUX_bact_mod))
qqPlot(resid(soil_betadisp_LUX_bact_mod))
shapiro.test(resid(soil_betadisp_LUX_bact_mod))
#W = 0.99328, p-value = 0.1169
simulateResiduals(fittedModel = soil_betadisp_LUX_bact_mod, plot = T)


anova(soil_betadisp_LUX_bact_mod, type = 3)
#collectionDate            3.08433 0.220309    14 312.06  2.8424 0.0004906 ***
#FertStatus                0.22953 0.229533     1 312.04  2.9614 0.0862655 .  
#collectionDate:FertStatus 1.16370 0.083122    14 312.03  1.0724 0.3820224     

emmeans(soil_betadisp_LUX_bact_mod,pairwise~FertStatus|collectionDate)

emmeans(soil_betadisp_LUX_bact_mod,pairwise~collectionDate)

#Roots
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis_root=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Root")

root_betadisp_LUX_fung_mod=lmer(log(betadisp)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis_root)
plot(root_betadisp_LUX_fung_mod)
hist(resid(root_betadisp_LUX_fung_mod))
qqPlot(resid(root_betadisp_LUX_fung_mod))
shapiro.test(resid(root_betadisp_LUX_fung_mod))
#W = 0.99043, p-value = 0.4506
simulateResiduals(fittedModel = root_betadisp_LUX_fung_mod, plot = T)

anova(root_betadisp_LUX_fung_mod)
#collectionDate            0.45116 0.090231     5 118.26  2.3862 0.04217 *

emmeans(root_betadisp_LUX_fung_mod,pairwise~collectionDate)

#Soils
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis_soil=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Soil")

soil_betadisp_LUX_fung_mod=lmer(log(betadisp)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis_soil)
plot(soil_betadisp_LUX_fung_mod)
hist(resid(soil_betadisp_LUX_fung_mod))
qqPlot(resid(soil_betadisp_LUX_fung_mod))
shapiro.test(resid(soil_betadisp_LUX_fung_mod))
#W = 0.99311, p-value = 0.1103
simulateResiduals(fittedModel = soil_betadisp_LUX_fung_mod, plot = T)


anova(soil_betadisp_LUX_fung_mod, type = 3)
#collectionDate            0.67084 0.047917    14 308.08  2.2575 0.006172 **


emmeans(soil_betadisp_LUX_fung_mod,pairwise~collectionDate)


####Lux Beta dispersion Plot####

#Graphing colors
Date_F <- function(x){
  format(as.Date(x, origin = '1970-01-01',"%Y-%m-%d"),"%b")
}



(lux_bact_betadisp_p3=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,
                                                                              fill=FertStatus,color=FertStatus))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
    scale_fill_manual(values = c("black","grey"), name=NULL)+scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+ylab("Betadispersion")+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(axis.text.x = element_blank()))

(lux_root_bact_betadisp_p2=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis, Root_soil=="Root"), 
                                  aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,
                                      fill=as.Date(collectionDate, format="%m/%d/%Y"),
                                      alpha=FertStatus))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+
    ylab("Bacterial\nBetadispersion")+ggtitle("Root")+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+ylim(c(0.25,0.48))+
    theme(axis.text.x = element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32)))

(lux_soil_bact_betadisp_p2=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis, Root_soil=="Soil"), 
                                  aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,
                                      fill=as.Date(collectionDate, format="%m/%d/%Y"),
                                      alpha=FertStatus))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_alpha_manual(values = c(1,0.4))+ggtitle("Soil")+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+ylim(c(0.25,0.48))+
    theme(axis.text.x = element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.title.y = element_blank(),axis.text.y = element_blank()))



(lux_root_fung_betadisp_p2=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Root"), 
                                  aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,
                                      fill=as.Date(collectionDate, format="%m/%d/%Y"),
                                      alpha=FertStatus))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+ylab("Fungal\nBetadispersion")+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+ylim(c(0.26,0.8))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none"))


(lux_soil_fung_betadisp_p2=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Soil"), 
                                  aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,
                                      fill=as.Date(collectionDate, format="%m/%d/%Y"),
                                      alpha=FertStatus))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_alpha_manual(values = c(1,0.4))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+ylim(c(0.26,0.8))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none",
          axis.title.y = element_blank(),axis.text.y = element_blank()))



(beta_lux_2panel=plot_grid(lux_root_bact_betadisp_p2,lux_soil_bact_betadisp_p2,
                           lux_root_fung_betadisp_p2,lux_soil_fung_betadisp_p2,ncol = 2,labels = c('a)','', 'b)',''), 
                           label_size = 30))

plot_grid(beta_lux_2panel,get_legend(lux_bact_betadisp_p3),ncol = 2,rel_widths = c(4,1))

#NOT INCLUDED IN REPOSITORY
ggsave(plot_grid(beta_lux_2panel,get_legend(lux_bact_betadisp_p3),ncol = 2,rel_widths = c(4,1)),
       filename = "Betadisp_boxplot_Lux_Arbor_date_color_p.png",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height =10)

ggsave(plot_grid(beta_lux_2panel,get_legend(lux_bact_betadisp_p3),ncol = 2,rel_widths = c(4,1)),
       filename = "Betadisp_boxplot_Lux_Arbor_date_color_p.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height =10)
#NOT INCLUDED IN REPOSITORY
