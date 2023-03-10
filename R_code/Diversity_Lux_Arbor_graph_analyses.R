## ---------------------------
##
## Script name: Community richness and diversity of the switchgrass microbiomes across 
## the Lux Arbor Growing Season
##
## Purpose of script: Estimating and graphing the richness and diversity of the bacteria and fungi
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

####Diversity ####

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div=estimate_richness(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div)
#495
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,by="row.names")


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div=estimate_richness(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map,by="row.names")


#####Richness####


#Stats

#Roots
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_rich_LUX_bact_mod=lmer(Observed~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root)
plot(root_rich_LUX_bact_mod)
hist(resid(root_rich_LUX_bact_mod))
qqPlot(resid(root_rich_LUX_bact_mod))
shapiro.test(resid(root_rich_LUX_bact_mod))
#W = 0.99436, p-value = 0.8554
simulateResiduals(fittedModel = root_rich_LUX_bact_mod, plot = T)

anova(root_rich_LUX_bact_mod)
#collectionDate            524869  104974     5 127.05  5.4371 0.0001444 ***
#FertStatus                 46751   46751     1 127.04  2.4215 0.1221689    
#collectionDate:FertStatus 120034   24007     5 127.05  1.2434 0.2927484    







#Soils
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_rich_LUX_bact_mod=lmer((Observed)^2~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_rich_LUX_bact_mod)
hist(resid(soil_rich_LUX_bact_mod))
qqPlot(resid(soil_rich_LUX_bact_mod))
shapiro.test(resid(soil_rich_LUX_bact_mod))
#W = 0.997, p-value = 0.7637
simulateResiduals(fittedModel = soil_rich_LUX_bact_mod, plot = T)


anova(soil_rich_LUX_bact_mod)
#collectionDate            9.2361e+13 6.5972e+12    14 3.3998e+10 17.6458 < 2e-16 ***
#FertStatus                1.5374e+12 1.5374e+12     1 8.6013e+11  4.1121 0.04258 *  
#collectionDate:FertStatus 7.6441e+12 5.4601e+11    14 9.5095e+10  1.4604 0.11670   

emmeans(soil_rich_LUX_bact_mod,pairwise~FertStatus|collectionDate)
#collectionDate = 8/20/2018:
#contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert   782623 260722 320  3.002  0.0029 
#collectionDate = 4/30/2018:
# contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert  -464653 255285 320 -1.820  0.0697 




#Roots
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_rich_LUX_mod=lmer(Observed~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root)
plot(root_rich_LUX_mod)
hist(resid(root_rich_LUX_mod))
qqPlot(resid(root_rich_LUX_mod))
shapiro.test(resid(root_rich_LUX_mod))
#W = 0.99478, p-value = 0.8941
simulateResiduals(fittedModel = root_rich_LUX_mod, plot = T)

anova(root_rich_LUX_mod)
#collectionDate            134575 26915.1     5 126.02 23.5489 < 2e-16 ***
#FertStatus                  4088  4087.9     1 126.03  3.5767 0.06089 .  
#collectionDate:FertStatus   5012  1002.4     5 126.02  0.8771 0.49861   

emmeans(root_rich_LUX_mod,pairwise~FertStatus|collectionDate)
#collectionDate = 10/3/2018:
#contrast      estimate   SE  df t.ratio p.value
#Fert - Unfert   -31.83 13.8 126 -2.306  0.0227  


#Soils
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_rich_LUX_mod=lmer((Observed)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_rich_LUX_mod)
hist(resid(soil_rich_LUX_mod))
qqPlot(resid(soil_rich_LUX_mod))
shapiro.test(resid(soil_rich_LUX_mod))
#W = 0.99539, p-value = 0.3919
simulateResiduals(fittedModel = soil_rich_LUX_mod, plot = T)


anova(soil_rich_LUX_mod, type = 3)
#collectionDate            175880 12562.9    14 316.01  8.5704 9.729e-16 ***
#FertStatus                   731   730.9     1 316.01  0.4986    0.4806    
#collectionDate:FertStatus  21023  1501.7    14 316.01  1.0244    0.4283  




####Inverse Simpson####



#Stats

#Roots
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_invsimp_LUX_bact_mod=lmer(sqrt(InvSimpson)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root)
plot(root_invsimp_LUX_bact_mod)
hist(resid(root_invsimp_LUX_bact_mod))
qqPlot(resid(root_invsimp_LUX_bact_mod))
shapiro.test(resid(root_invsimp_LUX_bact_mod))
#W = 0.99026, p-value = 0.43
simulateResiduals(fittedModel = root_invsimp_LUX_bact_mod, plot = T)


anova(root_invsimp_LUX_bact_mod)
#collectionDate            84.356 16.8711     5 127.01  9.1702 1.802e-07 ***
#FertStatus                 4.706  4.7060     1 127.01  2.5579    0.1122    
#collectionDate:FertStatus  4.996  0.9993     5 127.01  0.5432    0.7433   





#Soils
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_invsimp_LUX_bact_mod=lmer((InvSimpson)^2~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_invsimp_LUX_bact_mod)
hist(resid(soil_invsimp_LUX_bact_mod))
qqPlot(resid(soil_invsimp_LUX_bact_mod))
shapiro.test(resid(soil_invsimp_LUX_bact_mod))
#W = 0.99024, p-value = 0.0191
simulateResiduals(fittedModel = soil_invsimp_LUX_bact_mod, plot = T)

anova(soil_invsimp_LUX_bact_mod)
#collectionDate            8.0301e+10 5735804522    14 320.01  5.1549 9.056e-09 ***
#FertStatus                4.3726e+09 4372598865     1 320.00  3.9298   0.04829 *   
#collectionDate:FertStatus 1.3600e+10  971425975    14 320.00  0.8730   0.58878

emmeans(soil_invsimp_LUX_bact_mod,pairwise~FertStatus|collectionDate)

emmeans(soil_invsimp_LUX_bact_mod,pairwise~FertStatus)

#Roots
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_invsimp_LUX_mod=lmer(log(InvSimpson)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root)
plot(root_invsimp_LUX_mod)
hist(resid(root_invsimp_LUX_mod))
qqPlot(resid(root_invsimp_LUX_mod))
shapiro.test(resid(root_invsimp_LUX_mod))
#W = 0.98585, p-value = 0.1569
simulateResiduals(fittedModel = root_invsimp_LUX_mod, plot = T)


anova(root_invsimp_LUX_mod)
#collectionDate            6.8094 1.36188     5 126.04  7.9468 1.548e-06 ***
#FertStatus                1.8709 1.87090     1 126.06 10.9170  0.001241 ** 
#collectionDate:FertStatus 0.5566 0.11131     5 126.06  0.6495  0.662365   

emmeans(root_invsimp_LUX_mod,pairwise~FertStatus|collectionDate)
#$contrasts
#collectionDate = 6/25/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.4082 0.169 126 -2.416  0.0171 

#collectionDate = 10/3/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.3184 0.169 126 -1.884  0.0619 


#Soils
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_invsimp_LUX_mod=lmer((InvSimpson)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_invsimp_LUX_mod)
hist(resid(soil_invsimp_LUX_mod))
qqPlot(resid(soil_invsimp_LUX_mod))
shapiro.test(resid(soil_invsimp_LUX_mod))
#W = 0.99593, p-value = 0.5087
simulateResiduals(fittedModel = soil_invsimp_LUX_mod, plot = T)

anova(soil_invsimp_LUX_mod)
#collectionDate            2914.02 208.144    14 316.02  1.6894 0.05651 .
#FertStatus                  82.39  82.393     1 316.02  0.6688 0.41410  
#collectionDate:FertStatus 1674.90 119.636    14 316.02  0.9710 0.48289  



#Lux Arbor Diversity plot####


pairwise_fert_bact_rich_lux=data.frame(y_bot=c(0,0,2310,2500),
                                       x_min=as.Date(c(NA,NA,"4/25/2018","8/15/2018"), format="%m/%d/%Y"),
                                       x_max=as.Date(c(NA,NA,"5/5/2018","8/25/2018"), format="%m/%d/%Y"),
                                       annot_text=c("","","#","**"),
                                       Root_soil=c(rep("Root",2),rep("Soil",2)))

pairwise_fert_fung_invsimp_lux=data.frame(y_bot=c(27,18.5),
                                          x_min=as.Date(c("6/20/2018","9/29/2018"), format="%m/%d/%Y"),
                                          x_max=as.Date(c("6/30/2018","10/8/2018"), format="%m/%d/%Y"),
                                          annot_text=c("*","#"),
                                          Root_soil=c(rep("Root",2)))

(lux_bact_invsimp_p3.1=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=InvSimpson,
                                                                           color=FertStatus,fill=FertStatus))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
    facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+ylab("Inverse Simpson")+
    scale_fill_manual(values = c("black","grey"), name=NULL)+scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+
    theme(axis.text.x = element_blank(),legend.key.size = unit(2, 'cm'),legend.text = element_text(size=28)))



Date_F <- function(x){
  format(as.Date(x, origin = '1970-01-01',"%Y-%m-%d"),"%b")
}
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div$Observed)
(lux_bact_root_rich_p2=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Root"), 
                              aes(y=Observed))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     group=interaction(collectionDate,FertStatus),
                     alpha=FertStatus),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+
    scale_y_continuous(name="Bacterial\nRichness",limits = c(930,2780))+ggtitle("Root")+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.text.x = element_blank()))

(lux_bact_soil_rich_p2=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil"))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=Observed,
                     fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     alpha=FertStatus,group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_y_continuous(limits = c(930,2780))+
    scale_alpha_manual(values = c(1,0.4))+ggtitle("Soil")+
    geom_signif(data=pairwise_fert_bact_rich_lux,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.title.y = element_blank(),axis.text = element_blank()))

plot_grid(lux_bact_root_rich_p2,lux_bact_soil_rich_p2)

summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div$Observed)
(lux_fung_root_rich_p2=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Root"), 
                              aes(y=Observed))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     group=interaction(collectionDate,FertStatus),
                     alpha=FertStatus),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+
    scale_y_continuous(name="Fungal\nRichness",limits = c(70,520))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.text.x = element_blank()))

(lux_fung_soil_rich_p2=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil"))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=Observed,
                     fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     alpha=FertStatus,group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_y_continuous(name="Fungal\nRichness",limits = c(70,520))+
    scale_alpha_manual(values = c(1,0.4))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.title.y = element_blank(),axis.text = element_blank()))


(lux_bact_root_invSimp_p2=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Root"), 
                                 aes(y=InvSimpson))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     group=interaction(collectionDate,FertStatus),
                     alpha=FertStatus),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+
    scale_y_continuous(name="Bacterial\nInverse Simpson",limits = c(8,465))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.text.x = element_blank()))

(lux_bact_soil_invSimp_p2=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil"))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=InvSimpson,
                     fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     alpha=FertStatus,group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_y_continuous(limits = c(8,465))+
    scale_alpha_manual(values = c(1,0.4))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.title.y = element_blank(),axis.text = element_blank()))


(lux_fung_root_invSimp_p2=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Root"), 
                                 aes(y=InvSimpson))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     group=interaction(collectionDate,FertStatus),
                     alpha=FertStatus),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+
    geom_signif(data=pairwise_fert_fung_invsimp_lux,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    scale_y_continuous(name="Fungi\nInverse Simpson",limits = c(1,72))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32)))

(lux_fung_soil_invSimp_p2=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil"))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=InvSimpson,
                     fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     alpha=FertStatus,group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_y_continuous(limits = c(1,72))+
    scale_alpha_manual(values = c(1,0.4))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.title.y = element_blank(),axis.text.y = element_blank()))

plot_grid(lux_fung_root_invSimp_p2,lux_fung_soil_invSimp_p2)


#
(rich_invSimp_LUX_4panel=plot_grid(plot_grid(lux_bact_root_rich_p2,lux_bact_soil_rich_p2),
                                   plot_grid(lux_fung_root_rich_p2,lux_fung_soil_rich_p2),
                                   plot_grid(lux_bact_root_invSimp_p2,lux_bact_soil_invSimp_p2),
                                   plot_grid(lux_fung_root_invSimp_p2,lux_fung_soil_invSimp_p2),ncol = 1,
                                   labels = c('a)', 'b)','c)','d)'), label_size = 30,align = "v",
                                   rel_heights = c(1,0.8,0.8,1)))

plot_grid(rich_invSimp_LUX_4panel,get_legend(lux_bact_invsimp_p3.1),ncol = 2,rel_widths = c(4,1.5))


#NOT INCLUDED IN REPOSITORY
ggsave(plot_grid(rich_invSimp_LUX_4panel,get_legend(lux_bact_invsimp_p3.1),ncol = 2,rel_widths = c(4,1)),
       filename = "richness_invSimpson_boxplot_Lux_Arbor.png",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 15,height =15)

ggsave(plot_grid(rich_invSimp_LUX_4panel,get_legend(lux_bact_invsimp_p3.1),ncol = 2,rel_widths = c(4,1.5)),
       filename = "richness_invSimpson_boxplot_Lux_Arbor.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 15,height =15)
#NOT INCLUDED IN REPOSITORY



