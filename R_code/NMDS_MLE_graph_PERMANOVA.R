
## ---------------------------
##
## Script name: Community characterization of the switchgrass microbiomes of 
## the Marginal Lands Experiment
##
## Purpose of script: Graphing and characterization of the bacteria and fungi
## of roots and soils of switchgrass (Panicum virgatum)
##
## Author: Lukas Bell-Dereske
##
## Email: lukas.dereske@gmail.com
##


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


#####All Sites community composition####


#All Sites is defined as the July overlap in root and soil sampling

unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"&siteID=="LUX"))$collectionDate)
#"5/29/2018" "9/17/2018" "8/20/2018" "7/30/2018" "6/25/2018" "10/3/2018"
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"&siteID=="LC"))$collectionDate)
#"7/10/2018"


GLBRC018_OTU_bact_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,siteID!="LUX"|collectionDate=="7/30/2018")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites)
#513
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LUX"))$collectionDate)

GLBRC018_OTU_bact_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,siteID!="LC"|collectionDate=="7/10/2018")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites)
#323
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LC"))$collectionDate)

#Let's use only G5 since there is overlap in sampling between roots and soil

GLBRC018_OTU_bact_MMPRNT_All_sites_G5=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,plotType=="G5")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
#227
rm(GLBRC018_OTU_bact_MMPRNT_All_sites)
rm(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)

#I am going to define All Sites as the July overlap in root and soil sampling

unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&siteID=="LUX"))$collectionDate)
#"5/29/2018" "9/17/2018" "8/20/2018" "7/30/2018" "6/25/2018" "10/3/2018"
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&siteID=="LC"))$collectionDate)
#"7/10/2018"
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#1086

GLBRC018_OTU_fung_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,siteID!="LUX"|collectionDate=="7/30/2018")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites)
#501
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LUX"))$collectionDate)

GLBRC018_OTU_fung_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,siteID!="LC"|collectionDate=="7/10/2018")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites)
#319
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LC"))$collectionDate)

#Let's use only G5 since there is overlap in sampling between roots and soil

GLBRC018_OTU_fung_MMPRNT_All_sites_G5=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,plotType=="G5")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)
#223
rm(GLBRC018_OTU_fung_MMPRNT_All_sites)
rm(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)



#Separate in organs



#Roots only 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,Root_soil=="Root")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root=
  prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)>0,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)
#113

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,Root_soil=="Root")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root=
  prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)>0,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)
#111

#NMDS
set.seed(2021)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root,method = "NMDS")
#*** Solution reached
#0.1226669   

set.seed(2021)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root,method = "NMDS")
#*** Solution reached
#0.2188881  



#Soil only 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,Root_soil=="Soil")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil=
  prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)>0,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)
#114

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,Root_soil=="Soil")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil=
  prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)>0,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)
#112

#NMDS
set.seed(2021)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil,method = "NMDS")
#*** Solution reached
#0.1199976    

set.seed(2021)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil,method = "NMDS")
#*** Solution reached
#0.1712595   

(all_site_compart_fung_soil_legend=plot_ordination(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_NMDS_points=merge(sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root),GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_ord$points,
                                                             by="row.names")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_NMDS_points_sum=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

(all_site_compart_bact_root_mean_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+ggtitle("Bacteria")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+ylab("NMDS2")+annotate("text", x = -0.55, y = -0.54, label = "Stress = 0.123", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5)),axis.title.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))




GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_NMDS_points=merge(sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root),GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_ord$points,
                                                             by="row.names")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_NMDS_points_sum=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

(all_site_compart_fung_root_mean_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+ggtitle("Fungi")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+annotate("text", x = -0.495, y = -0.77, label = "Stress = 0.219", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5)), axis.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_NMDS_points=merge(sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil),GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_ord$points,
                                                             by="row.names")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_NMDS_points_sum=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

(all_site_compart_bact_soil_mean_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+xlab("NMDS1")+ylab("NMDS2")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+annotate("text", x = -0.51, y = -0.6, label = "Stress = 0.120", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_NMDS_points=merge(sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil),GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_ord$points,
                                                             by="row.names")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_NMDS_points_sum=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))


(all_site_compart_fung_soil_mean_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+xlab("NMDS1")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5))+scale_x_continuous(breaks = c(-1,0,1))+annotate("text", x = -1.02, y = -1.4, label = "Stress = 0.171", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",axis.title.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(root_title=ggplot()+ggtitle("Root")+
    theme(plot.title = (element_text(size = 30,hjust = 0.5, angle = 90))))

(soil_title=ggplot()+ggtitle("Soil")+
    theme(plot.title = (element_text(size = 30,hjust = 0.5, angle = 90))))
y_titles=plot_grid(get_title(root_title),get_title(soil_title),nrow = 2)

(all_site_mean_4_panel=plot_grid(all_site_compart_bact_root_mean_p2,all_site_compart_fung_root_mean_p2,
                                 all_site_compart_bact_soil_mean_p2,all_site_compart_fung_soil_mean_p2,ncol = 2,
                                 label_size = 26,rel_widths = c(1.05,1,1.05,1),
                                 rel_heights = c(1,1,1.05,1.05),labels = c('a)', 'b)','c)', 'd)'),
                                 label_x = c(0,-0.035,0,-0.035),label_y = c(0.92,0.92,0.98,0.98)))

plot_grid(y_titles,all_site_mean_4_panel,get_legend(all_site_compart_fung_soil_legend),ncol = 3,rel_widths = c(.15,4,1))


ggsave(plot_grid(y_titles,all_site_mean_4_panel,get_legend(all_site_compart_fung_soil_legend),ncol = 3,rel_widths = c(.15,4,1)), 
       filename = "NMDS_All_Sites_mean_p.png",path = here::here("Manuscript","MLE_comm_figs"),width = 15,height = 10)


ggsave(plot_grid(y_titles,all_site_mean_4_panel,get_legend(all_site_compart_fung_soil_legend),ncol = 3,rel_widths = c(.15,4,1)), 
       filename = "NMDS_All_Sites_mean_p.svg",path = here::here("Manuscript","MLE_comm_figs"),width = 15,height = 10)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_dis=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,method = "bray")



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_dis=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,method = "bray")


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)


#####All Sites Beta-disp####



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betamod=betadisper(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_dis, with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map,interaction(Root_soil,siteID,FertStatus)))
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis=as.data.frame(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betamod$distances)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis=merge(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map,by="row.names")
colnames(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis)[colnames(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis)==
                                                          "GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betamod$distances"]="betadisp"



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betamod=betadisper(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_dis, with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map,interaction(Root_soil,siteID,FertStatus)))
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis=as.data.frame(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betamod$distances)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis=merge(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map,by="row.names")
colnames(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis)[colnames(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis)==
                                                          "GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betamod$distances"]="betadisp"

head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis)
unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis$siteID)



#Stats 

#Roots Bacteria
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_root=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis,Root_soil=="Root")

root_betadisp_all_site_bact_mod=lmer(log(betadisp)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_root)
plot(root_betadisp_all_site_bact_mod)
hist(resid(root_betadisp_all_site_bact_mod))
qqPlot(resid(root_betadisp_all_site_bact_mod))
shapiro.test(resid(root_betadisp_all_site_bact_mod))
#W = 0.98237, p-value = 0.143
simulateResiduals(fittedModel = root_betadisp_all_site_bact_mod, plot = T)

anova(root_betadisp_all_site_bact_mod)
#siteID            0.37223 0.093057     4 101.18  6.3060 0.0001428 ***
#FertStatus        0.01573 0.015731     1 100.31  1.0660 0.3043382    
#siteID:FertStatus 0.06075 0.015187     4 100.31  1.0292 0.3960468 


emmeans(root_betadisp_all_site_bact_mod,pairwise~siteID)
#ESC - HAN   0.0412 0.0386 103   1.068  0.8222
#ESC - LC   -0.0910 0.0351 100  -2.594  0.0791
#ESC - LUX   0.0166 0.0355 100   0.467  0.9901
#ESC - RHN   0.0770 0.0351 100   2.197  0.1894
#HAN - LC   -0.1322 0.0386 103  -3.427  0.0077
#HAN - LUX  -0.0247 0.0388 103  -0.636  0.9688
#HAN - RHN   0.0358 0.0386 103   0.928  0.8853
#LC - LUX    0.1075 0.0355 100   3.031  0.0252
#LC - RHN    0.1680 0.0351 100   4.791  0.0001
#LUX - RHN   0.0605 0.0355 100   1.705  0.4359

#Soil
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_soil=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis,Root_soil=="Soil")

soil_betadisp_all_site_bact_mod=lmer((betadisp)^-3~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_soil)
plot(soil_betadisp_all_site_bact_mod)
hist(resid(soil_betadisp_all_site_bact_mod))
qqPlot(resid(soil_betadisp_all_site_bact_mod))
shapiro.test(resid(soil_betadisp_all_site_bact_mod))
#W = 0.97723, p-value = 0.04884
simulateResiduals(fittedModel = soil_betadisp_all_site_bact_mod, plot = T)

anova(soil_betadisp_all_site_bact_mod)
#siteID            18877.5  4719.4     4 101.90 28.7874 5.204e-16 ***
#FertStatus          168.8   168.8     1 101.03  1.0295    0.3127    
#siteID:FertStatus   167.3    41.8     4 101.03  0.2552    0.9059  



emmeans(soil_betadisp_all_site_bact_mod,pairwise~siteID)
#ESC - HAN    36.69 4.07 104  9.020  <.0001 
#ESC - LC     32.57 3.70 101  8.812  <.0001 
#ESC - LUX    28.23 3.70 101  7.636  <.0001 
#ESC - RHN    22.24 3.70 101  6.018  <.0001 
#HAN - LC     -4.12 4.07 104 -1.013  0.8487 
#HAN - LUX    -8.47 4.07 104 -2.082  0.2359 
#HAN - RHN   -14.45 4.07 104 -3.553  0.0051 
#LC - LUX     -4.35 3.70 101 -1.176  0.7651 
#LC - RHN    -10.33 3.70 101 -2.795  0.0478 
#LUX - RHN    -5.98 3.70 101 -1.619  0.4889 






#Roots
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_root=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis,Root_soil=="Root")

root_betadisp_all_site_fung_mod=lmer(sqrt(betadisp)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_root)
plot(root_betadisp_all_site_fung_mod)
hist(resid(root_betadisp_all_site_fung_mod))
qqPlot(resid(root_betadisp_all_site_fung_mod))
shapiro.test(resid(root_betadisp_all_site_fung_mod))
#W = 0.98262, p-value = 0.1585
simulateResiduals(fittedModel = root_betadisp_all_site_fung_mod, plot = T)

anova(root_betadisp_all_site_fung_mod)
#siteID            0.078001 0.0195002     4   101  6.5540 9.909e-05 ***
#FertStatus        0.002818 0.0028182     1   101  0.9472    0.3328    
#siteID:FertStatus 0.002614 0.0006535     4   101  0.2197    0.9269  



emmeans(root_betadisp_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate     SE    df t.ratio p.value
#ESC - HAN  0.02244 0.0174 100.8  1.290  0.6982 
#ESC - LC  -0.03331 0.0163  98.5 -2.044  0.2533 
#ESC - LUX  0.04587 0.0159  98.2  2.880  0.0384 
#ESC - RHN  0.00202 0.0159  98.2  0.126  0.9999 
#HAN - LC  -0.05574 0.0177 100.8 -3.150  0.0179 
#HAN - LUX  0.02344 0.0173 101.0  1.354  0.6584 
#HAN - RHN -0.02042 0.0173 101.0 -1.180  0.7630 
#LC - LUX   0.07918 0.0161  98.2  4.915  <.0001 
#LC - RHN   0.03532 0.0161  98.2  2.193  0.1912 
#LUX - RHN -0.04386 0.0157  98.0 -2.785  0.0492 



#Soil
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_soil=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis,Root_soil=="Soil")

soil_betadisp_all_site_fung_mod=lmer(betadisp~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_soil)
plot(soil_betadisp_all_site_fung_mod)
hist(resid(soil_betadisp_all_site_fung_mod))
qqPlot(resid(soil_betadisp_all_site_fung_mod))
shapiro.test(resid(soil_betadisp_all_site_fung_mod))
#W = 0.97726, p-value = 0.05284
simulateResiduals(fittedModel = soil_betadisp_all_site_fung_mod, plot = T)

anova(soil_betadisp_all_site_fung_mod)
#siteID            0.270100 0.067525     4 100.062 17.5100 6.387e-11 ***
#FertStatus        0.005119 0.005119     1  99.196  1.3274    0.2520    
#siteID:FertStatus 0.019024 0.004756     4  99.164  1.2333    0.3017    



emmeans(soil_betadisp_all_site_fung_mod,pairwise~siteID)
#ESC - HAN  -0.1529 0.0197 102.0 -7.748  <.0001 
#ESC - LC   -0.0861 0.0179  99.0 -4.804  0.0001 
#ESC - LUX  -0.1074 0.0181  99.1 -5.922  <.0001 
#ESC - RHN  -0.0668 0.0181  99.1 -3.685  0.0034 
#HAN - LC    0.0667 0.0197 102.0  3.383  0.0088 
#HAN - LUX   0.0455 0.0200 102.0  2.278  0.1605 
#HAN - RHN   0.0860 0.0198 101.7  4.339  0.0003 
#LC - LUX   -0.0213 0.0181  99.1 -1.173  0.7665 
#LC - RHN    0.0193 0.0181  99.1  1.063  0.8249 
#LUX - RHN   0.0406 0.0183  99.3  2.210  0.1844


#All Sites Beta-disp figures####
site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")


brewer.pal(5,"Set1")
site_colors=c("ESC"="#E41A1C", "HAN"="#377EB8", "LC"="#4DAF4A", 
              "LUX"="#984EA3", "RHN"="#FF7F00")

#Betadisp

root_bact_disp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("ab","a","b","a","a"),
                                  Root_soil=rep("Root"))

soil_bact_disp_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("AB","A","C","AB","B"),
                                  Root_soil=rep("Soil"))

root_soil_bact_disp_letters=rbind(root_bact_disp_letters,soil_bact_disp_letters)
bact_betadisp_max=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(betadisp))

bact_betadisp_max_disp_letters=merge(root_soil_bact_disp_letters,bact_betadisp_max,by=c("siteID","Root_soil"))




root_fung_disp_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("a","b","bc","ac","c"),
                                  Root_soil=rep("Root"))


soil_fung_disp_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("AB","A","C","B","A"),
                                  Root_soil=rep("Soil"))
root_soil_fung_disp_letters=rbind(root_fung_disp_letters,soil_fung_disp_letters)

fung_betadisp_max=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(betadisp))

fung_betadisp_max_disp_letters=merge(root_soil_fung_disp_letters,fung_betadisp_max,by=c("siteID","Root_soil"))







(bact_betadisp_p3=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),
                     fill=FertStatus,color=FertStatus))+
    geom_text(data = bact_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_fill_manual(values = c("black","grey"), name=NULL)+
    scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    scale_y_continuous(name = "Betadispersion",limits = c(0.22,0.57))+theme(axis.text.x = element_blank()))





(bact_betadisp_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = bact_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    scale_y_continuous(name = "Bacterial\nBetadispersion",limits = c(0.22,0.57))+
    theme(axis.text.x = element_blank(),legend.position = "none",
          strip.background = element_rect(fill="white",color = "black")))

(fung_betadisp_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = fung_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Fungal\nBetadispersion",limits = c(0.2,0.77))+ 
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none"))


(betadisp_2panel=plot_grid(bact_betadisp_p2,fung_betadisp_p2,ncol = 1,labels = c('a)', 'b)'), label_size = 30))


ggsave(plot_grid(betadisp_2panel,get_legend(bact_betadisp_p3),ncol = 2,rel_widths = c(4,1)),
       filename = "Betadisp_boxplot_All_Sites_p.png",path = here::here("Manuscript","MLE_comm_figs"),width = 20,height =10)

ggsave(plot_grid(betadisp_2panel,get_legend(bact_betadisp_p3),ncol = 2,rel_widths = c(4,1)),
       filename = "Betadisp_boxplot_All_Sites_p.svg",path = here::here("Manuscript","MLE_comm_figs"),width = 20,height =10)




