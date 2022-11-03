## ---------------------------
##
## Script name: Core Taxa from switchgrass microbiomes of the Marginal Lands Experiment
##
## Purpose of script: Identification and graphing of Core taxa from the bacterial and fungal communities 
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
library(vegan)
library(metacoder)

# Fitting neutral models through the sncm.fit.R function modified form 
# Burns et al. - 2016 ISME JOURNAL
source(here::here("R_functions","FitNeutral.R"))
# Example: ExtractCore(phyloseq_object, "Date", "elbow", Treatment", "Fertilized", "")

#ploting the neutral model from FitNeutral
source(here::here("R_functions","PlotNeutral.R"))

# Below function will filter the dataset by Fertilized samples only and calculate
# the core according to the different sampling dates.
#method is "total", "elbow", or number which sets the Bray-Curtis cut of a percent (whole number)
source(here::here("R_functions","ExtractCoreFlex.R"))

#Graphing functions for Bray-Curtis thresholds for diagnostic 
#PlotBCincreaseFlex()
#PlotBCThreshold_Rich()
#PlotBCThreshold_Abun()
source(here::here("R_functions","PlotBCincrease_Richness_Abundance.R"))

#Load CONSTAX created database
GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus = read.delim(here::here("Publish_data",
                                                                     "constax_taxonomy.txt"),
                                                          sep = "\t",header = T,fill=T, row.names = 1)

GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus[GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus==""]="Unknown"
#Load community data and make the phyloseq object


#load FunGuild Database based on CONSTAX2 database


GLBRC_CON_FunGuild_raw_F=read.csv(here::here("Publish_data","GLBRC_CON_FunGuild_raw_F.csv"),
                                  header = T)

dim(GLBRC_CON_FunGuild_raw_F)
#2432   18


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


#####Fitting Core Taxa Bacteria Roots ####

##Prioritizing core microbiome based on the SITE

#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5=ExtractCoreFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root, "siteID", 5, Group=NULL, Level=NULL)

#Plotting the effects of BC thresholds on OTU inclusion
PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,400,5)

#Fit the neutral model to the core taxa
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5_Neut=FitNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5)

PlotNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5_Neut)


#####Fitting Core Taxa Bacteria Soil All Sites####

##Prioritizing core microbiome based on the SITE


#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5=ExtractCoreFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil, "siteID", 5, Group=NULL, Level=NULL)


#Ploting the effects of BC thresholds on OTU incluesion
PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,400,5)


#Fit the nuetral model to the core taxa
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5_Neut=FitNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5)

PlotNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5_Neut)

plot_grid(PlotNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5_Neut),
          PlotNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5_Neut),
          PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,1000,5),
          PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,1000,5),ncol = 2,
          labels = c("a)","b)","c)","d)"))

plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,20,5),
          PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,20,5),ncol = 2,
          labels = c("a)","b)","c)","d)"))

ggsave(plot_grid(PlotNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5_Neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,1000,5),
                 PlotNeutral(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5_Neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,1000,5),ncol = 2,
                 labels = c("a)","b)","c)","d)")), filename = "All_sites_bacteria_core_fits.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)



ggsave(plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5,20,5),
                 PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5,20,5),ncol = 2,
                 labels = c("a)","b)","c)","d)")), filename = "All_sites_bacteria_core_diagnostic.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)





#####Fungi Root All Sites Core Taxa####

##Prioritizing core microbiome based on the SITE

#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5=ExtractCoreFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root, "siteID", 5, Group=NULL, Level=NULL)

#Ploting the effects of BC thresholds on OTU incluesion
PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,1000,5)

#Fit the nuetral model to the core taxa
FitNeutral(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5)
#Waiting for profiling to be done...
#Error in optim(start, f, method = method, hessian = TRUE, ...) : 
#  non-finite finite-difference value [1]

ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],aes(x=log(otu_rel), y=otu_occ,color=fill ))+labs(x="Log10(mean abundance)", y="Occupancy") + 
  geom_point()+theme_cowplot()

#####Fungi Soil All Sites Core Taxa####

##Prioritizing core microbiome based on the SITE


#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5=ExtractCoreFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil, "siteID", 5, Group=NULL, Level=NULL)


#Ploting the effects of BC thresholds on OTU incluesion
PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,1000,5)

#Fit the nuetral model to the core taxa
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5_Neut=FitNeutral(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5)

PlotNeutral(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5_Neut)




plot_grid(ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],aes(x=log(otu_rel), y=otu_occ,color=fill ))+
            labs(x="Log10(mean abundance)", y="Occupancy") + 
            geom_point()+theme_cowplot(),
          PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,1000,5),ncol = 2,
          PlotNeutral(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5_Neut),
          PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,1000,5),
          labels = c("a)","b)","c)","d)"))





plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,20,5),
          PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,20,5),ncol = 2,
          labels = c("a)","b)","c)","d)"))


ggsave(plot_grid(ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],aes(x=log(otu_rel), y=otu_occ,color=fill ))+
                   labs(x="Log10(mean abundance)", y="Occupancy") + 
                   geom_point()+theme_cowplot(),
                 PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,1000,5),ncol = 2,
                 PlotNeutral(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5_Neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,1000,5),
                 labels = c("a)","b)","c)","d)")), filename = "All_sites_fungi_core_fit.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)


ggsave(plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5,20,5),
                 PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5,20,5),ncol = 2,
                 labels = c("a)","b)","c)","d)")), filename = "All_sites_fungi_core_diagnostic.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)


#####Stack Soil Bacteria Core Taxa plot####


#Need to turn the subsetted community into a phyloseq obj
subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5[[4]],fill=="core")$otu
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5[[4]],fill=="core")$otu,
                                                                GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)

ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl)
#133
ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl)/ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)
#0.005701548

sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl))
#342556
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil))
#0.3004877

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl))
#4077
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl))
#2230
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl)
#114

#Calculated the relative abundance of the core community
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl_map=(sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl))
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl)
bact_All_sites_soil_core_freq=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl_map%>%group_by(siteID)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)

#merge OTUs by site
GLBRC018_OTU_bact_soil_core_facet=merge_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_phyl, "siteID")
sample_names(GLBRC018_OTU_bact_soil_core_facet)     

#Order Level
#combine the reads at Class level
get_taxa_unique(GLBRC018_OTU_bact_soil_core_facet, taxonomic.rank="Order")
#[1] "o:Rhizobiales"         "o:Sphingomonadales"    "o:Streptomycetales"    "o:Corynebacteriales"   "o:Caulobacterales"     "o:Pseudomonadales"     "o:Burkholderiales"     "o:Anaerolineales"      "o:Subgroup_7"          "o:Micromonosporales"   "o:Micrococcales"    [12] "o:Myxococcales"        "o:Pseudonocardiales"   "o:Acidimicrobiales"    "o:Bacillales"          "Unknown"               "o:Opitutales"          "o:Gemmatimonadales"    "o:Nitrosomonadales"    "o:Solirubrobacterales" "o:Subgroup_6"          "o:Xanthomonadales"  [23] "o:Frankiales"          "o:Chthoniobacterales"  "o:SC-I-84"             "o:Sphingobacteriales"  "o:Subgroup_3"          "o:Oligoflexales"       "o:Planctomycetales"    "o:Gaiellales"          "o:Acidobacteriales"    "o:Subgroup_4"          "o:TRA3-20"        [34] "o:Rhodospirillales"    "o:Subgroup_17"         "o:Lineage_IIa"         "o:Nitrospirales"       "o:Subgroup_5"          "o:WD2101_soil_group" 
(GLBRC018_OTU_bact_soil_core_facet.order<-tax_glom(GLBRC018_OTU_bact_soil_core_facet, taxrank="Order"))
#46



GLBRC018_OTU_bact_soil_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_bact_soil_core_facet.order))$Phylum,
                                                             "Order"=data.frame(tax_table(GLBRC018_OTU_bact_soil_core_facet.order))$Order,
                                                             "OTU"=taxa_names(GLBRC018_OTU_bact_soil_core_facet.order))

GLBRC018_OTU_bact_soil_core_facet.order_Names_raw$P_O=with(GLBRC018_OTU_bact_soil_core_facet.order_Names_raw,interaction(Phylum,Order))
#Let's take the top 10 taxa 

AllSite_soil_bact_TopORDER = names(sort(taxa_sums(GLBRC018_OTU_bact_soil_core_facet.order), TRUE)[1:10])
AllSite_soil_bact_TOP10_Order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(AllSite_soil_bact_TopORDER, GLBRC018_OTU_bact_soil_core_facet.order)))$Phylum,
                                               "Order"=data.frame(tax_table(prune_taxa(AllSite_soil_bact_TopORDER, GLBRC018_OTU_bact_soil_core_facet.order)))$Order,
                                               "OTU"=taxa_names(prune_taxa(AllSite_soil_bact_TopORDER, GLBRC018_OTU_bact_soil_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_bact_soil_core_facet.order_Names=merge(GLBRC018_OTU_bact_soil_core_facet.order_Names_raw,AllSite_soil_bact_TOP10_Order_Names, by="OTU", all.x = T)

GLBRC018_OTU_bact_soil_core_facet.order_Names$Real_P_O=ifelse(is.na(GLBRC018_OTU_bact_soil_core_facet.order_Names$Order.y),
                                                              paste(GLBRC018_OTU_bact_soil_core_facet.order_Names$Phylum.x,
                                                                    "Other_spp",sep="."),
                                                              ifelse(GLBRC018_OTU_bact_soil_core_facet.order_Names$Order.x=="Unknown",
                                                                     paste(GLBRC018_OTU_bact_soil_core_facet.order_Names$Phylum.x,
                                                                           "Other_spp",sep="."),                              
                                                                     as.character(GLBRC018_OTU_bact_soil_core_facet.order_Names$P_O)
                                                              ))




#Transform the read counts to prop of total reads

GLBRC018_OTU_bact_soil_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_bact_soil_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_bact_soil_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_bact_soil_core_facet.order.prop)))
GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2=merge(GLBRC018_OTU_bact_soil_core_facet.order_Names[,c("Real_P_O","OTU")],
                                                        GLBRC018_OTU_bact_soil_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2$OTU=NULL


GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_trunc=
  GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2[order(-rowSums(
    GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2)])),][1:10,]

GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_trunc1=GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))

GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M=pivot_longer(data.frame(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_trunc1),!Real_P_O, 
                                                                     names_to = "variable")
sort(unique(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M$Real_P_O))

GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M2=merge(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M,
                                                               bact_All_sites_soil_core_freq,by.x="variable",
                                                               by.y="siteID")
GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M2=
  GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M2%>%mutate(variable_std=value*samp_sum_mean)
sort(unique(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M2$Real_P_O))

All_site_bact_soil_C_order_order=c("p:Acidobacteria.o:Subgroup_4","p:Acidobacteria.o:Subgroup_6",
                                   "p:Actinobacteria.o:Micrococcales","p:Actinobacteria.o:Solirubrobacterales","p:Bacteroidetes.o:Sphingobacteriales",
                                   "p:Firmicutes.o:Bacillales","p:Proteobacteria.o:Pseudomonadales","p:Proteobacteria.o:Rhizobiales",
                                   "p:Proteobacteria.o:Sphingomonadales","p:Verrucomicrobia.o:Chthoniobacterales","Other_Taxa")

#Names



All_site_bact_soil_root_C_order_names=c("p:Acidobacteria.o:Subgroup_4"="Subgroup_4",
                                        "p:Acidobacteria.o:Subgroup_6"="Subgroup_6",
                                        "p:Actinobacteria.o:Kineosporiales"="Kineosporiales",
                                        "p:Actinobacteria.o:Micrococcales"="Micrococcales",
                                        "p:Actinobacteria.o:Micromonosporales"="Micromonosporales",
                                        "p:Actinobacteria.o:Pseudonocardiales"="Pseudonocardiales",
                                        "p:Actinobacteria.o:Solirubrobacterales"="Solirubrobacterales",
                                        "p:Bacteroidetes.o:Sphingobacteriales"="Sphingobacteriales",
                                        "p:Firmicutes.o:Bacillales"="Bacillales",
                                        "p:Proteobacteria.o:Burkholderiales"="Burkholderiales",
                                        "p:Proteobacteria.o:Caulobacterales"="Caulobacterales",
                                        "p:Proteobacteria.o:Myxococcales"="Myxococcales",
                                        "p:Proteobacteria.o:Pseudomonadales"="Pseudomonadales",
                                        "p:Proteobacteria.o:Rhizobiales"="Rhizobiales",
                                        "p:Proteobacteria.o:Sphingomonadales"="Sphingomonadales",
                                        "p:Proteobacteria.o:Xanthomonadales"="Xanthomonadales",
                                        "p:Verrucomicrobia.o:Chthoniobacterales"="Chthoniobacterales",
                                        "Other_Taxa"="Other taxa"
)






All_site_bact_soil_root_C_order_color=c("p:Acidobacteria.o:Subgroup_4"="#1B9E77",
                                        "p:Acidobacteria.o:Subgroup_6"="#147759",
                                        "p:Actinobacteria.o:Kineosporiales"="#D95F02",
                                        "p:Actinobacteria.o:Micrococcales"="#A34702",
                                        "p:Actinobacteria.o:Micromonosporales"="#6D3001",
                                        "p:Actinobacteria.o:Pseudonocardiales"="#361801",
                                        "p:Actinobacteria.o:Solirubrobacterales"="#1B0C00",
                                        "p:Bacteroidetes.o:Sphingobacteriales"="#7570B3",
                                        "p:Firmicutes.o:Bacillales"="#E7298A",
                                        "p:Proteobacteria.o:Burkholderiales"="#66A61E",
                                        "p:Proteobacteria.o:Caulobacterales"="#59911A",
                                        "p:Proteobacteria.o:Myxococcales"="#4D7D16",
                                        "p:Proteobacteria.o:Pseudomonadales"="#406813",
                                        "p:Proteobacteria.o:Rhizobiales"="#33530F",
                                        "p:Proteobacteria.o:Sphingomonadales"="#263E0B",
                                        "p:Proteobacteria.o:Xanthomonadales"="#1A2A07",
                                        "p:Verrucomicrobia.o:Chthoniobacterales"="#E6AB02",
                                        "Other_Taxa"="grey"
)


site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")

(o_All_Sites__bact_soil_color=ggplot(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M2,aes(x=factor(variable,levels = site_order),
                                                                                                  y=variable_std))+
    geom_bar(stat = "identity",aes( fill=factor(Real_P_O,levels = All_site_bact_soil_C_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_fill_manual(values = All_site_bact_soil_root_C_order_color,labels=All_site_bact_soil_root_C_order_names)+
    guides(fill=guide_legend(title="Order")))

#1300x650





#####Stack Root Bacteria Core Taxa plot####

#Need to turn the subsetted community into a phyloseq obj

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5[[4]],fill=="core")$otu,
                                                                GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)

ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl)
#33
ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl)/ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)
#0.002315789
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl))
#387474
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root))
#0.3428973

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl))
#5207
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl))
#1879
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl)
#113


#Calculated the relative abundance of the core community
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl_map=(sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl))
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl)
bact_All_sites_root_core_freq=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl_map%>%group_by(siteID)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)

#merge OTUs by site
GLBRC018_OTU_bact_root_core_facet=merge_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_phyl, "siteID")
sample_names(GLBRC018_OTU_bact_root_core_facet)     


#Order Level
#combine the reads at Class level
get_taxa_unique(GLBRC018_OTU_bact_root_core_facet, taxonomic.rank="Order")
#[1] "o:Burkholderiales"     "o:Sphingomonadales"    "o:Streptomycetales"    "o:Caulobacterales"     "o:Corynebacteriales"   "o:Pseudonocardiales"  
#[7] "o:Rhizobiales"         "o:Pseudomonadales"     "o:Micromonosporales"   "o:Bacillales"          "o:Solirubrobacterales" "o:Micrococcales"      
#[13] "o:Kineosporiales"      "o:Myxococcales"        "o:Xanthomonadales"
(GLBRC018_OTU_bact_root_core_facet.order<-tax_glom(GLBRC018_OTU_bact_root_core_facet, taxrank="Order"))
#15



GLBRC018_OTU_bact_root_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_bact_root_core_facet.order))$Phylum,
                                                             "Order"=data.frame(tax_table(GLBRC018_OTU_bact_root_core_facet.order))$Order,
                                                             "OTU"=taxa_names(GLBRC018_OTU_bact_root_core_facet.order))

GLBRC018_OTU_bact_root_core_facet.order_Names_raw$P_O=with(GLBRC018_OTU_bact_root_core_facet.order_Names_raw,interaction(Phylum,Order))
#Let's take the top 10 taxa 

AllSite_root_bact_TopORDER = names(sort(taxa_sums(GLBRC018_OTU_bact_root_core_facet.order), TRUE)[1:10])
AllSite_root_bact_TOP10_Order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(AllSite_root_bact_TopORDER, GLBRC018_OTU_bact_root_core_facet.order)))$Phylum,
                                               "Order"=data.frame(tax_table(prune_taxa(AllSite_root_bact_TopORDER, GLBRC018_OTU_bact_root_core_facet.order)))$Order,
                                               "OTU"=taxa_names(prune_taxa(AllSite_root_bact_TopORDER, GLBRC018_OTU_bact_root_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_bact_root_core_facet.order_Names=merge(GLBRC018_OTU_bact_root_core_facet.order_Names_raw,AllSite_root_bact_TOP10_Order_Names, by="OTU", all.x = T)

GLBRC018_OTU_bact_root_core_facet.order_Names$Real_P_O=ifelse(is.na(GLBRC018_OTU_bact_root_core_facet.order_Names$Order.y),
                                                              paste(GLBRC018_OTU_bact_root_core_facet.order_Names$Phylum.x,
                                                                    "Other_spp",sep="."),
                                                              ifelse(GLBRC018_OTU_bact_root_core_facet.order_Names$Order.x=="Unknown",
                                                                     paste(GLBRC018_OTU_bact_root_core_facet.order_Names$Phylum.x,
                                                                           "Other_spp",sep="."),                              
                                                                     as.character(GLBRC018_OTU_bact_root_core_facet.order_Names$P_O)
                                                              ))

#Transform the read counts to prop of total reads

GLBRC018_OTU_bact_root_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_bact_root_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_bact_root_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_bact_root_core_facet.order.prop)))
GLBRC018_OTU_bact_root_core_facet.order.prop_otu2=merge(GLBRC018_OTU_bact_root_core_facet.order_Names[,c("Real_P_O","OTU")],
                                                        GLBRC018_OTU_bact_root_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

GLBRC018_OTU_bact_root_core_facet.order.prop_otu2$OTU=NULL

rowSums(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2)])
GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_trunc=
  GLBRC018_OTU_bact_root_core_facet.order.prop_otu2[order(-rowSums(
    GLBRC018_OTU_bact_root_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2)])),][1:10,]

GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_trunc1=GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))

GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M=pivot_longer(data.frame(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_trunc1),
                                                                     !Real_P_O, names_to = "variable")
sort(unique(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M$Real_P_O))

GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M2=merge(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M,
                                                               bact_All_sites_root_core_freq,by.x="variable",
                                                               by.y="siteID")
GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M2=
  GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M2%>%mutate(variable_std=value*samp_sum_mean)
sort(unique(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M2$Real_P_O))

All_site_bact_root_C_order_order=c("p:Actinobacteria.o:Kineosporiales","p:Actinobacteria.o:Micromonosporales","p:Actinobacteria.o:Pseudonocardiales",
                                   "p:Proteobacteria.o:Burkholderiales","p:Proteobacteria.o:Caulobacterales","p:Proteobacteria.o:Myxococcales",
                                   "p:Proteobacteria.o:Pseudomonadales","p:Proteobacteria.o:Rhizobiales","p:Proteobacteria.o:Sphingomonadales",
                                   "p:Proteobacteria.o:Xanthomonadales","Other_Taxa")
length(All_site_bact_root_C_order_order)

#Names



All_site_bact_soil_root_C_order_names=c("p:Acidobacteria.o:Subgroup_4"="Subgroup_4",
                                        "p:Acidobacteria.o:Subgroup_6"="Subgroup_6",
                                        "p:Actinobacteria.o:Kineosporiales"="Kineosporiales",
                                        "p:Actinobacteria.o:Micrococcales"="Micrococcales",
                                        "p:Actinobacteria.o:Micromonosporales"="Micromonosporales",
                                        "p:Actinobacteria.o:Pseudonocardiales"="Pseudonocardiales",
                                        "p:Actinobacteria.o:Solirubrobacterales"="Solirubrobacterales",
                                        "p:Bacteroidetes.o:Sphingobacteriales"="Sphingobacteriales",
                                        "p:Firmicutes.o:Bacillales"="Bacillales",
                                        "p:Proteobacteria.o:Burkholderiales"="Burkholderiales",
                                        "p:Proteobacteria.o:Caulobacterales"="Caulobacterales",
                                        "p:Proteobacteria.o:Myxococcales"="Myxococcales",
                                        "p:Proteobacteria.o:Pseudomonadales"="Pseudomonadales",
                                        "p:Proteobacteria.o:Rhizobiales"="Rhizobiales",
                                        "p:Proteobacteria.o:Sphingomonadales"="Sphingomonadales",
                                        "p:Proteobacteria.o:Xanthomonadales"="Xanthomonadales",
                                        "p:Verrucomicrobia.o:Chthoniobacterales"="Chthoniobacterales",
                                        "Other_Taxa"="Other taxa"
)



#Let's set our color pallete 







All_site_bact_soil_root_C_order_color=c("p:Acidobacteria.o:Subgroup_4"="#1B9E77",
                                        "p:Acidobacteria.o:Subgroup_6"="#147759",
                                        "p:Actinobacteria.o:Kineosporiales"="#D95F02",
                                        "p:Actinobacteria.o:Micrococcales"="#A34702",
                                        "p:Actinobacteria.o:Micromonosporales"="#6D3001",
                                        "p:Actinobacteria.o:Pseudonocardiales"="#361801",
                                        "p:Actinobacteria.o:Solirubrobacterales"="#1B0C00",
                                        "p:Bacteroidetes.o:Sphingobacteriales"="#7570B3",
                                        "p:Firmicutes.o:Bacillales"="#E7298A",
                                        "p:Proteobacteria.o:Burkholderiales"="#66A61E",
                                        "p:Proteobacteria.o:Caulobacterales"="#59911A",
                                        "p:Proteobacteria.o:Myxococcales"="#4D7D16",
                                        "p:Proteobacteria.o:Pseudomonadales"="#406813",
                                        "p:Proteobacteria.o:Rhizobiales"="#33530F",
                                        "p:Proteobacteria.o:Sphingomonadales"="#263E0B",
                                        "p:Proteobacteria.o:Xanthomonadales"="#1A2A07",
                                        "p:Verrucomicrobia.o:Chthoniobacterales"="#E6AB02",
                                        "Other_Taxa"="grey"
)


length(All_site_bact_root_C_order_color)

site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")

(o_All_Sites__bact_root_color=ggplot(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M2,aes(x=factor(variable,levels = site_order),
                                                                                                  y=variable_std))+
    geom_bar(stat = "identity",aes( fill=factor(Real_P_O,levels = All_site_bact_root_C_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank(),legend.position = "none")+xlab(NULL)+ylab("Proportion")+
    scale_fill_manual(values = All_site_bact_soil_root_C_order_color,labels=All_site_bact_soil_root_C_order_names)+
    guides(fill=guide_legend(title="Order")))

#1300x650

(o_All_Sites__bact_root_color2=ggplot(GLBRC018_OTU_bact_root_core_facet.order.prop_otu2_sum_M2,aes(x=factor(variable,levels = site_order),
                                                                                                   y=variable_std*100))+
    geom_bar(stat = "identity",aes( fill=factor(Real_P_O,levels = All_site_bact_root_C_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_blank(),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = All_site_bact_soil_root_C_order_color,labels=All_site_bact_soil_root_C_order_names)+
    guides(fill=guide_legend(title="Order")))


(o_All_Sites__bact_soil_color2=ggplot(GLBRC018_OTU_bact_soil_core_facet.order.prop_otu2_sum_M2,
                                      aes(x=factor(variable,levels = site_order,
                                                   labels =site_labels),
                                          y=variable_std*100))+
    geom_bar(stat = "identity",aes( fill=factor(Real_P_O,levels = All_site_bact_soil_C_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = All_site_bact_soil_root_C_order_color,labels=All_site_bact_soil_root_C_order_names)+
    guides(fill=guide_legend(title="Order")))



plot_grid(o_All_Sites__bact_root_color2,o_All_Sites__bact_soil_color2,nrow = 2,labels = c("a)","b)"),
          align = "v")

#1200x1100

#####Stack Soil Fungi Core Taxa plot####

#Need to turn the subsetted community into a phyloseq obj

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5[[4]],fill=="core")$otu,
                                                                GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)
#I am getting a lot of Unknown phyla with sintax so I am going to use CONSTAX created dataset


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl=phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl),
                                                              sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl),
                                                              tax_table(as.matrix(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)))


ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl)
#60
ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl)/ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)
#0.02259036
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl))
# 522566
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil))
#0.4665768

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl))
#6933
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl))
#1887
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl)
#112

#Calculated the relative abundance of the core community
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl_map=(sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl))
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl)
fung_All_sites_soil_core_freq=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl_map%>%group_by(siteID)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)

#merge OTUs by site
GLBRC018_OTU_fung_soil_core_facet=merge_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_phyl, "siteID")
sample_names(GLBRC018_OTU_fung_soil_core_facet)     


#Order Level
#combine the reads at Order level
get_taxa_unique(GLBRC018_OTU_fung_soil_core_facet, taxonomic.rank="Order")
#[1] "Unknown"             "Glomerales"          "Hypocreales"         "Mortierellales"      "Agaricales"          "Cystofilobasidiales" "Helotiales"         
#[8] "Chaetosphaeriales"   "Dothideales"         "Filobasidiales"      "Glomerellales"       "Melanosporales"      "Chaetothyriales"     "Pleosporales"       
#[15] "Polyporales"         "Sordariales" 
(GLBRC018_OTU_fung_soil_core_facet.order<-tax_glom(GLBRC018_OTU_fung_soil_core_facet, taxrank="Order"))
#17



GLBRC018_OTU_fung_soil_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_fung_soil_core_facet.order))$Phylum,
                                                             "Order"=data.frame(tax_table(GLBRC018_OTU_fung_soil_core_facet.order))$Order,
                                                             "OTU"=taxa_names(GLBRC018_OTU_fung_soil_core_facet.order))

GLBRC018_OTU_fung_soil_core_facet.order_Names_raw$P_C=with(GLBRC018_OTU_fung_soil_core_facet.order_Names_raw,interaction(Phylum,Order))
#Let's take the top 10 taxa 

AllSite_soil_fung_Toporder = names(sort(taxa_sums(GLBRC018_OTU_fung_soil_core_facet.order), TRUE)[1:10])
AllSite_soil_fung_TOP10_order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(AllSite_soil_fung_Toporder, GLBRC018_OTU_fung_soil_core_facet.order)))$Phylum,
                                               "Order"=data.frame(tax_table(prune_taxa(AllSite_soil_fung_Toporder, GLBRC018_OTU_fung_soil_core_facet.order)))$Order,
                                               "OTU"=taxa_names(prune_taxa(AllSite_soil_fung_Toporder, GLBRC018_OTU_fung_soil_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_fung_soil_core_facet.order_Names=merge(GLBRC018_OTU_fung_soil_core_facet.order_Names_raw,AllSite_soil_fung_TOP10_order_Names, by="OTU", all.x = T)



#Transform the read counts to prop of total reads

GLBRC018_OTU_fung_soil_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_fung_soil_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_fung_soil_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_fung_soil_core_facet.order.prop)))
GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2=merge(GLBRC018_OTU_fung_soil_core_facet.order_Names[,c("P_C","OTU")],
                                                        GLBRC018_OTU_fung_soil_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2$OTU=NULL

GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_trunc=
  GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2[order(-rowSums(
    GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2)])),][1:10,]

GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_trunc1=GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))


GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M=pivot_longer(data.frame(GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_trunc1),
                                                                 !P_C, names_to = "variable")
sort(unique(GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M$P_C))

GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M2=merge(GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M,
                                                           fung_All_sites_soil_core_freq,by.x="variable",
                                                           by.y="siteID")
GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M2=GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M2%>%
  mutate(value_std=value*samp_sum_mean)


All_site_fung_soil_O_order_order=c("Ascomycota.Chaetosphaeriales","Ascomycota.Chaetothyriales",
                                   "Ascomycota.Dothideales","Ascomycota.Helotiales",
                                   "Ascomycota.Hypocreales","Ascomycota.Sordariales",
                                   "Ascomycota.Unknown","Basidiomycota.Agaricales",
                                   "Mortierellomycota.Mortierellales", "Other_Taxa","Unknown.Unknown")
length(All_site_fung_soil_O_order_order)

#Names



All_site_fung_soil_root_O_order_names=c("Ascomycota.Chaetosphaeriales"="Chaetosphaeriales",
                                        "Ascomycota.Chaetothyriales"="Chaetothyriales",
                                        "Ascomycota.Dothideales"="Dothideales",
                                        "Ascomycota.Helotiales"="Helotiales",
                                        "Ascomycota.Hypocreales"="Hypocreales",
                                        "Ascomycota.Pleosporales"="Pleosporales",
                                        "Ascomycota.Sordariales"="Sordariales",
                                        "Ascomycota.Unknown"="Unclassified Ascomycota",
                                        "Basidiomycota.Agaricales"="Agaricales",
                                        "Basidiomycota.Auriculariales"="Auriculariales",
                                        "Glomeromycota.Glomerales"="Glomerales",
                                        "Mortierellomycota.Mortierellales"="Mortierellales", 
                                        "Other_Taxa"="Other taxa",
                                        "Unknown.Unknown"="Unclassified")



#Color


All_site_fung_soil_root_O_order_color=c("Ascomycota.Chaetosphaeriales"="#D95F02",
                                        "Ascomycota.Chaetothyriales"="#BE5302",
                                        "Ascomycota.Dothideales"="#A34702",
                                        "Ascomycota.Helotiales"="#883B01",
                                        "Ascomycota.Hypocreales"="#6D3001",
                                        "Ascomycota.Pleosporales"="#512401",
                                        "Ascomycota.Sordariales"="#361801",
                                        "Ascomycota.Unknown"="#1B0C00",
                                        "Basidiomycota.Agaricales"="#7570B3",
                                        "Basidiomycota.Auriculariales"="#3A3859",
                                        "Glomeromycota.Glomerales"="#E7298A",
                                        "Mortierellomycota.Mortierellales"="#66A61E", 
                                        "Other_Taxa"="darkgrey",
                                        "Unknown.Unknown"="lightgrey")

site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")


(o_All_Sites__fung_soil_color=ggplot(GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M2,
                                     aes(x=factor(variable,levels = site_order),y=value_std))+
    geom_bar(stat = "identity",aes( fill=factor(P_C,levels = All_site_fung_soil_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = All_site_fung_soil_root_O_order_color,labels=All_site_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))

#1300x650








#####Stack Root Fungi Core Taxa plot####

#Need to turn the subsetted community into a phyloseq obj
subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],fill=="core")$otu
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],fill=="core")$otu,
                                                                GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl=phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl),
                                                              sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl),
                                                              tax_table(as.matrix(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)))


ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl)
#123
ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl)/ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)
#0.05336226
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl))
#658508
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root))
#0.5932505

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl))
#9319
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl))
#1423
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl)
#111


#Calculated the relative abundance of the core community
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl_map=(sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl))
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl)
fung_All_sites_root_core_freq=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl_map%>%group_by(siteID)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)

#merge OTUs by site
GLBRC018_OTU_fung_root_core_facet=merge_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_phyl, "siteID")
sample_names(GLBRC018_OTU_fung_root_core_facet)     

#Order Level
#combine the reads at Order level
get_taxa_unique(GLBRC018_OTU_fung_root_core_facet, taxonomic.rank="Order")
#"Unknown"               "Glomerales"            "Hypocreales"           "Mortierellales"        "Auriculariales"        "Archaeorhizomycetales"
#[7] "Helotiales"            "Sordariales"           "Sebacinales"           "Chaetothyriales"       "Xylariales"            "Pleosporales"         
#[13] "Paraglomerales"        "Dothideales"           "Agaricales"            "Capnodiales"           "Glomerellales"         "Chaetosphaeriales" 
(GLBRC018_OTU_fung_root_core_facet.order<-tax_glom(GLBRC018_OTU_fung_root_core_facet, taxrank="Order"))
#22



GLBRC018_OTU_fung_root_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_fung_root_core_facet.order))$Phylum,
                                                             "Order"=data.frame(tax_table(GLBRC018_OTU_fung_root_core_facet.order))$Order,
                                                             "OTU"=taxa_names(GLBRC018_OTU_fung_root_core_facet.order))

GLBRC018_OTU_fung_root_core_facet.order_Names_raw$P_C=with(GLBRC018_OTU_fung_root_core_facet.order_Names_raw,interaction(Phylum,Order))
#Let's take the top 10 taxa 

AllSite_root_fung_Toporder = names(sort(taxa_sums(GLBRC018_OTU_fung_root_core_facet.order), TRUE)[1:10])
AllSite_root_fung_TOP10_order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(AllSite_root_fung_Toporder, GLBRC018_OTU_fung_root_core_facet.order)))$Phylum,
                                               "Order"=data.frame(tax_table(prune_taxa(AllSite_root_fung_Toporder, GLBRC018_OTU_fung_root_core_facet.order)))$Order,
                                               "OTU"=taxa_names(prune_taxa(AllSite_root_fung_Toporder, GLBRC018_OTU_fung_root_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_fung_root_core_facet.order_Names=merge(GLBRC018_OTU_fung_root_core_facet.order_Names_raw,AllSite_root_fung_TOP10_order_Names, by="OTU", all.x = T)



#Transform the read counts to prop of total reads

GLBRC018_OTU_fung_root_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_fung_root_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_fung_root_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_fung_root_core_facet.order.prop)))
GLBRC018_OTU_fung_root_core_facet.order.prop_otu2=merge(GLBRC018_OTU_fung_root_core_facet.order_Names[,c("P_C","OTU")],
                                                        GLBRC018_OTU_fung_root_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

GLBRC018_OTU_fung_root_core_facet.order.prop_otu2$OTU=NULL

GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_trunc=
  GLBRC018_OTU_fung_root_core_facet.order.prop_otu2[order(-rowSums(
    GLBRC018_OTU_fung_root_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2)])),][1:10,]

GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_trunc1=GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))


GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M=pivot_longer(data.frame(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_trunc1),
                                                         !P_C,names_to = "variable")
sort(unique(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M$P_C))

GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M2=merge(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M,
                                                           fung_All_sites_root_core_freq,by.x="variable",
                                                           by.y="siteID")
GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M2=GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M2%>%
  mutate(value_std=value*samp_sum_mean)
sort(unique(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M$P_C))

All_site_fung_root_O_order_order=sort(unique(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M$P_C))

length(All_site_fung_root_O_order_order)

#Names



All_site_fung_soil_root_O_order_names=c("Ascomycota.Chaetosphaeriales"="Chaetosphaeriales",
                                        "Ascomycota.Chaetothyriales"="Chaetothyriales",
                                        "Ascomycota.Dothideales"="Dothideales",
                                        "Ascomycota.Helotiales"="Helotiales",
                                        "Ascomycota.Hypocreales"="Hypocreales",
                                        "Ascomycota.Pleosporales"="Pleosporales",
                                        "Ascomycota.Sordariales"="Sordariales",
                                        "Ascomycota.Unknown"="Unclassified Ascomycota",
                                        "Basidiomycota.Agaricales"="Agaricales",
                                        "Basidiomycota.Auriculariales"="Auriculariales",
                                        "Glomeromycota.Glomerales"="Glomerales",
                                        "Mortierellomycota.Mortierellales"="Mortierellales", 
                                        "Other_Taxa"="Other taxa",
                                        "Unknown.Unknown"="Unclassified")



#Color

All_site_fung_soil_root_O_order_color=c("Ascomycota.Chaetosphaeriales"="#D95F02",
                                        "Ascomycota.Chaetothyriales"="#BE5302",
                                        "Ascomycota.Dothideales"="#A34702",
                                        "Ascomycota.Helotiales"="#883B01",
                                        "Ascomycota.Hypocreales"="#6D3001",
                                        "Ascomycota.Pleosporales"="#512401",
                                        "Ascomycota.Sordariales"="#361801",
                                        "Ascomycota.Unknown"="#1B0C00",
                                        "Basidiomycota.Agaricales"="#7570B3",
                                        "Basidiomycota.Auriculariales"="#3A3859",
                                        "Glomeromycota.Glomerales"="#E7298A",
                                        "Mortierellomycota.Mortierellales"="#66A61E", 
                                        "Other_Taxa"="darkgrey",
                                        "Unknown.Unknown"="lightgrey")

site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")

(o_All_Sites__fung_root_color=ggplot(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M2,
                                     aes(x=factor(variable,levels = site_order),y=value_std))+
    geom_bar(stat = "identity",aes( fill=factor(P_C,levels = All_site_fung_root_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_fill_manual(values = All_site_fung_soil_root_O_order_color,labels=All_site_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))

#1300x650



(o_All_Sites__fung_root_color2=ggplot(GLBRC018_OTU_fung_root_core_facet.order.prop_otu2_M2,
                                      aes(x=factor(variable,levels = site_order),y=value_std*100))+
    geom_bar(stat = "identity",aes( fill=factor(P_C,levels = All_site_fung_root_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_blank(),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = All_site_fung_soil_root_O_order_color,labels=All_site_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))


(o_All_Sites__fung_soil_color2=ggplot(GLBRC018_OTU_fung_soil_core_facet.order.prop_otu2_M2,
                                      aes(x=factor(variable,levels = site_order,
                                                   labels = site_labels),y=value_std*100))+
    geom_bar(stat = "identity",aes( fill=factor(P_C,levels = All_site_fung_soil_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = All_site_fung_soil_root_O_order_color,labels=All_site_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))


plot_grid(o_All_Sites__fung_root_color2,o_All_Sites__fung_soil_color2,nrow = 2,labels = c("a)","b)"),label_size = 28,
          align = "v")




plot_grid(o_All_Sites__bact_root_color2,o_All_Sites__fung_root_color2,o_All_Sites__bact_soil_color2,o_All_Sites__fung_soil_color2,
          nrow = 2,labels = c("a)","b)"),label_size = 28,
          align = "v")
ggsave(plot_grid(o_All_Sites__bact_root_color2,o_All_Sites__fung_root_color2,o_All_Sites__bact_soil_color2,o_All_Sites__fung_soil_color2,
                 nrow = 2,labels = c("a)","b)","c)","d)"),label_size = 48,label_x = c(-0.02,-0.02,-0.02,-0.02),
                 align = "v"), filename = "All_sites_bacteria_fungi_root_soil_core_stack_order_p.png",path = here::here("Manuscript","Core_comm_figs"),width = 35,height = 15)


ggsave(plot_grid(o_All_Sites__bact_root_color2,o_All_Sites__fung_root_color2,o_All_Sites__bact_soil_color2,o_All_Sites__fung_soil_color2,
                 nrow = 2,labels = c("a)","b)","c)","d)"),label_size = 48,label_x = c(-0.02,-0.02,-0.02,-0.02),
                 align = "v"), filename = "All_sites_bacteria_fungi_root_soil_core_stack_order_p.svg",path = here::here("Manuscript","Core_comm_figs"),width = 35,height = 15)



#####metacoder Soil versus Root Bacteria  Core Taxa####


#Root Core taxa in Soil community community


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_core_5[[4]],fill=="core")$otu,
                                                           GLBRC018_OTU_bact_MMPRNT_All_sites_G5)

ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT)
#33
ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT)/ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
#0.001036465

sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT))
#471430
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5))
#0.2076784

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT))
#5207
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT))
#464
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT)
#227

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC=data.frame(tax_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT))
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC$Species=
  with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC,
       paste(Genus,row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC),sep = ""))


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC=as.matrix(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC|>mutate_all(str_remove,"^(.+):"))
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT,function(x)(x/10000))),
           tax_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC),
           sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT))
ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC_phyl)
#33
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC=parse_phyloseq(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_taxC_phyl)
print(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC)

site_order=c("LUX","LC","ESC", "HAN","RHN")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$siteID=
  factor(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$siteID,
         levels = site_order,ordered = T)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$siteID),]

set.seed(2021) # This makes the plot appear the same each time it is run 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC |> 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations




GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID,
                   groups = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Project,
                   out_names = "Mean_across_samples")


#Need the mean taxon abundance
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID)




GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                           log2_median_ratio))

summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  heat_tree(node_label = taxon_names,
            node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_color = log2_median_ratio_Fam, # A column from `obj$data$diff_table`
            node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median proportions",
            node_color_range = diverging_palette(), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            layout = "davidson-harel", # The primary layout algorithm
            title = "Root core: Root and Soil Communities",
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam<0)

#Known Taxa only 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                              log2_median_ratio))



subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table,log2_median_ratio_Fam_KN<0)




#Soil Core taxa in Root community


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_core_5[[4]],fill=="core")$otu,
                                                           GLBRC018_OTU_bact_MMPRNT_All_sites_G5)

ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL)
#133
ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL)/ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
#0.004177267

sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL))
#726632
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5))
#0.3201022

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL))
#4834
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL))
#2230
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL)
#227

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC=data.frame(tax_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL))
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC$Species=
  with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC,
       paste(Genus,row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC),sep = ""))


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC=as.matrix(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC|>mutate_all(str_remove,"^(.+):"))
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL,function(x)(x/10000))),
           tax_table(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC),
           sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL))
ntaxa(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC_phyl)
#133
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC=parse_phyloseq(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_taxC_phyl)
print(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC)

site_order=c("LUX","LC","ESC", "HAN","RHN")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$siteID=
  factor(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$siteID,
         levels = site_order,ordered = T)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$siteID),]

set.seed(2021) # This makes the plot appear the same each time it is run 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC |> 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations




GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID,
                   groups = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Project,
                   out_names = "Mean_across_samples")



#Need the mean taxon abundance
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID)




GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                           log2_median_ratio))

summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  heat_tree(node_label = taxon_names,
            node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_color = log2_median_ratio_Fam, # A column from `obj$data$diff_table`
            node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median proportions",
            node_color_range = diverging_palette(), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            layout = "davidson-harel", # The primary layout algorithm
            title = "Soil core: Root and Soil Communities",
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations


#Known Taxa only 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                              log2_median_ratio))






(Bact_core_ROOT_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              title = "Root core: Root and Soil Communities",
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

(Bact_core_SOIL_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              title = "Soil core: Root and Soil Communities",
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

plot_grid(Bact_core_ROOT_metaHeat_P2,Bact_core_SOIL_metaHeat_P2,
          align = "h")





#####metacoder Soil versus Root Fungi Core Taxa####

GLBRC018_OTU_fung_MMPRNT_All_sites_G5=phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5),
                                               sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5),
                                               tax_table(as.matrix(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)))



#Root Core taxa in Soil community community

subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],fill=="core")$otu
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_core_5[[4]],fill=="core")$otu,
                                                           GLBRC018_OTU_fung_MMPRNT_All_sites_G5)

ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT)
#123
ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT)/ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)
#0.03980583

sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT))
#1129417
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5))
#0.506465

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT))
#9319
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT))
#1423
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT)
#223

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC=data.frame(tax_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT))
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC$Species=
  with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC,
       paste(Species,row.names(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC),sep = ""))



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT,function(x)(x/10000))),
           tax_table(as.matrix(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC)),
           sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT))
ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC_phyl)
#123
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC=parse_phyloseq(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_taxC_phyl)
print(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC)



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil=
  factor(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil,
         levels = c("Root","Soil"),ordered = T)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil),]

site_order=c("LUX","LC","ESC", "HAN","RHN")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$siteID=
  factor(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$siteID,
         levels = site_order,ordered = T)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$siteID),]

set.seed(2021) # This makes the plot appear the same each time it is run 

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq,
                   groups = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$project,
                   out_names = "Mean_across_samples")


#Need the mean taxon abundance
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq)




GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                           log2_median_ratio))

summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  heat_tree(node_label = taxon_names,
            node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_color = log2_median_ratio_Fam, # A column from `obj$data$diff_table`
            node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median proportions",
            node_color_range = diverging_palette(), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            layout = "davidson-harel", # The primary layout algorithm
            title = "Root core: Root and Soil Communities",
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam<0)

#Known Taxa only 
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                              log2_median_ratio))




#Soil Core taxa in Root community



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_core_5[[4]],fill=="core")$otu,
                                                           GLBRC018_OTU_fung_MMPRNT_All_sites_G5)

ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL)
#60
ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL)/ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)
#0.02103787

sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL))
#649536
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5))
#0.2912717

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL))
#6933

min(sample_sums(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL))
#252
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL)
#223

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC=data.frame(tax_table(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL))
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC$Species=
  with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC,
       paste(Species,row.names(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC),sep = ""))



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL,function(x)(x/10000))),
           tax_table(as.matrix(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC)),
           sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL))
ntaxa(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC_phyl)
#60
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC=parse_phyloseq(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_taxC_phyl)
print(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil=
  factor(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil,
         levels = c("Root","Soil"),ordered = T)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil),]

site_order=c("LUX","LC","ESC", "HAN","RHN")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$siteID=
  factor(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$siteID,
         levels = site_order,ordered = T)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$siteID),]


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq,
                   groups = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$project,
                   out_names = "Mean_across_samples")



#Need the mean taxon abundance
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq)




GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                           log2_median_ratio))

summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  heat_tree(node_label = taxon_names,
            node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_color = log2_median_ratio_Fam, # A column from `obj$data$diff_table`
            node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median proportions",
            node_color_range = diverging_palette(), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            layout = "davidson-harel", # The primary layout algorithm
            title = "Soil core: Root and Soil Communities",
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam>0)

#Known Taxa only 
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                              log2_median_ratio))



(Fung_core_ROOT_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              title = "Root core: Root and Soil Communities",
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

(Fung_core_SOIL_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              title = "Soil core: Root and Soil Communities",
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

plot_grid(Fung_core_ROOT_metaHeat_P2,Fung_core_SOIL_metaHeat_P2,
          align = "h")




#####Figures metacoder Soil versus Root Bacteria Fungi Core Taxa#####
set.seed(2021)
(Bact_core_ROOT_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
   filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
   filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
   heat_tree(node_label = taxon_names,
             node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
             node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
             node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
             node_size_axis_label = "OTU count",
             node_color_axis_label = "Log 2 ratio of median proportions",
             node_color_range = diverging_palette(), # The built-in palette for diverging data
             node_color_trans = "linear", # The default is scaled by circle area
             layout = "davidson-harel", # The primary layout algorithm
             initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

(Bact_core_SOIL_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations


(Fung_core_ROOT_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_ROOT_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

(Fung_core_SOIL_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_CORE_SOIL_phyl_MC|>
    filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
    filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
    heat_tree(node_label = taxon_names,
              node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
              node_color = log2_median_ratio_Fam_KN, # A column from `obj$data$diff_table`
              node_color_interval = c(-4, 4), # The range of `log2_median_ratio` to display
              node_size_axis_label = "OTU count",
              node_color_axis_label = "Log 2 ratio of median proportions",
              node_color_range = diverging_palette(), # The built-in palette for diverging data
              node_color_trans = "linear", # The default is scaled by circle area
              layout = "davidson-harel", # The primary layout algorithm
              initial_layout = "reingold-tilford")) # The layout algorithm that initializes node locations

Bact_title=ggplot()+ggtitle("Bacterial core community")+theme(plot.title = element_text(size = 36,hjust = 0.5))
Fung_title=ggplot()+ggtitle("Fungi core community")+theme(plot.title = element_text(size = 36,hjust = 0.5))
Root_90_title=ggplot()+ggtitle("Root")+theme(plot.title = element_text(angle = 90,size = 36))
Soil_90_title=ggplot()+ggtitle("Soil")+theme(plot.title = element_text(angle = 90,size = 36))
plot_grid(ggplot(),get_title(Bact_title),get_title(Fung_title),
          get_title(Root_90_title),Bact_core_ROOT_metaHeat_P2,Fung_core_ROOT_metaHeat_P2,
          get_title(Soil_90_title),Bact_core_SOIL_metaHeat_P2,Fung_core_SOIL_metaHeat_P2,
          rel_widths = c(0.1,1,1),rel_heights = c(0.1,1,1),ncol = 3)



ggsave(plot_grid(ggplot()+theme_void(),get_title(Bact_title),get_title(Fung_title),
                 get_title(Root_90_title),Bact_core_ROOT_metaHeat_P2,Fung_core_ROOT_metaHeat_P2,
                 get_title(Soil_90_title),Bact_core_SOIL_metaHeat_P2,Fung_core_SOIL_metaHeat_P2,
                 rel_widths = c(0.1,1,1),rel_heights = c(0.2,1,1),ncol = 3), 
       filename = "MLE_Heat_tree_Bact_Fung_Root_Soil_core_comm.png",path = here::here("Manuscript","Core_comm_figs"),width = 15,height = 10)





######MLE CORE Funguild Fig#### 



#Funguild Database



head(GLBRC_CON_FunGuild_raw_F)




#Extract OTU Table
GLBRC018_OTU_fung_soil_core_facet_otus=data.frame(t(otu_table(GLBRC018_OTU_fung_soil_core_facet)))
GLBRC018_OTU_fung_soil_core_facet_otus[1:10,1:4]



#Extract OTU Table
GLBRC018_OTU_fung_root_core_facet_otus=data.frame(t(otu_table(GLBRC018_OTU_fung_root_core_facet)))
GLBRC018_OTU_fung_root_core_facet_otus[1:10,1:4]

#Combine OTU and FungalTraits

GLBRC018_OTU_fung_soil_core_facet_otus_FunG=merge(GLBRC018_OTU_fung_soil_core_facet_otus,
                                                  GLBRC_CON_FunGuild_raw_F,
                                                  by.x = "row.names", by.y = "OTU",
                                                  all.x = T)

GLBRC018_OTU_fung_root_core_facet_otus_FunG=merge(GLBRC018_OTU_fung_root_core_facet_otus,
                                                  GLBRC_CON_FunGuild_raw_F,
                                                  by.x = "row.names", by.y = "OTU",
                                                  all.x = T)

summary(GLBRC018_OTU_fung_root_core_facet_otus_FunG)
unique(GLBRC018_OTU_fung_root_core_facet_otus_FunG$simp_guild)
GLBRC018_OTU_fung_soil_core_facet_otus_FunG[is.na(GLBRC018_OTU_fung_soil_core_facet_otus_FunG)]="Unknown"
GLBRC018_OTU_fung_root_core_facet_otus_FunG[is.na(GLBRC018_OTU_fung_root_core_facet_otus_FunG)]="Unknown"
#Simplified Guild 

GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum=GLBRC018_OTU_fung_soil_core_facet_otus_FunG%>%
  group_by(simp_guild)%>%summarise(across(ESC:RHN,~sum(.)))

GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M=GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum%>%
  pivot_longer(!simp_guild, names_to = "siteID", values_to = "count")
GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M=GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M%>%
  group_by(siteID)%>%mutate(prop=count/sum(count))

GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum=GLBRC018_OTU_fung_root_core_facet_otus_FunG%>%
  group_by(simp_guild)%>%summarise(across(ESC:RHN,~sum(.)))

GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M=GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum%>%
  pivot_longer(!simp_guild, names_to = "siteID", values_to = "count")
GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M=GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M%>%
  group_by(siteID)%>%mutate(prop=count/sum(count))

GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M2=merge(
  GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M,fung_All_sites_soil_core_freq, by="siteID")

dim(GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M2)
#45  5

GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M2=merge(
  GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M,fung_All_sites_root_core_freq, by="siteID")

dim(GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M2)
#55   5



GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M2=GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M2%>%mutate(value_std=prop*samp_sum_mean)
GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M2=GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M2%>%mutate(value_std=prop*samp_sum_mean)
sort(unique(GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M2$simp_guild))



All_site_fung_soil_FGuild_order=c("Arbuscular Mycorrhizal",
                                  "Animal Pathogen",
                                  "Multiple Saprotroph","Undefined Saprotroph",
                                  "Symbiotroph-Pathogen",
                                  "Symbiotroph-Saprotroph",
                                  "Pathogen-Saprotroph",
                                  "Symbiotroph-Pathogen-Saprotroph",
                                  "Unknown")


All_site_fung_soil_FGuild_color=c("Arbuscular Mycorrhizal"="#1B9E77",
                                  "Animal Pathogen"="#D95F02",
                                  "Multiple Saprotroph"="#7570B3","Undefined Saprotroph"="#3A3859",
                                  "Symbiotroph-Pathogen"="#E7298A",
                                  "Symbiotroph-Saprotroph"="#66A61E",
                                  "Pathogen-Saprotroph"="#E6AB02",
                                  "Symbiotroph-Pathogen-Saprotroph"="#A6761D",
                                  "Unknown"="lightgrey")



site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")


(FG_All_Sites_FG_fung_soil_color2=ggplot(GLBRC018_OTU_fung_soil_core_facet_otus_FunG_sum_M2,
                                         aes(x=factor(siteID,levels = site_order),y=(value_std*100)))+
    geom_bar(stat = "identity",aes( fill=factor(simp_guild,levels = All_site_fung_soil_FGuild_order)),color="black")+
    theme_bw()+scale_y_continuous(limits = c(0,70),name = "Soil\nPercentage of reads")+
    scale_x_discrete(labels=site_labels)+
    theme(axis.text.y=element_text(size=28,color = "black"),axis.text.x=element_text(size=30,color = "black"),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+
    scale_fill_manual(values = All_site_fung_soil_FGuild_color)+
    guides(fill=guide_legend(title="FunGuild")))



All_site_fung_root_FGuild_order=c("Arbuscular Mycorrhizal","Ectomycorrhizal","Endophyte", 
                                  "Plant Pathogen",
                                  "Multiple Saprotroph","Undefined Saprotroph",
                                  "Symbiotroph-Pathogen",
                                  "Symbiotroph-Saprotroph",
                                  "Pathogen-Saprotroph",
                                  "Symbiotroph-Pathogen-Saprotroph",
                                  "Unknown")


All_site_fung_root_FGuild_color=c("Arbuscular Mycorrhizal"="#1B9E77","Ectomycorrhizal"="#11634A","Endophyte"="#07281E", 
                                  "Plant Pathogen"="#D95F02",
                                  "Multiple Saprotroph"="#7570B3","Undefined Saprotroph"="#3A3859",
                                  "Symbiotroph-Pathogen"="#E7298A",
                                  "Symbiotroph-Saprotroph"="#66A61E",
                                  "Pathogen-Saprotroph"="#E6AB02",
                                  "Symbiotroph-Pathogen-Saprotroph"="#A6761D",
                                  "Unknown"="lightgrey")






(FG_All_Sites_FG_fung_root_color2=ggplot(GLBRC018_OTU_fung_root_core_facet_otus_FunG_sum_M2,
                                         aes(x=factor(siteID,levels = site_order),y=(value_std*100)))+
    geom_bar(stat = "identity",aes( fill=factor(simp_guild,levels = All_site_fung_root_FGuild_order)),color="black")+
    theme_bw()+scale_y_continuous(limits = c(0,70),name = "Root\nPercentage of reads")+
    scale_x_discrete(labels=NULL,name=NULL)+
    theme(axis.text.y=element_text(size=28,color = "black"),axis.text.x=element_text(size=30,color = "black"),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+
    scale_fill_manual(values = All_site_fung_root_FGuild_color)+
    guides(fill=guide_legend(title="FunGuild")))



plot_grid(FG_All_Sites_FG_fung_root_color2,FG_All_Sites_FG_fung_soil_color2,
          label_size = 48,labels = c("a)","b)"),
          nrow = 2,
          axis = "lr",align = "v",
          label_x = c(-0.015,-0.015))


ggsave(plot_grid(FG_All_Sites_FG_fung_root_color2,FG_All_Sites_FG_fung_soil_color2,
                 label_size = 48,labels = c("a)","b)"),
                 nrow = 2,
                 axis = "lr",align = "v",
                 label_x = c(-0.015,-0.015)),
       filename = "All_sites_FunGuild_fungi_root_soil_core_stack_order_sep_color.png",
       path = here::here("Manuscript","Core_comm_figs"),width = 19,height = 15)



