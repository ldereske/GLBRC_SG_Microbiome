library(here)
#here::i_am("R_code/OTU_bacterial_fungal_community_analyses_202110011.R")
library(phyloseq)
library(plyr); library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(vegan)
library(metacoder)
library(stringr)


"%w/o%" <- function(x,y)!('%in%'(x,y))
Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "English")

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


#Roots only 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,Root_soil=="Root")
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)
#142
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root) > 0, GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root=subset_samples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,Root_soil=="Root")
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)
#141
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root) > 0, GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)


#Soil only 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,Root_soil=="Soil")
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)
#353
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil) > 0, GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,Root_soil=="Soil")
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil)
#349
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil) > 0, GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil)

#####Lux Temporal Core Taxa####

#####Fitting Core Taxa Bacteria Soil Lux Temporal ####

##Prioritizing core microbiome based on the DATE

#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5=ExtractCoreFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil, "collectionDate", 5, Group=NULL, Level=NULL)

#Plotting the effects of BC thresholds on OTU inclusion
PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_2,1000,5)

#Fit the neutral model to the core taxa
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5_Neut=FitNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5)

PlotNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5_Neut)

#####Bacteria Root Lux Temporal Core Taxa####

##Prioritizing core microbiome based on the Collection Date


#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5=ExtractCoreFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root, "collectionDate", 5, Group=NULL, Level=NULL)

#Plotting the effects of BC thresholds on OTU inclusion
PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,1000,5)

#Fit the neutral model to the core taxa
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5_Neut=FitNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5)

PlotNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5_Neut)




plot_grid(PlotNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5_Neut),
          PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,1000,5),
          PlotNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5_Neut),
          PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5,1000,5),
          ncol = 2,
          labels = c("a)","b)","c)","d)"))

ggsave(plot_grid(PlotNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5_Neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,1000,5),
                 PlotNeutral(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5_Neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5,1000,5),
                 ncol = 2,
                 labels = c("a)","b)","c)","d)")), 
       filename = "Lux_Arbor_bacteria_core_fit.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)

plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,20,5),
          PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5,20,5),
          ncol = 2,
          labels = c("a)","b)","c)","d)"))


ggsave(plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5,20,5),
                 PlotBCThreshold_Rich(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5,20,5),
                 ncol = 2,
                 labels = c("a)","b)","c)","d)")), filename = "Lux_Arbor_bacteria_core_diagnostic.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)



#####Fungi Lux Temporal Core Taxa####


#####Fungi Soil Lux Temporal Core Taxa####


##Prioritizing core microbiome based on the Collectiondate

#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5=ExtractCoreFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil, "collectionDate", 5, Group=NULL, Level=NULL)

#Plotting the effects of BC thresholds on OTU inclusion
PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,1000,5)

#Fit the neutral model to the core taxa
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5_Neut=FitNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5)
 
PlotNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5_Neut)




####Fungi Root Lux Temporal Core Taxa####

##Prioritizing core microbiome based on the Collectiondate

#Extracting the OTUs and formatting the dataset
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5=ExtractCoreFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root, "collectionDate", 5, Group=NULL, Level=NULL)

#Plotting the effects of BC thresholds on OTU inclusion
PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_2,1000,5)

#Fit the neutral model to the core taxa
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5_neut=FitNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5)
 
PlotNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5_neut)


plot_grid(PlotNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5_neut),
          PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5,1000,5),
          PlotNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5_Neut),
          PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,1000,5),
          ncol = 2,
          labels = c("a)","b)","c)","d)"))


ggsave(plot_grid(PlotNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5_neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5,1000,5),
                 PlotNeutral(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5_Neut),
                 PlotBCincreaseFlex(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,1000,5),
                 ncol = 2,
                 labels = c("a)","b)","c)","d)")), 
       filename = "Lux_Arbor_fungi_core_fit.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)

plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5,20,5),
          PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,20,5),
          PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,20,5),
          ncol = 2,
          labels = c("a)","b)","c)","d)"))


ggsave(plot_grid(PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5,20,5),
                 PlotBCThreshold_Rich(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,20,5),
                 PlotBCThreshold_Abun(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5,20,5),
                 ncol = 2,
                 labels = c("a)","b)","c)","d)")), filename = "Lux_Arbor_fungi_core_diagnostic.png",path = here::here("Manuscript","Core_comm_figs"),width = 20,height = 15)









#####Bacteria Lux Temporal Core Taxon Barplots####

#####Stack Soil Bacteria Core Taxa plot####


#Need to turn the subsetted community into a phyloseq obj

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5[[4]],fill=="core")$otu,
                                                          GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)

ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl)
#74
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl)/ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)
#0.00326178

sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl))
#942843
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil))
#0.2670943

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl))
#3541
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl))
#1893
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl)
#353

#Calculated the relative abundance of the core community
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl_map=(sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl))
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl)

bact_LUX_soil_CD_core_freq=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl_map%>%group_by(collectionDate)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)
bact_LUX_soil_CD_core_freq$raw_date=str_replace_all(bact_LUX_soil_CD_core_freq$collectionDate,"/",".")
#merge OTUs by site
GLBRC018_OTU_LUX_bact_soil_core_facet=merge_samples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_phyl, "collectionDate")
sample_names(GLBRC018_OTU_LUX_bact_soil_core_facet)     

#Order Level
#combine the reads at Class level
get_taxa_unique(GLBRC018_OTU_LUX_bact_soil_core_facet, taxonomic.rank="Order")
#[1] "o:Corynebacteriales"   "o:Sphingomonadales"    "o:Rhizobiales"         "o:Burkholderiales"     "o:Subgroup_7"          "o:Micrococcales"      
#[7] "o:Myxococcales"        "o:Bacillales"          "Unknown"               "o:Chthoniobacterales"  "o:Gemmatimonadales"    "o:Subgroup_6"         
#[13] "o:Acidobacteriales"    "o:Subgroup_3"          "o:Xanthomonadales"     "o:Nitrosomonadales"    "o:Rhodospirillales"    "o:Frankiales"         
#[19] "o:Gaiellales"          "o:Solirubrobacterales" "o:Subgroup_2"          "o:TRA3-20"             "o:SC-I-84"    
(GLBRC018_OTU_LUX_bact_soil_core_facet.order<-tax_glom(GLBRC018_OTU_LUX_bact_soil_core_facet, taxrank="Order"))
#27 



GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_LUX_bact_soil_core_facet.order))$Phylum,
                                                                 "Class"=data.frame(tax_table(GLBRC018_OTU_LUX_bact_soil_core_facet.order))$Class,
                                                                 "Order"=data.frame(tax_table(GLBRC018_OTU_LUX_bact_soil_core_facet.order))$Order,
                                                                 "OTU"=taxa_names(GLBRC018_OTU_LUX_bact_soil_core_facet.order))

GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names_raw$P_C=with(GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names_raw,interaction(Phylum,Class,Order))
#Let's take the top 10 taxa 

LUX_soil_bact_TopCLASS = names(sort(taxa_sums(GLBRC018_OTU_LUX_bact_soil_core_facet.order), TRUE)[1:10])
LUX_soil_bact_TOP10_Order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(LUX_soil_bact_TopCLASS, 
                                                                                    GLBRC018_OTU_LUX_bact_soil_core_facet.order)))$Phylum,
                                           "Class"=data.frame(tax_table(prune_taxa(LUX_soil_bact_TopCLASS, 
                                                                                   GLBRC018_OTU_LUX_bact_soil_core_facet.order)))$Class,
                                           "Order"=data.frame(tax_table(prune_taxa(LUX_soil_bact_TopCLASS, GLBRC018_OTU_LUX_bact_soil_core_facet.order)))$Order,
                                           "OTU"=taxa_names(prune_taxa(LUX_soil_bact_TopCLASS, GLBRC018_OTU_LUX_bact_soil_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names=merge(GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names_raw,LUX_soil_bact_TOP10_Order_Names, by="OTU", all.x = T)

GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names$Real_P_C=ifelse(is.na(GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names$Order.y),
                                                                  paste(GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names$Phylum.x,
                                                                        "Other_spp",sep="."),                             
                                                                  as.character(GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names$P_C)
)

#Transform the read counts to prop of total reads

GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_LUX_bact_soil_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop)))
GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2=merge(GLBRC018_OTU_LUX_bact_soil_core_facet.order_Names[,c("Real_P_C","OTU")],
                                                            GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

row.names(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2)=GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2$OTU
GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2$OTU=NULL
GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum=GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))
colSums(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2)])
GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_trunc=
  GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2[order(-rowSums(
    GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2)])),][1:10,]
summary(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_trunc)
GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_trunc1=GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))

GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M=pivot_longer(data.frame(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_trunc1),
                                                                         !Real_P_C,names_to = "variable")

GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M$raw_date=str_replace(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M$variable,"X","")

GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2=merge(
  GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M,
  bact_LUX_soil_CD_core_freq, by.x="raw_date",by.y="raw_date")
summary(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2)
GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2=
  GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2%>%mutate(value_std=value*samp_sum_mean)

sort(unique(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2$Real_P_C))

Lux_Time_bact_soil_O_order_order=c("p:Acidobacteria.c:Acidobacteria.o:Acidobacteriales",           
                                   "p:Acidobacteria.c:Acidobacteria.o:Subgroup_3",
                                   "p:Acidobacteria.c:Acidobacteria.o:Subgroup_6",
                                   "p:Actinobacteria.c:Actinobacteria.o:Frankiales", 
                                   "p:Firmicutes.c:Bacilli.o:Bacillales",
                                   "p:Proteobacteria.c:Betaproteobacteria.o:Burkholderiales",
                                   "p:Proteobacteria.c:Alphaproteobacteria.o:Rhizobiales",
                                   "p:Proteobacteria.c:Betaproteobacteria.Unknown",
                                   "p:Proteobacteria.Unknown.Unknown",
                                   "p:Verrucomicrobia.c:Spartobacteria.o:Chthoniobacterales",
                                   "Other_Taxa" )

length(Lux_Time_bact_soil_O_order_order)




Lux_Time_bact_soil_root_O_order_names=c("p:Acidobacteria.c:Acidobacteria.o:Acidobacteriales"="Acidobacteriales",           
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_3"="Subgroup_3",
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_6"="Subgroup_6",
                                        "p:Actinobacteria.c:Actinobacteria.o:Frankiales"="Frankiales",
                                        "p:Actinobacteria.c:Actinobacteria.o:Kineosporiales"="Kineosporiales",
                                        "p:Actinobacteria.c:Actinobacteria.o:Micromonosporales"="Micromonosporales",
                                        "p:Actinobacteria.c:Actinobacteria.o:Pseudonocardiales"="Pseudonocardiales",
                                        "p:Firmicutes.c:Bacilli.o:Bacillales"="Bacillales",
                                        "p:Proteobacteria.c:Betaproteobacteria.o:Burkholderiales"="Burkholderiales",
                                        "p:Proteobacteria.c:Deltaproteobacteria.o:Myxococcales"="Myxococcales",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Pseudomonadales"="Pseudomonadales",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Rhizobiales"="Rhizobiales",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Sphingomonadales"="Sphingomonadales",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Xanthomonadales"="Xanthomonadales",
                                        "p:Proteobacteria.c:Betaproteobacteria.Unknown"="Unclassified Betaproteobacteria",
                                        "p:Proteobacteria.Unknown.Unknown"="Unclassified Proteobacteria",
                                        "p:Verrucomicrobia.c:Spartobacteria.o:Chthoniobacterales"="Chthoniobacterales",
                                        "Other_Taxa" ="Other taxa")


length(Lux_Time_bact_soil_root_O_order_names)


Lux_Time_bact_soil_root_O_order_color=c("p:Acidobacteria.c:Acidobacteria.o:Acidobacteriales"="#1B9E77",           
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_3"="#11634A",
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_6"="#07281E",
                                        "p:Actinobacteria.c:Actinobacteria.o:Frankiales"="#D95F02", 
                                        "p:Actinobacteria.c:Actinobacteria.o:Kineosporiales"="#A34702",
                                        "p:Actinobacteria.c:Actinobacteria.o:Micromonosporales"="#6D3001",
                                        "p:Actinobacteria.c:Actinobacteria.o:Pseudonocardiales"="#361801",
                                        "p:Firmicutes.c:Bacilli.o:Bacillales"="#E7298A", 
                                        "p:Proteobacteria.c:Betaproteobacteria.o:Burkholderiales"="#66A61E",
                                        "p:Proteobacteria.c:Deltaproteobacteria.o:Myxococcales"="#59911A",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Pseudomonadales"="#4D7D16",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Rhizobiales"="#406813",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Sphingomonadales"="#33530F",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Xanthomonadales"="#263E0B",
                                        "p:Proteobacteria.c:Betaproteobacteria.Unknown"="#1A2A07",
                                        "p:Proteobacteria.Unknown.Unknown"="#0D1504",
                                        "p:Verrucomicrobia.c:Spartobacteria.o:Chthoniobacterales"="#735501",
                                        "Other_Taxa" ="grey")
length(Lux_Time_bact_soil_root_O_order_color)

summary(as.Date(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2$raw_date,format="%m.%d.%Y"))


(o_LUX__bact_soil_color=ggplot(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2,
                               aes(x=(as.Date(raw_date,format="%m.%d.%Y")),y=value_std))+
    geom_area(aes( fill=factor(Real_P_C,levels = Lux_Time_bact_soil_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_fill_manual(values = Lux_Time_bact_soil_root_O_order_color,labels=Lux_Time_bact_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))

#1300x650









#####Stack Root Bacteria Core Taxa plot####

#Need to turn the subsetted community into a phyloseq obj

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5[[4]],fill=="core")$otu,
                                                          GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)

ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl)
#66
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl)/ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)
#0.004982636
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl))
#647679
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root))
#0.456112

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl))
#5834
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl))
#2437
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl)
#142



#Calculated the relative abundance of the core community
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl_map=(sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl))
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl)
bact_LUX_root_core_freq=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl_map%>%group_by(collectionDate)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)
bact_LUX_root_core_freq$raw_date=str_replace_all(bact_LUX_root_core_freq$collectionDate,"/",".")
#merge OTUs by site
GLBRC018_OTU_bact_LUX_root_core_facet=merge_samples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_phyl, "collectionDate")
sample_names(GLBRC018_OTU_bact_LUX_root_core_facet)     


#order Level
#combine the reads at Order level
get_taxa_unique(GLBRC018_OTU_bact_LUX_root_core_facet, taxonomic.rank="Order")
#[1] "o:Rhizobiales"         "o:Sphingomonadales"    "o:Caulobacterales"     "o:Corynebacteriales"   "o:Pseudonocardiales"   "o:Streptomycetales"   
#[7] "o:Pseudomonadales"     "o:Micromonosporales"   "o:Burkholderiales"     "o:Acidimicrobiales"    "Unknown"               "o:Bacillales"         
#[13] "o:Micrococcales"       "o:Kineosporiales"      "o:Xanthomonadales"     "o:Myxococcales"        "o:Chthoniobacterales"  "o:Gemmatimonadales"   
#[19] "o:Streptosporangiales" "o:Solirubrobacterales" "o:Rhodospirillales"    "o:Frankiales"          "o:Subgroup_10"         "o:Acidobacteriales"   
#[25] "o:Sphingobacteriales"  "o:Subgroup_3"          "o:Subgroup_6"          "o:Ktedonobacterales"   "o:Gaiellales" 
(GLBRC018_OTU_bact_LUX_root_core_facet.order<-tax_glom(GLBRC018_OTU_bact_LUX_root_core_facet, taxrank="Order"))
#31



GLBRC018_OTU_bact_LUX_root_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_bact_LUX_root_core_facet.order))$Phylum,
                                                                 "Class"=data.frame(tax_table(GLBRC018_OTU_bact_LUX_root_core_facet.order))$Class,
                                                                 "Order"=data.frame(tax_table(GLBRC018_OTU_bact_LUX_root_core_facet.order))$Order,
                                                                 "OTU"=taxa_names(GLBRC018_OTU_bact_LUX_root_core_facet.order))

GLBRC018_OTU_bact_LUX_root_core_facet.order_Names_raw$P_C=with(GLBRC018_OTU_bact_LUX_root_core_facet.order_Names_raw,interaction(Phylum,Class,Order))
#Let's take the top 10 taxa 

AllSite_root_bact_LUX_TopCLASS = names(sort(taxa_sums(GLBRC018_OTU_bact_LUX_root_core_facet.order), TRUE)[1:10])
AllSite_root_bact_LUX_TOP10_Order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(AllSite_root_bact_LUX_TopCLASS,
                                                                                            GLBRC018_OTU_bact_LUX_root_core_facet.order)))$Phylum,
                                                   "Class"=data.frame(tax_table(prune_taxa(AllSite_root_bact_LUX_TopCLASS,
                                                                                           GLBRC018_OTU_bact_LUX_root_core_facet.order)))$Class,
                                                   "Order"=data.frame(tax_table(prune_taxa(AllSite_root_bact_LUX_TopCLASS,
                                                                                           GLBRC018_OTU_bact_LUX_root_core_facet.order)))$Order,
                                                   "OTU"=taxa_names(prune_taxa(AllSite_root_bact_LUX_TopCLASS, GLBRC018_OTU_bact_LUX_root_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_bact_LUX_root_core_facet.order_Names=merge(GLBRC018_OTU_bact_LUX_root_core_facet.order_Names_raw,AllSite_root_bact_LUX_TOP10_Order_Names, by="OTU", all.x = T)

GLBRC018_OTU_bact_LUX_root_core_facet.order_Names$Real_P_C=ifelse(is.na(GLBRC018_OTU_bact_LUX_root_core_facet.order_Names$Order.y),
                                                                  paste(GLBRC018_OTU_bact_LUX_root_core_facet.order_Names$Phylum.x,
                                                                        "Other_spp",sep="."),                               
                                                                  as.character(GLBRC018_OTU_bact_LUX_root_core_facet.order_Names$P_C)
)

#Transform the read counts to prop of total reads

GLBRC018_OTU_bact_LUX_root_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_bact_LUX_root_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop)))
GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2=merge(GLBRC018_OTU_bact_LUX_root_core_facet.order_Names[,c("Real_P_C","OTU")],
                                                            GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2$OTU=NULL
GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum=GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))
colSums(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2)])
GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_trunc=
  GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2[order(-rowSums(
    GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2[,2:ncol(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2)])),][1:10,]

GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_trunc1=GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))

GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M=pivot_longer(data.frame(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_trunc1),
                                                                         !Real_P_C, names_to = "variable")
GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M$raw_date=str_replace(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M$variable,"X","")
sort(unique(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M$Real_P_C))

GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M2=merge(
  GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M,bact_LUX_root_core_freq,
  by="raw_date")
summary(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M2)
GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M2=GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M2%>%mutate(value_std=value*samp_sum_mean)


sort(unique(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M2$Real_P_C))

LUX_bact_root_O_order_order=c("p:Actinobacteria.c:Actinobacteria.o:Kineosporiales",
                              "p:Actinobacteria.c:Actinobacteria.o:Micromonosporales",
                              "p:Actinobacteria.c:Actinobacteria.o:Pseudonocardiales",
                              "p:Proteobacteria.c:Betaproteobacteria.o:Burkholderiales",
                              "p:Proteobacteria.c:Deltaproteobacteria.o:Myxococcales",
                              "p:Proteobacteria.c:Gammaproteobacteria.o:Pseudomonadales",
                              "p:Proteobacteria.c:Alphaproteobacteria.o:Rhizobiales",
                              "p:Proteobacteria.c:Alphaproteobacteria.o:Sphingomonadales",
                              "p:Proteobacteria.c:Gammaproteobacteria.o:Xanthomonadales",
                              "p:Proteobacteria.Unknown.Unknown",
                              "Other_Taxa")
length(LUX_bact_root_O_order_order)

Lux_Time_bact_soil_root_O_order_names=c("p:Acidobacteria.c:Acidobacteria.o:Acidobacteriales"="Acidobacteriales",           
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_3"="Subgroup_3",
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_6"="Subgroup_6",
                                        "p:Actinobacteria.c:Actinobacteria.o:Frankiales"="Frankiales",
                                        "p:Actinobacteria.c:Actinobacteria.o:Kineosporiales"="Kineosporiales",
                                        "p:Actinobacteria.c:Actinobacteria.o:Micromonosporales"="Micromonosporales",
                                        "p:Actinobacteria.c:Actinobacteria.o:Pseudonocardiales"="Pseudonocardiales",
                                        "p:Firmicutes.c:Bacilli.o:Bacillales"="Bacillales",
                                        "p:Proteobacteria.c:Betaproteobacteria.o:Burkholderiales"="Burkholderiales",
                                        "p:Proteobacteria.c:Deltaproteobacteria.o:Myxococcales"="Myxococcales",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Pseudomonadales"="Pseudomonadales",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Rhizobiales"="Rhizobiales",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Sphingomonadales"="Sphingomonadales",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Xanthomonadales"="Xanthomonadales",
                                        "p:Proteobacteria.c:Betaproteobacteria.Unknown"="Unclassified Betaproteobacteria",
                                        "p:Proteobacteria.Unknown.Unknown"="Unclassified Proteobacteria",
                                        "p:Verrucomicrobia.c:Spartobacteria.o:Chthoniobacterales"="Chthoniobacterales",
                                        "Other_Taxa" ="Other taxa")


length(Lux_Time_bact_soil_root_O_order_names)


Lux_Time_bact_soil_root_O_order_color=c("p:Acidobacteria.c:Acidobacteria.o:Acidobacteriales"="#1B9E77",           
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_3"="#11634A",
                                        "p:Acidobacteria.c:Acidobacteria.o:Subgroup_6"="#07281E",
                                        "p:Actinobacteria.c:Actinobacteria.o:Frankiales"="#D95F02", 
                                        "p:Actinobacteria.c:Actinobacteria.o:Kineosporiales"="#A34702",
                                        "p:Actinobacteria.c:Actinobacteria.o:Micromonosporales"="#6D3001",
                                        "p:Actinobacteria.c:Actinobacteria.o:Pseudonocardiales"="#361801",
                                        "p:Firmicutes.c:Bacilli.o:Bacillales"="#E7298A", 
                                        "p:Proteobacteria.c:Betaproteobacteria.o:Burkholderiales"="#66A61E",
                                        "p:Proteobacteria.c:Deltaproteobacteria.o:Myxococcales"="#59911A",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Pseudomonadales"="#4D7D16",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Rhizobiales"="#406813",
                                        "p:Proteobacteria.c:Alphaproteobacteria.o:Sphingomonadales"="#33530F",
                                        "p:Proteobacteria.c:Gammaproteobacteria.o:Xanthomonadales"="#263E0B",
                                        "p:Proteobacteria.c:Betaproteobacteria.Unknown"="#1A2A07",
                                        "p:Proteobacteria.Unknown.Unknown"="#0D1504",
                                        "p:Verrucomicrobia.c:Spartobacteria.o:Chthoniobacterales"="#735501",
                                        "Other_Taxa" ="grey")
length(Lux_Time_bact_soil_root_O_order_color)


(o_LUX__bact_soil_color2=ggplot(GLBRC018_OTU_LUX_bact_soil_core_facet.order.prop_otu2_sum_M2,
                                aes(x=(as.Date(raw_date,format="%m.%d.%Y")),y=value_std*100))+
    geom_area(aes( fill=factor(Real_P_C,levels = Lux_Time_bact_soil_O_order_order)),color="black")+
    theme_bw()+scale_x_date(limits = as.Date(c("3.19.2018","11.05.2018"),format="%m.%d.%Y"))+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = Lux_Time_bact_soil_root_O_order_color,labels=Lux_Time_bact_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))


(o_LUX__bact_root_color2=ggplot(GLBRC018_OTU_bact_LUX_root_core_facet.order.prop_otu2_sum_M2,aes(x=(as.Date(raw_date,format="%m.%d.%Y")),y=value_std*100))+
    geom_area(aes( fill=factor(Real_P_C,levels = LUX_bact_root_O_order_order)),color="black")+
    theme_bw()+scale_x_date(limits = as.Date(c("3.19.2018","11.05.2018"),format="%m.%d.%Y"))+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_blank(),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = Lux_Time_bact_soil_root_O_order_color,labels=Lux_Time_bact_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))

plot_grid(o_LUX__bact_root_color2,o_LUX__bact_soil_color2,nrow = 2,labels = c("a)","b)"),
          align = "v")

#1200x1100


#####Stack Soil Fungi Core Taxa plot####

#Need to turn the subsetted community into a phyloseq obj

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5[[4]],fill=="core")$otu,
                                                          GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil)
#I am getting a lot of Unknown phyla with sintax so I am going to use CONSTAX 


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl=phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl),
                                                        sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl),
                                                        tax_table(as.matrix(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)))



ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl)
#154
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl)/ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil)
#0.0606538
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl))
#2431168
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil))
#0.6966097

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl))
#9022
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl))
#2719
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl)
#349

#Calculated the relative abundance of the core community
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl_map=(sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl))
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl)
fung_LUX_soil_core_freq=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl_map%>%group_by(collectionDate)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)
fung_LUX_soil_core_freq$raw_date=str_replace_all(fung_LUX_soil_core_freq$collectionDate,"/",".")
#merge OTUs by site
GLBRC018_OTU_LUX_fung_soil_core_facet=merge_samples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_phyl, "collectionDate")
sample_names(GLBRC018_OTU_LUX_fung_soil_core_facet)     


#order Level
#combine the reads at Order level
get_taxa_unique(GLBRC018_OTU_LUX_fung_soil_core_facet, taxonomic.rank="Order")
#[1] "Unknown"               "Glomerales"            "Agaricales"            "Hypocreales"           "Mortierellales"        "Archaeorhizomycetales"
#[7] "Tremellales"           "Sordariales"           "Helotiales"            "Pezizales"             "Capnodiales"           "Paraglomerales"       
#[13] "Rhizophydiales"        "Sebacinales"           "Chaetothyriales"       "Pleosporales"          "Eurotiales"            "Chaetosphaeriales"    
#[19] "Dothideales"           "Filobasidiales"        "Thelephorales"         "Glomerellales"         "Polyporales"           "Gs07"                 
#[25] "Gs20"                  "Orbiliales"            "Venturiales"           "Cantharellales"        "Xylariales"  
(GLBRC018_OTU_LUX_fung_soil_core_facet.order<-tax_glom(GLBRC018_OTU_LUX_fung_soil_core_facet, taxrank="Order"))
#37



GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_LUX_fung_soil_core_facet.order))$Phylum,
                                                                 "Order"=data.frame(tax_table(GLBRC018_OTU_LUX_fung_soil_core_facet.order))$Order,
                                                                 "OTU"=taxa_names(GLBRC018_OTU_LUX_fung_soil_core_facet.order))

GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names_raw$P_C=with(GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names_raw,interaction(Phylum,Order))
#Let's take the top 10 taxa 

LUX_soil_fung_TopORDER = names(sort(taxa_sums(GLBRC018_OTU_LUX_fung_soil_core_facet.order), TRUE)[1:10])
LUX_soil_fung_TOP10_Order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(LUX_soil_fung_TopORDER, GLBRC018_OTU_LUX_fung_soil_core_facet.order)))$Phylum,
                                           "Order"=data.frame(tax_table(prune_taxa(LUX_soil_fung_TopORDER, GLBRC018_OTU_LUX_fung_soil_core_facet.order)))$Order,
                                           "OTU"=taxa_names(prune_taxa(LUX_soil_fung_TopORDER, GLBRC018_OTU_LUX_fung_soil_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names=merge(GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names_raw,LUX_soil_fung_TOP10_Order_Names, by="OTU", all.x = T)

GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names$Real_P_C=ifelse(is.na(GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names$Order.y),
                                                                  paste(GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names$Phylum.x,
                                                                        "Other_spp",sep="."),                    
                                                                  as.character(GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names$P_C)
)

#Transform the read counts to prop of total reads

GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_LUX_fung_soil_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop)))
GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2=merge(GLBRC018_OTU_LUX_fung_soil_core_facet.order_Names[,c("Real_P_C","OTU")],
                                                            GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

row.names(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2)=GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2$OTU
GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2$OTU=NULL
GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum=GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))
colSums(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum[,2:ncol(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum)])
GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_trunc=GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2[row.names(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2)%in%LUX_soil_fung_TopORDER,]

GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_trunc1=GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))

GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M=pivot_longer(data.frame(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_trunc1),
                                                                         !Real_P_C, names_to = "variable")
GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M$raw_date=str_replace(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M$variable,"X","")
sort(unique(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M$Real_P_C))

GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M2=merge(
  GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M,fung_LUX_soil_core_freq, by="raw_date")
GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M2=GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M2%>%mutate(value_std=value*samp_sum_mean)

LUX_fung_soil_O_order_order=sort(unique(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M2$Real_P_C))

LUX_fung_soil_O_order_names=c("Ascomycota.Archaeorhizomycetales"="Archaeorhizomycetales","Ascomycota.Chaetothyriales"="Chaetothyriales",
                              "Ascomycota.Helotiales"="Helotiales","Ascomycota.Hypocreales"="Hypocreales",
                              "Ascomycota.Unknown"="Unclassified Ascomycota",
                              "Basidiomycota.Agaricales"="Agaricales","Chytridiomycota.Rhizophydiales"="Rhizophydiales",
                              "Glomeromycota.Paraglomerales"="Paraglomerales",
                              "Mortierellomycota.Mortierellales"="Mortierellales", "Other_Taxa"="Other taxa","Unknown.Unknown"="Unclassified taxa")


LUX_fung_soil_O_order_color=c("Ascomycota.Archaeorhizomycetales"="#D95F02","Ascomycota.Chaetothyriales"="#A34702",
                              "Ascomycota.Helotiales"="#6D3001","Ascomycota.Hypocreales"="#361801",
                              "Ascomycota.Unknown"="#1B0C00",
                              "Basidiomycota.Agaricales"="#7570B3","Chytridiomycota.Rhizophydiales"="#1B9E77",
                              "Glomeromycota.Paraglomerales"="#E7298A",
                              "Mortierellomycota.Mortierellales"="#66A61E", "Other_Taxa"="darkgrey","Unknown.Unknown"="lightgrey")




(o_LUX__fung_soil_color=ggplot(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M2,aes(x=as.Date(raw_date,format="%m.%d.%Y"),y=value_std))+
    geom_area(aes( fill=factor(Real_P_C,levels = LUX_fung_soil_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_fill_manual(values = LUX_fung_soil_O_order_color,labels=LUX_fung_soil_O_order_names)+
    guides(fill=guide_legend(title="Order")))



#####Stack Root Fungi Core Taxa plot####

#Need to turn the subsetted community into a phyloseq obj

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5[[4]],fill=="core")$otu,
                                                          GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl=phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl),
                                                        sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl),
                                                        tax_table(as.matrix(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)))


ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl)
#60
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl)/ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)
#0.02806361
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl))
#929974
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root))
#0.659556

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl))
#9347
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl))
#223
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl)
#141

#Calculated the relative abundance of the core community
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl_map=(sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl))
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl_map$samp_sum=sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl)

fung_LUX_root_core_freq=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl_map%>%group_by(collectionDate)%>%summarise(samp_sum_mean=mean(samp_sum)/10000)
fung_LUX_root_core_freq$raw_date=str_replace_all(fung_LUX_root_core_freq$collectionDate,"/",".")
#merge OTUs by site
GLBRC018_OTU_LUX_fung_root_core_facet=merge_samples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_phyl, "collectionDate")
sample_names(GLBRC018_OTU_LUX_fung_root_core_facet)     


#order Level
#combine the reads at Order level
get_taxa_unique(GLBRC018_OTU_LUX_fung_root_core_facet, taxonomic.rank="Order")
#[1] "Glomerales"            "Hypocreales"           "Mortierellales"        "Helotiales"            "Sordariales"           "Unknown"              
#[7] "Chaetothyriales"       "Archaeorhizomycetales" "Agaricales"            "Sebacinales"           "Paraglomerales"        "Chaetosphaeriales"    
#[13] "Pleosporales" 
(GLBRC018_OTU_LUX_fung_root_core_facet.order<-tax_glom(GLBRC018_OTU_LUX_fung_root_core_facet, taxrank="Order"))
#17



GLBRC018_OTU_LUX_fung_root_core_facet.order_Names_raw=data.frame("Phylum"=data.frame(tax_table(GLBRC018_OTU_LUX_fung_root_core_facet.order))$Phylum,
                                                                 "Order"=data.frame(tax_table(GLBRC018_OTU_LUX_fung_root_core_facet.order))$Order,
                                                                 "OTU"=taxa_names(GLBRC018_OTU_LUX_fung_root_core_facet.order))

GLBRC018_OTU_LUX_fung_root_core_facet.order_Names_raw$P_C=with(GLBRC018_OTU_LUX_fung_root_core_facet.order_Names_raw,interaction(Phylum,Order))
#Let's take the top 10 taxa 

LUX_root_fung_TopORDER = names(sort(taxa_sums(GLBRC018_OTU_LUX_fung_root_core_facet.order), TRUE)[1:10])
LUX_root_fung_TOP10_Order_Names=data.frame("Phylum"=data.frame(tax_table(prune_taxa(LUX_root_fung_TopORDER, GLBRC018_OTU_LUX_fung_root_core_facet.order)))$Phylum,
                                           "Order"=data.frame(tax_table(prune_taxa(LUX_root_fung_TopORDER, GLBRC018_OTU_LUX_fung_root_core_facet.order)))$Order,
                                           "OTU"=taxa_names(prune_taxa(LUX_root_fung_TopORDER, GLBRC018_OTU_LUX_fung_root_core_facet.order)))



#Need to make an other section

GLBRC018_OTU_LUX_fung_root_core_facet.order_Names=merge(GLBRC018_OTU_LUX_fung_root_core_facet.order_Names_raw,LUX_root_fung_TOP10_Order_Names, by="OTU", all.x = T)

GLBRC018_OTU_LUX_fung_root_core_facet.order_Names$Real_P_C=ifelse(is.na(GLBRC018_OTU_LUX_fung_root_core_facet.order_Names$Order.y),
                                                                  paste(GLBRC018_OTU_LUX_fung_root_core_facet.order_Names$Phylum.x,
                                                                        "Other_spp",sep="."),                    
                                                                  as.character(GLBRC018_OTU_LUX_fung_root_core_facet.order_Names$P_C)
)

#Transform the read counts to prop of total reads

GLBRC018_OTU_LUX_fung_root_core_facet.order.prop=transform_sample_counts(GLBRC018_OTU_LUX_fung_root_core_facet.order, function(x)x/sum(x))



GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu=as.data.frame(t(otu_table(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop)))
GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2=merge(GLBRC018_OTU_LUX_fung_root_core_facet.order_Names[,c("Real_P_C","OTU")],
                                                            GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu,by.y = "row.names",by.x = "OTU")

row.names(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2)=GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2$OTU
GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2$OTU=NULL
GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum=GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2%>%group_by(Real_P_C)%>%
  summarise_all(~sum(.))
colSums(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum[,2:ncol(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum)])
GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_trunc=GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2[row.names(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2)%in%LUX_root_fung_TopORDER,]

GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_trunc1=GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_trunc %>% 
  bind_rows(summarise_all(., list(~if(is.numeric(.)) 1-sum(.) else "Other_Taxa")))

GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M=pivot_longer(data.frame(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_trunc1),
                                                                         !Real_P_C,names_to = "variable")
GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M$raw_date=str_replace(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M$variable,"X","")
sort(unique(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M$Real_P_C))

GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2=merge(
  GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M,fung_LUX_root_core_freq, by="raw_date")
GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2=GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2%>%mutate(value_std=value*samp_sum_mean)

#Ascomycota.Unknown is actually a Sordariomycetes

GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2$Real_P_C=str_replace(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2$Real_P_C,
                                                                                  "Ascomycota.Unknown","Ascomycota.Sod.Unknown")
sort(unique(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2$Real_P_C))

LUX_fung_root_O_order_order=c("Ascomycota.Chaetosphaeriales",
                              "Ascomycota.Chaetothyriales",
                              "Ascomycota.Helotiales",
                              "Ascomycota.Hypocreales",
                              "Ascomycota.Pleosporales",
                              "Ascomycota.Sordariales",
                              "Ascomycota.Sod.Unknown",
                              "Basidiomycota.Agaricales",
                              "Basidiomycota.Unknown",
                              "Glomeromycota.Glomerales",
                              "Other_Taxa")

LUX_fung_soil_root_O_order_names=c("Ascomycota.Archaeorhizomycetales"="Archaeorhizomycetales",
                                   "Ascomycota.Chaetosphaeriales"="Chaetosphaeriales",
                                   "Ascomycota.Chaetothyriales"="Chaetothyriales",
                                   "Ascomycota.Helotiales"="Helotiales",
                                   "Ascomycota.Hypocreales"="Hypocreales",
                                   "Ascomycota.Pleosporales"="Pleosporales",
                                   "Ascomycota.Sordariales"="Sordariales",
                                   "Ascomycota.Sod.Unknown"="Unclassified Sordariomycetes",
                                   "Ascomycota.Unknown"="Unclassified Ascomycota",
                                   "Basidiomycota.Agaricales"="Agaricales",
                                   "Basidiomycota.Unknown"="Unclassified Agaricomycetes",
                                   "Chytridiomycota.Rhizophydiales"="Rhizophydiales",
                                   "Glomeromycota.Glomerales"="Glomerales",
                                   "Glomeromycota.Paraglomerales"="Paraglomerales",
                                   "Mortierellomycota.Mortierellales"="Mortierellales", 
                                   "Other_Taxa"="Other taxa",
                                   "Unknown.Unknown"="Unclassified taxa")
length(LUX_fung_soil_root_O_order_names)

LUX_fung_soil_root_O_order_color=c("Ascomycota.Archaeorhizomycetales"="#D95F02",
                                   "Ascomycota.Chaetosphaeriales"="#BE5302",
                                   "Ascomycota.Chaetothyriales"="#A34702",
                                   "Ascomycota.Helotiales"="#883B01",
                                   "Ascomycota.Hypocreales"="#6D3001",
                                   "Ascomycota.Pleosporales"="#512401",
                                   "Ascomycota.Sordariales"="#361801",
                                   "Ascomycota.Sod.Unknown"="#1B0C00",
                                   "Ascomycota.Unknown"="black",
                                   "Basidiomycota.Agaricales"="#7570B3",
                                   "Basidiomycota.Unknown"="#3A3859",
                                   "Chytridiomycota.Rhizophydiales"="#1B9E77",
                                   "Glomeromycota.Glomerales"="#E7298A",
                                   "Glomeromycota.Paraglomerales"="#731545",
                                   "Mortierellomycota.Mortierellales"="#66A61E", 
                                   "Other_Taxa"="darkgrey",
                                   "Unknown.Unknown"="lightgrey")



length(LUX_fung_soil_root_O_order_color)


(o_LUX__fung_root_color=ggplot(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2,aes(x=as.Date(raw_date,format="%m.%d.%Y"),y=value_std))+
    geom_area(aes( fill=factor(Real_P_C,levels = LUX_fung_root_O_order_order)),color="black")+
    theme_bw()+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_fill_manual(values = LUX_fung_soil_root_O_order_color,labels=LUX_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))


(o_LUX__fung_soil_color2=ggplot(GLBRC018_OTU_LUX_fung_soil_core_facet.order.prop_otu2_sum_M2,aes(x=as.Date(raw_date,format="%m.%d.%Y"),y=value_std*100))+
    geom_area(aes( fill=factor(Real_P_C,levels = LUX_fung_soil_O_order_order)),color="black")+
    theme_bw()+scale_x_date(limits = as.Date(c("3.19.2018","11.05.2018"),format="%m.%d.%Y"))+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = LUX_fung_soil_root_O_order_color,labels=LUX_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))

(o_LUX__fung_root_color2=ggplot(GLBRC018_OTU_LUX_fung_root_core_facet.order.prop_otu2_sum_M2,aes(x=as.Date(raw_date,format="%m.%d.%Y"),y=value_std*100))+
    geom_area(aes( fill=factor(Real_P_C,levels = LUX_fung_root_O_order_order)),color="black")+
    theme_bw()+scale_x_date(limits = as.Date(c("3.19.2018","11.05.2018"),format="%m.%d.%Y"))+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_blank(),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Percentage of reads")+
    scale_fill_manual(values = LUX_fung_soil_root_O_order_color,labels=LUX_fung_soil_root_O_order_names)+
    guides(fill=guide_legend(title="Order")))

plot_grid(o_LUX__fung_root_color2,o_LUX__fung_soil_color2,nrow = 2,labels = c("a)","b)"),
          align = "v")




plot_grid(o_LUX__bact_root_color2,o_LUX__fung_root_color2,o_LUX__bact_soil_color2,o_LUX__fung_soil_color2,nrow = 2,labels = c("a)","b)","c)","d)"),
          align = "v")

ggsave(plot_grid(o_LUX__bact_root_color2,o_LUX__fung_root_color2,o_LUX__bact_soil_color2,o_LUX__fung_soil_color2,
                 nrow = 2,labels = c("a)","b)","c)","d)"),label_size = 48,label_x = c(-0.02,-0.02,-0.02,-0.02),
                 align = "v"), filename = "Lux_temporal_bacteria_fungi_root_soil_core_stack_order_p.png",path = here::here("Manuscript","Core_comm_figs"),width = 30,height = 15)

ggsave(plot_grid(o_LUX__bact_root_color2,o_LUX__fung_root_color2,o_LUX__bact_soil_color2,o_LUX__fung_soil_color2,
                 nrow = 2,labels = c("a)","b)","c)","d)"),label_size = 48,label_x = c(-0.02,-0.02,-0.02,-0.02),
                 align = "v"), filename = "Lux_temporal_bacteria_fungi_root_soil_core_stack_order_p.svg",path = here::here("Manuscript","Core_comm_figs"),width = 30,height = 15)


####LUX metacoder Soil versus Root Bacteria  Core Taxa####



#Root Core taxa in Soil community community


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_core_5[[4]],fill=="core")$otu,
                                                     GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)

ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT)
#66
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT)/ntaxa(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
#0.002751032

sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT))
#1176000
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil))
#0.2375758

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT))
#5834
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT))
#991
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT)
#495

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC=data.frame(tax_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT))
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC$Species=
  with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC,
       paste(Genus,row.names(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC),sep = ""))


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC=as.matrix(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC|>mutate_all(str_remove,"^(.+):"))
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT,function(x)(x/10000))),
           tax_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC),
           sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT))
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC_phyl)
#66
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC=parse_phyloseq(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_taxC_phyl)
print(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC)


set.seed(2021) # This makes the plot appear the same each time it is run 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC |> 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations




GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID,
                   groups = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Project,
                   out_names = "Mean_across_samples")


#Need the mean taxon abundance
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID)




GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                     log2_median_ratio))




GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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
subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam<0)

#Known Taxa only 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                        log2_median_ratio))




(Bact_LUX_core_ROOT_metaHeat_P=GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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
subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table,log2_median_ratio_Fam_KN<0)




#Soil Core taxa in Root community


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL=prune_taxa(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_core_5[[4]],fill=="core")$otu,
                                                     GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)

ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL)
#74
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL)/ntaxa(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
#0.00308449

sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL))
#1293608
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL))/sum(otu_table(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil))
#0.2613349

max(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL))
#4093
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL))
#1416
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL)
#495

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC=data.frame(tax_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL))
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC$Species=
  with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC,
       paste(Genus,row.names(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC),sep = ""))


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC=as.matrix(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC|>mutate_all(str_remove,"^(.+):"))
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL,function(x)(x/10000))),
           tax_table(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC),
           sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL))
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC_phyl)
#74
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC=parse_phyloseq(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_taxC_phyl)
print(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC)


set.seed(2021) # This makes the plot appear the same each time it is run 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC |> 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations




GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID,
                   groups = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Project,
                   out_names = "Mean_across_samples")



#Need the mean taxon abundance
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID)




GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                     log2_median_ratio))




GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
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
subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam>0)

#Known Taxa only 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$SampleID, # What columns of sample data to use
                 groups = GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                        log2_median_ratio))





(Bact_LUX_core_ROOT_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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

(Bact_LUX_core_SOIL_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
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

plot_grid(Bact_LUX_core_ROOT_metaHeat_P2,Bact_LUX_core_SOIL_metaHeat_P2,
          align = "h")





#####LUX metacoder Soil versus Root Fungi Core Taxa####


GLBRC018_OTU_fung_MMPRNT_LUX_G5=phyloseq(otu_table(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5),
                                         sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5),
                                         tax_table(as.matrix(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)))








#Root Core taxa in Soil community community


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_core_5[[4]],fill=="core")$otu,
                                                     GLBRC018_OTU_fung_MMPRNT_LUX_G5)

ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT)
#60
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT)/ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5)
#0.02203452
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT))
#1842455
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5))
#0.3760112

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT))
#9347
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT))
#223
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT)
#490

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC=data.frame(tax_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT))
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC$Species=
  with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC,
       paste(Species,row.names(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC),sep = ""))



GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT,function(x)(x/10000))),
           tax_table(as.matrix(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC)),
           sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT))
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC_phyl)
#60
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC=parse_phyloseq(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_taxC_phyl)
print(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC)



GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil=
  factor(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil,
         levels = c("Root","Soil"),ordered = T)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil),]

set.seed(2021) # This makes the plot appear the same each time it is run 
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC |> 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations




GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq,
                   groups = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$project,
                   out_names = "Mean_across_samples")



#Need the mean taxon abundance
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq)




GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                     log2_median_ratio))




GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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
subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam<0)

#Known Taxa only 
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                        log2_median_ratio))




(Fung_LUX_core_ROOT_metaHeat_P=GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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
subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC$data$Fam_KN_diff_table,log2_median_ratio_Fam_KN<0)




#Soil Core taxa in Root community


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL=prune_taxa(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_core_5[[4]],fill=="core")$otu,
                                                     GLBRC018_OTU_fung_MMPRNT_LUX_G5)

ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL)
#154
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL)/ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5)
# 0.05655527

sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL))
#3072923
sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL))/sum(otu_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5))
# 0.6271271

max(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL))
#9022

min(sample_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL))
#243
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL)
#490

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC=data.frame(tax_table(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL))
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC$Species=
  with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC,
       paste(Species,row.names(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC),sep = ""))



GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC_phyl=
  phyloseq(otu_table(transform_sample_counts(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL,function(x)(x/10000))),
           tax_table(as.matrix(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC)),
           sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL))
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC_phyl)
#154
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC=parse_phyloseq(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_taxC_phyl)
print(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil=
  factor(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil,
         levels = c("Root","Soil"),ordered = T)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data[order(
    GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil),]


set.seed(2021) # This makes the plot appear the same each time it is run 
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC |> 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$All_tax_abund <- 
  calc_taxon_abund(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC, "otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq,
                   groups = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$project,
                   out_names = "Mean_across_samples")


#Need the mean taxon abundance
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$All_tax_abund$Mean_across_samples/length(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq)




GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  compare_groups(data = "Fam_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value_Fam<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table$log2_median_ratio_Fam= 
  with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table, ifelse(wilcox_p_value_Fam > 0.05,0,
                                                                                     log2_median_ratio))



GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
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
subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_diff_table,log2_median_ratio_Fam>0)

#Known Taxa only 
GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_tax_abund<- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  calc_taxon_abund("otu_table",
                   cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table <- 
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
  filter_taxa(taxon_ranks == "Family", supertaxa = TRUE)|>
  filter_taxa(taxon_names!="Unknown", supertaxa = TRUE)|>
  compare_groups(data = "Fam_KN_tax_abund",
                 cols = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$sampleID_seq, # What columns of sample data to use
                 groups = GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$sample_data$Root_soil) # What category each sample is assigned to
print(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value_Fam_KN<-
  p.adjust(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$wilcox_p_value,
           method = "fdr")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table$log2_median_ratio_Fam_KN= 
  with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table, ifelse(wilcox_p_value_Fam_KN > 0.05,0,
                                                                                        log2_median_ratio))

summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC$data$Fam_KN_diff_table)



(Fung_LUX_core_ROOT_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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

(Fung_LUX_core_SOIL_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
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

plot_grid(Fung_LUX_core_ROOT_metaHeat_P2,Fung_LUX_core_SOIL_metaHeat_P2,
          align = "h")


#####Figures metacoder Soil versus Root Bacteria Fungi Core Taxa#####

(Bact_LUX_core_ROOT_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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

(Bact_LUX_core_SOIL_metaHeat_P2=GLBRC018_OTU_bact_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
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


(Fung_LUX_core_ROOT_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_ROOT_phyl_MC|>
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

(Fung_LUX_core_SOIL_metaHeat_P2=GLBRC018_OTU_fung_MMPRNT_LUX_G5_CORE_SOIL_phyl_MC|>
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
          get_title(Root_90_title),Bact_LUX_core_ROOT_metaHeat_P2,Fung_LUX_core_ROOT_metaHeat_P2,
          get_title(Soil_90_title),Bact_LUX_core_SOIL_metaHeat_P2,Fung_LUX_core_SOIL_metaHeat_P2,
          rel_widths = c(0.1,1,1),rel_heights = c(0.1,1,1),ncol = 3)



ggsave(plot_grid(ggplot()+theme_void(),get_title(Bact_title),get_title(Fung_title),
                 get_title(Root_90_title),Bact_LUX_core_ROOT_metaHeat_P2,Fung_LUX_core_ROOT_metaHeat_P2,
                 get_title(Soil_90_title),Bact_LUX_core_SOIL_metaHeat_P2,Fung_LUX_core_SOIL_metaHeat_P2,
                 rel_widths = c(0.1,1,1),rel_heights = c(0.2,1,1),ncol = 3), 
       filename = "LUX_Heat_tree_Bact_Fung_Root_Soil_core_comm.png",path = here::here("Manuscript","Core_comm_figs"),width = 15,height = 10)










######Lux Arbor CORE FunGuild Fig#### 


#Extract OTU Table
GLBRC018_OTU_LUX_fung_soil_core_facet_otus=data.frame(t(otu_table(GLBRC018_OTU_LUX_fung_soil_core_facet)))
GLBRC018_OTU_LUX_fung_soil_core_facet_otus[1:10,1:6]

#Combine OTU and FunGuild


GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG=merge(GLBRC018_OTU_LUX_fung_soil_core_facet_otus,
                                                      GLBRC_CON_FunGuild_raw_F,
                                                      by.x = "row.names", by.y = "OTU",
                                                      all.x = T)

summary(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG)
unique(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG$simp_guild)
GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG[is.na(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG)]="Unknown"
#Simplified Guild 

GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum=GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG%>%
  group_by(simp_guild)%>%summarise(across(starts_with("X"),~sum(.)))

GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M=GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum%>%
  pivot_longer(!simp_guild, names_to = "X.date", values_to = "count")
GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M=GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M%>%
  group_by(X.date)%>%mutate(prop=count/sum(count))

GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M$raw_date=str_replace(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M$X.date,"X","")
sort(unique(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M$simp_guild))

GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2=merge(
  GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M,fung_LUX_soil_core_freq, by="raw_date")

dim(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2)
#210   8
head(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2)



GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2=GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2%>%mutate(value_std=prop*samp_sum_mean)
sort(unique(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2$simp_guild))


#Extract OTU Table
GLBRC018_OTU_LUX_fung_root_core_facet_otus=data.frame(t(otu_table(GLBRC018_OTU_LUX_fung_root_core_facet)))
GLBRC018_OTU_LUX_fung_root_core_facet_otus[1:10,1:6]

#Combine OTU and FunGuild


GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG=merge(GLBRC018_OTU_LUX_fung_root_core_facet_otus,
                                                      GLBRC_CON_FunGuild_raw_F,
                                                      by.x = "row.names", by.y = "OTU",
                                                      all.x = T)

summary(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG)
unique(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG$simp_guild)
GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG[is.na(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG)]="Unknown"
#Simplified Guild 

GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum=GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG%>%
  group_by(simp_guild)%>%summarise(across(starts_with("X"),~sum(.)))

GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M=GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum%>%
  pivot_longer(!simp_guild, names_to = "X.date", values_to = "count")
GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M=GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M%>%
  group_by(X.date)%>%mutate(prop=count/sum(count))

GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M$raw_date=str_replace(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M$X.date,"X","")
sort(unique(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M$simp_guild))

GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2=merge(
  GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M,fung_LUX_root_core_freq, by="raw_date")

dim(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2)
#42  7
head(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2)



GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2=GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2%>%mutate(value_std=prop*samp_sum_mean)
sort(unique(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2$simp_guild))


LUX_fung_soil_FGuild_order=c("Arbuscular Mycorrhizal","Ectomycorrhizal","Endophyte",
                             "Animal Pathogen","Plant Pathogen",
                             "Multiple Saprotroph","Undefined Saprotroph","Soil Saprotroph","Wood Saprotroph",
                             "Symbiotroph-Pathogen",
                             "Symbiotroph-Saprotroph",
                             "Pathogen-Saprotroph",
                             "Symbiotroph-Pathogen-Saprotroph",
                             "Unknown")


LUX_fung_soil_FGuild_color=c("Arbuscular Mycorrhizal"="#1B9E77","Ectomycorrhizal"="#11634A","Endophyte"="#07281E",
                             "Animal Pathogen"="#D95F02","Plant Pathogen"="#6D3001",
                             "Multiple Saprotroph"="#7570B3","Undefined Saprotroph"="#585486",
                             "Soil Saprotroph"="#3A3859","Wood Saprotroph"="#1D1C2D",
                             "Symbiotroph-Pathogen"="#E7298A",
                             "Symbiotroph-Saprotroph"="#66A61E",
                             "Pathogen-Saprotroph"="#E6AB02",
                             "Symbiotroph-Pathogen-Saprotroph"="#A6761D",
                             "Unknown"="lightgrey")




(FunGuild_LUX__fung_soil_color2=ggplot(GLBRC018_OTU_LUX_fung_soil_core_facet_otus_FunG_sum_M2,aes(x=as.Date(raw_date,format="%m.%d.%Y"),
                                                                                                  y=(value_std*100)))+
    geom_area(aes( fill=factor(simp_guild,levels = LUX_fung_soil_FGuild_order)),color="black")+
    theme_bw()+scale_x_date(limits = as.Date(c("3.19.2018","11.05.2018"),format="%m.%d.%Y"))+
    scale_y_continuous(name = "Soil\nPercentage of reads",limits = c(0,80))+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+
    scale_fill_manual(values = LUX_fung_soil_FGuild_color)+
    guides(fill=guide_legend(title="FunGuild")))



LUX_fung_root_FGuild_order=c("Arbuscular Mycorrhizal","Endophyte",
                             "Multiple Saprotroph","Undefined Saprotroph",
                             "Symbiotroph-Saprotroph",
                             "Pathogen-Saprotroph",
                             "Unknown")


LUX_fung_root_FGuild_color=c("Arbuscular Mycorrhizal"="#1B9E77","Endophyte"="#0E4F3B",
                             "Multiple Saprotroph"="#7570B3","Undefined Saprotroph"="#3A3859",
                             "Symbiotroph-Saprotroph"="#66A61E",
                             "Pathogen-Saprotroph"="#E6AB02",
                             "Unknown"="lightgrey")




(FunGuild_LUX__fung_root_color2=ggplot(GLBRC018_OTU_LUX_fung_root_core_facet_otus_FunG_sum_M2,aes(x=as.Date(raw_date,format="%m.%d.%Y"),
                                                                                                  y=(value_std*100)))+
    geom_area(aes( fill=factor(simp_guild,levels = LUX_fung_root_FGuild_order)),color="black")+
    theme_bw()+scale_x_date(limits = as.Date(c("3.19.2018","11.05.2018"),format="%m.%d.%Y"),labels = NULL)+
    scale_y_continuous(name = "Root\nPercentage of reads",limits = c(0,80))+
    theme(axis.text.y=element_text(size=28),axis.text.x=element_text(size=30),
          axis.title=element_text(size=36),panel.grid.major=element_blank(),legend.text = element_text(size=24), legend.title = element_text(size=28),
          panel.grid.minor=element_blank())+xlab(NULL)+
    scale_fill_manual(values = LUX_fung_root_FGuild_color)+
    guides(fill=guide_legend(title="FunGuild")))



plot_grid(FunGuild_LUX__fung_root_color2,FunGuild_LUX__fung_soil_color2,
          nrow = 2,labels = c("a)","b)"),axis = "lr",
          align = "v",label_size = 48,label_x = c(-0.02,-0.02))

ggsave(plot_grid(FunGuild_LUX__fung_root_color2,FunGuild_LUX__fung_soil_color2,
                 nrow = 2,labels = c("a)","b)"),axis = "lr",
                 align = "v",label_size = 48,label_x = c(-0.02,-0.02)), 
       filename = "Lux_temporal_fungi_FunGuils_root_soil_core_stack_order_sep_color.png",path = here::here("Manuscript","Core_comm_figs"),width = 16,height = 15)

