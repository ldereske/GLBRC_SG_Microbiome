
## ---------------------------
##
## Script name: Community characterization of the switchgrass microbiomes of 
## the Marginal Lands Experiment (MLE)
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
library(missForest)
library(stringr)
library(grid)
library(gdm)
library(eulerr)
library(ggplotify)

#Load code for calculating importance and p-values for GDMs with 3 or less factors
source(here::here("R_functions","gdm.varImp_MOD.R"))


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


#####MLE community composition####


#MLE is defined as the July overlap in root and soil sampling

unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"&siteID=="LUX"))$collectionDate)
#"5/29/2018" "9/17/2018" "8/20/2018" "7/30/2018" "6/25/2018" "10/3/2018"
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"&siteID=="LC"))$collectionDate)
#"7/10/2018"


GLBRC018_OTU_bact_MMPRNT_MLE=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,siteID!="LUX"|collectionDate=="7/30/2018")
nsamples(GLBRC018_OTU_bact_MMPRNT_MLE)
#513
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE,Root_soil=="Root"&siteID=="LUX"))$collectionDate)

GLBRC018_OTU_bact_MMPRNT_MLE=subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE,siteID!="LC"|collectionDate=="7/10/2018")
nsamples(GLBRC018_OTU_bact_MMPRNT_MLE)
#323
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE,Root_soil=="Root"&siteID=="LC"))$collectionDate)

#Let's use only G5 since there is overlap in sampling between roots and soil

GLBRC018_OTU_bact_MMPRNT_MLE_G5=subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE,plotType=="G5")
nsamples(GLBRC018_OTU_bact_MMPRNT_MLE_G5)
#227
rm(GLBRC018_OTU_bact_MMPRNT_MLE)
rm(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)

#I am going to define MLE as the July overlap in root and soil sampling

unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&siteID=="LUX"))$collectionDate)
#"5/29/2018" "9/17/2018" "8/20/2018" "7/30/2018" "6/25/2018" "10/3/2018"
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&siteID=="LC"))$collectionDate)
#"7/10/2018"
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#1086

GLBRC018_OTU_fung_MMPRNT_MLE=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,siteID!="LUX"|collectionDate=="7/30/2018")
nsamples(GLBRC018_OTU_fung_MMPRNT_MLE)
#501
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE,Root_soil=="Root"&siteID=="LUX"))$collectionDate)

GLBRC018_OTU_fung_MMPRNT_MLE=subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE,siteID!="LC"|collectionDate=="7/10/2018")
nsamples(GLBRC018_OTU_fung_MMPRNT_MLE)
#319
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE,Root_soil=="Root"&siteID=="LC"))$collectionDate)

#Let's use only G5 since there is overlap in sampling between roots and soil

GLBRC018_OTU_fung_MMPRNT_MLE_G5=subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE,plotType=="G5")
nsamples(GLBRC018_OTU_fung_MMPRNT_MLE_G5)
#223
rm(GLBRC018_OTU_fung_MMPRNT_MLE)
rm(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)



#Separate in organs



#Roots only 
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root=subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE_G5,Root_soil=="Root")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root=
  prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root)>0,GLBRC018_OTU_bact_MMPRNT_MLE_G5_root)
nsamples(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root)
#113

GLBRC018_OTU_fung_MMPRNT_MLE_G5_root=subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE_G5,Root_soil=="Root")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root=
  prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root)>0,GLBRC018_OTU_fung_MMPRNT_MLE_G5_root)
nsamples(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root)
#111

#NMDS
set.seed(2021)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root,method = "NMDS")
#*** Solution reached
#0.1226669   

set.seed(2021)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root,method = "NMDS")
#*** Solution reached
#0.2188881  



#Soil only 
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil=subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE_G5,Root_soil=="Soil")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil=
  prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil)>0,GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil)
nsamples(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil)
#114

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE_G5,Root_soil=="Soil")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil=
  prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)>0,GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)
nsamples(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)
#112

#NMDS
set.seed(2021)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil,method = "NMDS")
#*** Solution reached
#0.1199976    

set.seed(2021)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil,method = "NMDS")
#*** Solution reached
#0.1712595   




GLBRC018_OTU_bact_MMPRNT_MLE_G5_dis=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_MLE_G5,method = "bray")



GLBRC018_OTU_bact_MMPRNT_MLE_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_MLE_G5)



GLBRC018_OTU_fung_MMPRNT_MLE_G5_dis=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_MLE_G5,method = "bray")


GLBRC018_OTU_fung_MMPRNT_MLE_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_MLE_G5)


#####All Site Abiotic Metadata ####
#Metadata for community predictions


MMPRNT_2018.metadata_mapp_coverage= read.csv(here::here("Publish_data","MMPRNT_2018.metadata_mapp_coverage.csv"),
                                             header = T)

dim(MMPRNT_2018.metadata_mapp_coverage)
#870  73




MMPRNT_2018.metadata_mapp_coverage_MLE=subset(MMPRNT_2018.metadata_mapp_coverage,siteID!="LUX"|collectionDate=="7/30/2018")
nrow(MMPRNT_2018.metadata_mapp_coverage_MLE)
#486


MMPRNT_2018.metadata_mapp_coverage_MLE=subset(MMPRNT_2018.metadata_mapp_coverage_MLE,siteID!="LC"|collectionDate=="7/10/2018")
nrow(MMPRNT_2018.metadata_mapp_coverage_MLE)
#294

#matching to the actual sequence samples

MMPRNT_2018.metadata_mapp_coverage_MLE_G5=subset(MMPRNT_2018.metadata_mapp_coverage_MLE,plotType=="G5")
nrow(MMPRNT_2018.metadata_mapp_coverage_MLE_G5)
#114
colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_MLE_G5))


#Remove the redundant and not useful factors
MMPRNT_2018.metadata_mapp_MLE_G5=MMPRNT_2018.metadata_mapp_coverage_MLE_G5
MMPRNT_2018.metadata_mapp_MLE_G5[,c("event_age","pH_raw_MMPRNT","pH_MMPRNT_subplot","rain_mm_tot_sum")]=NULL

#Check number of NAs
colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_MLE_G5))

#Correlations in the data
cor(MMPRNT_2018.metadata_mapp_coverage_MLE_G5[,c("percent_soil_moisture_dry_weight","vwc_avg_mean","dry_matter_yield_mg_ha_mean",
                                                       "total.shoot.dry.weight","plant.height","specific.leaf.area","core_root_mass_subplot",
                                                       "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na","soil_temp_1_avg_mean","airtc_avg_mean","past_7d_rain",
                                                       "ph","p_ppm","k_ppm","ca_ppm","mg_ppm","ugN_fixed_gdrysoil_day","TOC","TON","TOC_TON_ratio"
)],use = "na.or.complete")






ggpairs(MMPRNT_2018.metadata_mapp_coverage_MLE_G5, columns = c("percent_soil_moisture_dry_weight","dry_matter_yield_mg_ha_mean","total.shoot.dry.weight","core_root_mass_subplot","ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na","past_7d_rain","ph","k_ppm","ca_ppm","TOC","TON","TOC_TON_ratio"), 
        upper = list(continuous = wrap("cor",size = 3)),
        lower = list(continuous = wrap("smooth",
                                       alpha = 0.3,
                                       size = 0.1)))




#Non-correlated Metadata Extraction 





MMPRNT_2018.metadata_mapp_MLE_G5_lim=
  data.frame(MMPRNT_2018.metadata_mapp_MLE_G5[,c("SampleID","siteID","percent_soil_moisture_dry_weight","dry_matter_yield_mg_ha_mean",
                                                       "total.shoot.dry.weight","core_root_mass_subplot","ugN_NH4_g_dry_soil_na",
                                                       "ugN_NO3_g_dry_soil_na","past_7d_rain","ph","k_ppm","ca_ppm","TOC","TON","TOC_TON_ratio")])
#Check number of NAs
colSums(is.na(MMPRNT_2018.metadata_mapp_MLE_G5_lim))

MMPRNT_2018.metadata_mapp_MLE_G5_lim|>
  group_by(siteID)|>
  summarise(total.shoot.dry.weight=sum(is.na(total.shoot.dry.weight)),
            ugN_NH4_g_dry_soil_na=sum(is.na(ugN_NH4_g_dry_soil_na)),
            ugN_NO3_g_dry_soil_na=sum(is.na(ugN_NO3_g_dry_soil_na)),
            TOC=sum(is.na(TOC)),
            TON=sum(is.na(TON)),
            TOC_TON_ratio=sum(is.na(TOC_TON_ratio)),
            percent_soil_moisture_dry_weight=sum(is.na(percent_soil_moisture_dry_weight)),
            core_root_mass_subplot=sum(is.na(core_root_mass_subplot)),
            past_7d_rain =sum(is.na(past_7d_rain)),
            ph=sum(is.na(ph)),
            k_ppm=sum(is.na(k_ppm)),
            ca_ppm=sum(is.na(ca_ppm)),
            k_ppm=sum(is.na(k_ppm)) )




#Impute the missing data
summary(MMPRNT_2018.metadata_mapp_MLE_G5_lim)
row.names(MMPRNT_2018.metadata_mapp_MLE_G5_lim)=MMPRNT_2018.metadata_mapp_MLE_G5_lim$SampleID

MMPRNT_2018.metadata_mapp_MLE_G5_lim$siteID=as.factor(MMPRNT_2018.metadata_mapp_MLE_G5_lim$siteID)

set.seed(2022)
MMPRNT_2018.metadata_mapp_MLE_G5_lim.imp <- missForest(MMPRNT_2018.metadata_mapp_MLE_G5_lim[,-1])


#check imputation error
MMPRNT_2018.metadata_mapp_MLE_G5_lim.imp$OOBerror
#     NRMSE        PFC 
#0.05325305 0.00000000  

MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp=MMPRNT_2018.metadata_mapp_MLE_G5_lim.imp$ximp
head(MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp)

MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp$sampleID_long=
  paste("MMPRNT",str_sub(row.names(MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp),start = 7,end = 11),sep = "-")

rm(list=c("MMPRNT_2018.metadata_mapp_MLE_G5_lim","MMPRNT_2018.metadata_mapp_coverage_MLE_G5",
          "MMPRNT_2018.metadata_mapp_coverage_MLE","MMPRNT_2018.metadata_mapp_coverage",
          "MMPRNT_2018.metadata_mapp_MLE_G5"))


#GDM MLE Bacterial Roots####



GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data=sample_data(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root)[,c("sampleID_long","SampleID","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data)
#113
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data)

colnames(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data)[colnames(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata=
  data.frame(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic=merge(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic)
#113  16



GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root,method = "bray")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF),GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic)
nrow(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic)

GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   61.339




#Remove factors with less than 0.01
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","past_7d_rain","k_ppm","ca_ppm",
                                                                       "ugN_NO3_g_dry_soil_na","core_root_mass_subplot",
                                                                       "ugN_NH4_g_dry_soil_na","TOC_TON_ratio",
                                                                       "dry_matter_yield_mg_ha_mean","TON")]
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub1_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub1_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub1_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   61.339

gdm.varImp(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub1_TabForm, geo = T, predSelect = T,nPerm = 100)

#Final set of predictors returned: 
#Geographic
#ph
#past_7d_rain
#k_ppm


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub=
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "past_7d_rain","ph","k_ppm")]
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   59.682

gdm.varImp(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm, geo = T, nPerm = 100)
#Model deviance                     92.283
#Percent deviance explained         59.682
#Model p-value                       0.000
#Fitted permutations               100.000



bact_R_meta_varSet <- vector("list", 2)

names(bact_R_meta_varSet)=c("MET","soil")
bact_R_meta_varSet$MET=c("past_7d_rain")
bact_R_meta_varSet$soil=c("ph","k_ppm")
summary(bact_R_meta_varSet)
#Variance partitioning
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm,varSets = bact_R_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_x_long=GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_x_long,
        GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_y_long)






VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")





bact_MLE_root_meta_color=c("past_7d_rain"="#4575B4",
                                 "ph"="#7C001D",
                                 "k_ppm"="#F46D43",
                                 "Geographic"="black")



bact_MLE_root_meta_names=c("past_7d_rain"="Seven day rain",
                                 "ph"="Soil pH",
                                 "k_ppm"="Soil K",
                                 "Geographic"="Spatial factor")



ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',linewidth=2)+
  scale_color_manual(values = bact_MLE_root_meta_color,name=NULL,labels=bact_MLE_root_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")




#GDM MLE Bacterial Soils####




GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data=sample_data(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil)[,c("sampleID_long","SampleID","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data)
#114
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data)

colnames(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data)[colnames(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata=
  data.frame(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic=merge(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)
#114  16



GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil,method = "bray")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF),GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)
nrow(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)

GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   82.883



#Remove factors with less than 0.01
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord","ph","ca_ppm",
                                                                       "TON","past_7d_rain","percent_soil_moisture_dry_weight",
                                                                       "k_ppm","ugN_NO3_g_dry_soil_na","dry_matter_yield_mg_ha_mean",
                                                                       "TOC_TON_ratio")]
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub1_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub1_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub1_gdmTabMod)
#82.781
gdm.varImp(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub1_TabForm, geo = T, predSelect = T, nPerm = 100)
#Final set of predictors returned: 
#Geographic
#ph
#ca_ppm
#TON
#past_7d_rain
#k_ppm






GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub=
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "past_7d_rain","ph","k_ppm","ca_ppm","TON")]
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   82.465

gdm.varImp(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm, geo = T, nPerm = 100)
#Model deviance                     56.484
#Percent deviance explained         82.465
#Model p-value                       0.000
#Fitted permutations               100.000



bact_S_meta_varSet <- vector("list", 2)

names(bact_S_meta_varSet)=c("MET","soil")
bact_S_meta_varSet$MET=c("past_7d_rain")
bact_S_meta_varSet$soil=c("ph","k_ppm","ca_ppm","TON")
summary(bact_S_meta_varSet)
#Variance partitioning
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm,
                                                                                   varSet = bact_S_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART)

GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_x_long=GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_x_long,
        GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_long)





VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")


bact_MLE_soil_meta_color=c("past_7d_rain"="#4575B4",
                                 "ph"="#7C001D",
                                 "k_ppm"="#F46D43",
                                 "ca_ppm"="#FEE090",
                                 "Geographic"="black",
                                 "TON"="#BC2A22")




bact_MLE_soil_meta_names=c("past_7d_rain"="Seven day rain",
                                 "ph"="Soil pH",
                                 "k_ppm"="Soil K",
                                 "ca_ppm"="Soil Ca",
                                 "Geographic"="Spatial factor",
                                 "TON"="Soil TON")

ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',linewidth=2)+
  scale_color_manual(values = bact_MLE_soil_meta_color,name=NULL,labels=bact_MLE_soil_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")

#GDM MLE Fungal Roots####




GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data=sample_data(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root)[,c("sampleID_long","sampleID_seq","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data)
#111
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data)

colnames(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data)[colnames(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata=
  data.frame(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic=merge(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic)
#111  16



GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root,method = "bray")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF),GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic)
nrow(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic)

GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   47.849




GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","core_root_mass_subplot","k_ppm","ca_ppm",
                                                                       "TOC","TOC_TON_ratio","percent_soil_moisture_dry_weight",
                                                                       "past_7d_rain","ugN_NO3_g_dry_soil_na")]
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub1_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub1_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub1_gdmTabMod)
#47.849

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub1_TabForm, geo = T, predSelect = T, nPerm = 100)



#Final set of predictors returned: 
#Geographic
#ph
#core_root_mass_subplot
#k_ppm


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub=
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","k_ppm","core_root_mass_subplot")]
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   46.709

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm, geo = T, nPerm = 100)
#Model deviance                    270.827
#Percent deviance explained         46.703
#Model p-value                       0.000
#Fitted permutations               100.000

fung_R_meta_varSet <- vector("list", 2)

names(fung_R_meta_varSet)=c("plant","soil")
fung_R_meta_varSet$plant=c("core_root_mass_subplot")
fung_R_meta_varSet$soil=c("ph","k_ppm")
summary(fung_R_meta_varSet)
#Variance partitioning
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm,varSets = fung_R_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART)



GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_x_long=GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_x_long,
        GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_y_long)


VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")

fung_MLE_root_meta_color=c("core_root_mass_subplot"="#006837",
                                 "ph"="#7C001D",
                                 "k_ppm"="#F46D43",
                                 "Geographic"="black")




fung_MLE_root_meta_names=c("core_root_mass_subplot"="Root biomass",
                                 "ph"="Soil pH",
                                 "k_ppm"="Soil K",
                                 "Geographic"="Spatial factor")




ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = fung_MLE_root_meta_color,name=NULL,labels=fung_MLE_root_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")


#GDM MLE Fungal Soils####

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data=sample_data(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)[,c("sampleID_long","sampleID_seq","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data)
#112
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data)

colnames(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data)[colnames(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata=
  data.frame(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic=merge(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_MLE_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)
#112  16



GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil,method = "bray")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF),GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)
nrow(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   70.146



#Remove factors with less than 0.01
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","dry_matter_yield_mg_ha_mean",
                                                                       "percent_soil_moisture_dry_weight",
                                                                       "TON","total.shoot.dry.weight","ca_ppm",
                                                                       "ugN_NH4_g_dry_soil_na","k_ppm","TOC",
                                                                       "past_7d_rain")]
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub1_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub1_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub1_gdmTabMod)
#70.146

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub1_TabForm, geo = T, predSelect = T, nPerm = 100)

#Final set of predictors returned: 
#Geographic
#ph
#dry_matter_yield_mg_ha_mean







GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub=
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "dry_matter_yield_mg_ha_mean","ph")]
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   67.55

gdm.varImp_MOD(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm, geo = T, nPerm = 100)
#                           All predictors
#Model deviance                     125.08
#Percent deviance explained          67.55
#Model p-value                        0.00
#Fitted permutations                 97.00


fung_S_meta_varSet <- vector("list", 2)

names(fung_S_meta_varSet)=c("plant","soil")
fung_S_meta_varSet$plant=c("dry_matter_yield_mg_ha_mean")
fung_S_meta_varSet$soil=c("ph")
summary(fung_S_meta_varSet)
#Variance partitioning
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm,
                                                                                   varSet = fung_S_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART)



GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_x_long=GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_x_long,
        GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long)


VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")

fung_MLE_soil_meta_color=c("ph"="#7C001D",
                                 "Geographic"="black",
                                 "dry_matter_yield_mg_ha_mean"="#66BD63")




fung_MLE_soil_meta_names=c("ph"="Soil pH",
                                 "Geographic"="Spatial factor",
                                 "dry_matter_yield_mg_ha_mean"="Subplot yield")




ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = fung_MLE_soil_meta_color,name=NULL,labels=fung_MLE_soil_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")





#MLE: Graphing the Eulerr diagrams of partitioning of deviance explained####



GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART$community=rep("Fungi")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART$niche=rep("Soil")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART2=
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART[9:15,]
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART$community=rep("Bacteria")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART$niche=rep("Soil")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART2=
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART[9:15,]
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART$community=rep("Fungi")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART$niche=rep("Root")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART2=
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART[9:15,]
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART$community=rep("Bacteria")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART$niche=rep("Root")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART2=
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART[9:15,]


GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART=
  rbind(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_TabForm_PART2,
        GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_TabForm_PART2,
        GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_TabForm_PART2,
        GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_TabForm_PART2)

unique(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART$VARIABLE_SET)
length(unique(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART$VARIABLE_SET))

GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART$VARIABLE_SET_fix=
  GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART$VARIABLE_SET

VARIABLE_SET_order_f=c("geo alone","soil alone","MET alone","plant alone",
                       "geo intersect soil, exclude MET",
                       "MET intersect geo, exclude soil",
                       "plant intersect geo, exclude soil",
                       "geo intersect soil, exclude plant",
                       "MET intersect soil, exclude geo",
                       "plant intersect soil, exclude geo",
                       "soil intersect MET intersect geo",
                       "soil intersect plant intersect geo")

euler_VARIABLE_SET_names=c("geo alone"="Spatial",
                           "soil alone"="Soil",
                           "MET alone"="MET",
                           "plant alone"="Plant",
                           "geo intersect soil, exclude MET"="Spatial&Soil",
                           "MET intersect geo, exclude soil"="Spatial&MET",
                           "plant intersect geo, exclude soil"="Spatial&Plant",
                           "geo intersect soil, exclude plant"="Spatial&Soil",
                           "MET intersect soil, exclude geo"="Soil&MET",
                           "plant intersect soil, exclude geo"="Soil&Plant",
                           "soil intersect MET intersect geo"="Spatial&Soil&MET",
                           "soil intersect plant intersect geo"="Spatial&Soil&Plant")

GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART$VARIABLE_SET_euler=
  GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART$VARIABLE_SET

GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART <- GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART %>% 
  mutate(VARIABLE_SET_euler=mapvalues(VARIABLE_SET_euler,
                                      from = VARIABLE_SET_order_f,
                                      to = euler_VARIABLE_SET_names))

GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART|>
  group_by(community,niche)|>
  summarise(sum(DEVIANCE))
MET_SOIL_SPAT_alone_color=c("MET"="#313695","Soil"="#A50026","Spatial"="darkgrey")
PLANT_SOIL_SPAT_alone_color=c("Plant"="#006837","Soil"="#A50026","Spatial"="darkgrey")
VARIABLE_SET_euler_alone_color=c("Spatial"="darkgrey",
                                 "Soil"="#A50026",
                                 "MET"="#313695",
                                 "Plant"="#006837")

MLE_Root_bact_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Root")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Root")$VARIABLE_SET_euler)

MLE_Root_bact_euler_fit=euler(MLE_Root_bact_euler_list)



error_plot(MLE_Root_bact_euler_fit)




(MLE_Root_bact_euler_p=plot(MLE_Root_bact_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = MET_SOIL_SPAT_alone_color))







MLE_Soil_bact_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Soil")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Soil")$VARIABLE_SET_euler)

MLE_Soil_bact_euler_fit=euler(MLE_Soil_bact_euler_list)

error_plot(MLE_Soil_bact_euler_fit)

(MLE_Soil_bact_euler_p=plot(MLE_Soil_bact_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = MET_SOIL_SPAT_alone_color))

plot_grid(MLE_Root_bact_euler_p,MLE_Soil_bact_euler_p,ncol = 1)


MLE_Root_fung_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Fungi" &niche=="Root")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Fungi" &niche=="Root")$VARIABLE_SET_euler)

MLE_Root_fung_euler_fit=euler(MLE_Root_fung_euler_list)



error_plot(MLE_Root_fung_euler_fit)




(MLE_Root_fung_euler_p=plot(MLE_Root_fung_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = PLANT_SOIL_SPAT_alone_color))







MLE_Soil_fung_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Fungi" &niche=="Soil")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_MLE_G5_sub_TabForm_PART,community=="Fungi" &niche=="Soil")$VARIABLE_SET_euler)

MLE_Soil_fung_euler_fit=euler(MLE_Soil_fung_euler_list)

error_plot(MLE_Soil_fung_euler_fit)

(MLE_Soil_fung_euler_p=plot(MLE_Soil_fung_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = PLANT_SOIL_SPAT_alone_color))



(MLE_compart_fung_soil_legend=plot_ordination(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil,GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))


GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_NMDS_points=merge(sample_data(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root),GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_ord$points,
                                                             by="row.names")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_NMDS_points_sum=GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

(MLE_compart_bact_root_mean_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+ggtitle("Bacteria")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+ylab("NMDS2")+annotate("text", x = -0.55, y = -0.54, label = "Stress = 0.123", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5)),axis.title.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))




GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_NMDS_points=merge(sample_data(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root),GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_ord$points,
                                                             by="row.names")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_NMDS_points_sum=GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

(MLE_compart_fung_root_mean_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+ggtitle("Fungi")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+annotate("text", x = -0.495, y = -0.77, label = "Stress = 0.219", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5)), axis.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_NMDS_points=merge(sample_data(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil),GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_ord$points,
                                                             by="row.names")
GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_NMDS_points_sum=GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

(MLE_compart_bact_soil_mean_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                          y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 20)+
    geom_errorbar(aes(color=siteID),width=0.05)+geom_errorbarh(aes(color=siteID),height=0.05)+xlab("NMDS1")+ylab("NMDS2")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+annotate("text", x = -0.51, y = -0.6, label = "Stress = 0.120", size=6)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_NMDS_points=merge(sample_data(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil),GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_ord$points,
                                                             by="row.names")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_NMDS_points_sum=GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_NMDS_points%>%group_by(siteID,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))


(MLE_compart_fung_soil_mean_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
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
y_titles=plot_grid(get_title(root_title),get_title(soil_title),ggplot+geom_blank(),nrow = 3)

(MLE_mean_8_panel=plot_grid(y_titles,plot_grid(plot_grid(MLE_compart_bact_root_mean_p2,MLE_compart_fung_root_mean_p2,
                                 MLE_compart_bact_soil_mean_p2,MLE_compart_fung_soil_mean_p2,
                                 ncol = 2,
                                 label_size = 36,rel_widths = c(1.05,1,1.05,1),
                                 rel_heights = c(1,1,1.05,1.05),
                                 labels = c('a)', 'b)','c)', 'd)'),
                                 label_x = c(0,-0.035,0,-0.035),
                                 label_y = c(0.92,0.92,0.98,0.98)),
                                 plot_grid(as.ggplot(MLE_Root_bact_euler_p)+ggtitle("Root deviance\nexplained = 59.7%")+
                                             theme(plot.background=element_rect(fill = "white",colour = "white"),
                                                   plot.title = element_text(size = 20,hjust = 0.85)),
                                           as.ggplot(MLE_Soil_bact_euler_p)+ggtitle("Soil deviance\nexplained = 82.5%")+
                                             theme(plot.background=element_rect(fill = "white",colour = "white"),
                                                   plot.title = element_text(size = 20,hjust = 0.85)),
                                           as.ggplot(MLE_Root_fung_euler_p)+ggtitle("Root deviance\nexplained = 46.7%")+
                                             theme(plot.background=element_rect(fill = "white",colour = "white"),
                                                   plot.title = element_text(size = 20,hjust = 0.85)),
                                           as.ggplot(MLE_Soil_fung_euler_p)+ggtitle("Soil deviance\nexplained = 68.3%")+
                                             theme(plot.background=element_rect(fill = "white",colour = "white"),
                                                   plot.title = element_text(size = 20,hjust = 0.85)),ncol = 4,
                                           label_size = 36,labels = c('e)', 'f)','g)', 'h)'),
                                           label_x = c(0,-0.035,0,-0.035),
                                           label_y = c(0.92,0.92,0.92,0.92)),
                                 nrow = 2,rel_heights = c(1,.45)),
                                 ncol = 2,rel_widths = c(.15,4)))

plot_grid(MLE_mean_8_panel,get_legend(MLE_compart_fung_soil_legend),ncol = 2,rel_widths = c(6,1))

#NOT INCLUDED IN REPOSITORY
ggsave(plot_grid(MLE_mean_8_panel,get_legend(MLE_compart_fung_soil_legend),ncol = 2,rel_widths = c(4,1)), 
       filename = "NMDS_Euler_All_Sites_mean_p.png",path = here::here("Manuscript","MLE_comm_figs"),width = 22,height = 22)


ggsave(plot_grid(MLE_mean_8_panel,get_legend(MLE_compart_fung_soil_legend),ncol = 2,rel_widths = c(4,1)), 
       filename = "NMDS_Euler_All_Sites_mean_p.svg",path = here::here("Manuscript","MLE_comm_figs"),width = 22,height = 22)
#NOT INCLUDED IN REPOSITORY

#MLE: Graphing the combined GDM graphs####

All_sample_GDM_color=c("Geographic"="black",
                       # "MET measurement"="#313695"
                       "past_7d_rain"="#4575B4",
                       "percent_soil_moisture_dry_weight"="#313695",
                       "soil_temp_1_avg_mean"="#74ADD1",
                       #"Soil nutrients"="#A50026"
                       "ph"="#7C001D",
                       "k_ppm"="#F46D43",
                       "ca_ppm"="#FEE090",
                       "TON"="#BC2A22",
                       "p_ppm"="#D73027",
                       #"Plant traits"="#006837"
                       "core_root_mass_subplot"="#006837",
                       "dry_matter_yield_mg_ha_mean"="#66BD63",
                       "collectionDate_N"="darkgrey"
)



All_sample_GDM_names=c("Geographic"="Spatial factor",
                       "past_7d_rain"="Seven day\nrain accumulation",
                       "percent_soil_moisture_dry_weight"="Gravimetric\nsoil moisture",
                       "soil_temp_1_avg_mean"="Soil temp (24h avg)",
                       "ph"="Soil pH",
                       "k_ppm"="Soil K",
                       "ca_ppm"="Soil Ca",
                       "TON"="Soil TON",
                       "p_ppm"="Soil P",
                       "core_root_mass_subplot"="Root biomass",
                       "dry_matter_yield_mg_ha_mean"="Subplot yield",
                       "collectionDate_N"="Temporal factor"
)




All_MLE_sample_GDM_color=All_sample_GDM_color[names(All_sample_GDM_color)%in%
                                                unique(c(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long$metadata))]

All_MLE_sample_GDM_names=All_sample_GDM_names[names(All_sample_GDM_names)%in%
                                                unique(c(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long$metadata))]



(MLE_bact_root_sub_GDM=ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    theme_cowplot(font_size = 24)+labs(x=NULL,y="Root\nPartial ecological distance\n(Bray-Curtis)")+
    ggtitle("Bacteria")+annotation_custom(grobTree(textGrob(paste("p < 0.001\n59.7%"), x=0.1,  y=0.90, hjust=0,
                                                            gp=gpar(col="black", fontsize=26))))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 36)))

(MLE_bact_soil_sub_GDM=ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n82.5%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Soil\nPartial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))

(MLE_fung_root_sub_GDM=ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n46.7%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x=NULL,y=NULL)+ggtitle("Fungi")+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 36)))


(MLE_fung_soil_sub_GDM=ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n67.6%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.position = "none"))


All_MLE_sample_GDM_order=c("Geographic","ph","TON","k_ppm","ca_ppm",
                           "past_7d_rain",
                           "core_root_mass_subplot","dry_matter_yield_mg_ha_mean")

(MLE_fung_soil_sub_GDM_legend=ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_sub_gdm_long,
                                      aes(x=x_predictors,y=y_values,color=metadata))+
    geom_smooth(linetype="solid",se=FALSE,method = 'loess',linewidth=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names,limits=All_MLE_sample_GDM_order)+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.key.size = unit(2, 'cm'),legend.position = "bottom"))

(four_panel_GDM_MLE=plot_grid(plot_grid(MLE_bact_root_sub_GDM,MLE_fung_root_sub_GDM,
                                        MLE_bact_soil_sub_GDM,MLE_fung_soil_sub_GDM,labels = c("a)","b)","c)","d)"),
                                        ncol = 2, rel_widths = c(1,1), label_size = 32,
                                        label_x = c(0,0.08,0,0.08),label_y = c(1,1,1.05,1.05),
                                        align = "v"),plot_grid(ggplot+geom_blank(),get_legend(MLE_fung_soil_sub_GDM_legend),
                                                               rel_widths = c(0.4,1)), rel_heights = c(1,0.2),ncol = 1))

#NOT INCLUDED IN REPOSITORY
ggsave(four_panel_GDM_MLE,
       filename = "GDM_loess_metadata_MLE_four_panel.png",path = here::here("Manuscript","GDM_metadata_figs"),width = 20,height =14)

ggsave(four_panel_GDM_MLE,
       filename = "GDM_loess_metadata_MLE_four_panel.svg",path = here::here("Manuscript","GDM_metadata_figs"),width = 20,height =14)
#NOT INCLUDED IN REPOSITORY



#####MLE Beta-disp####



GLBRC018_OTU_bact_MMPRNT_MLE_G5_betamod=betadisper(GLBRC018_OTU_bact_MMPRNT_MLE_G5_dis, with(GLBRC018_OTU_bact_MMPRNT_MLE_G5_map,interaction(Root_soil,siteID,FertStatus)))
GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis=as.data.frame(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betamod$distances)
GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis=merge(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis,GLBRC018_OTU_bact_MMPRNT_MLE_G5_map,by="row.names")
colnames(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis)[colnames(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis)==
                                                          "GLBRC018_OTU_bact_MMPRNT_MLE_G5_betamod$distances"]="betadisp"



GLBRC018_OTU_fung_MMPRNT_MLE_G5_betamod=betadisper(GLBRC018_OTU_fung_MMPRNT_MLE_G5_dis, with(GLBRC018_OTU_fung_MMPRNT_MLE_G5_map,interaction(Root_soil,siteID,FertStatus)))
GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis=as.data.frame(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betamod$distances)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis=merge(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis,GLBRC018_OTU_fung_MMPRNT_MLE_G5_map,by="row.names")
colnames(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis)[colnames(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis)==
                                                          "GLBRC018_OTU_fung_MMPRNT_MLE_G5_betamod$distances"]="betadisp"

head(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis)
unique(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis$siteID)



#Stats 

#Roots Bacteria
GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis_root=subset(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis,Root_soil=="Root")

root_betadisp_MLE_bact_mod=lmer(log(betadisp)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis_root)
plot(root_betadisp_MLE_bact_mod)
hist(resid(root_betadisp_MLE_bact_mod))
qqPlot(resid(root_betadisp_MLE_bact_mod))
shapiro.test(resid(root_betadisp_MLE_bact_mod))
#W = 0.98237, p-value = 0.143
simulateResiduals(fittedModel = root_betadisp_MLE_bact_mod, plot = T)

anova(root_betadisp_MLE_bact_mod)
#siteID            0.37223 0.093057     4 101.18  6.3060 0.0001428 ***
#FertStatus        0.01573 0.015731     1 100.31  1.0660 0.3043382    
#siteID:FertStatus 0.06075 0.015187     4 100.31  1.0292 0.3960468 


emmeans(root_betadisp_MLE_bact_mod,pairwise~siteID)
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
GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis_soil=subset(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis,Root_soil=="Soil")

soil_betadisp_MLE_bact_mod=lmer((betadisp)^-3~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis_soil)
plot(soil_betadisp_MLE_bact_mod)
hist(resid(soil_betadisp_MLE_bact_mod))
qqPlot(resid(soil_betadisp_MLE_bact_mod))
shapiro.test(resid(soil_betadisp_MLE_bact_mod))
#W = 0.97723, p-value = 0.04884
simulateResiduals(fittedModel = soil_betadisp_MLE_bact_mod, plot = T)

anova(soil_betadisp_MLE_bact_mod)
#siteID            18877.5  4719.4     4 101.90 28.7874 5.204e-16 ***
#FertStatus          168.8   168.8     1 101.03  1.0295    0.3127    
#siteID:FertStatus   167.3    41.8     4 101.03  0.2552    0.9059  



emmeans(soil_betadisp_MLE_bact_mod,pairwise~siteID)
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
GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis_root=subset(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis,Root_soil=="Root")

root_betadisp_MLE_fung_mod=lmer(sqrt(betadisp)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis_root)
plot(root_betadisp_MLE_fung_mod)
hist(resid(root_betadisp_MLE_fung_mod))
qqPlot(resid(root_betadisp_MLE_fung_mod))
shapiro.test(resid(root_betadisp_MLE_fung_mod))
#W = 0.98262, p-value = 0.1585
simulateResiduals(fittedModel = root_betadisp_MLE_fung_mod, plot = T)

anova(root_betadisp_MLE_fung_mod)
#siteID            0.078001 0.0195002     4   101  6.5540 9.909e-05 ***
#FertStatus        0.002818 0.0028182     1   101  0.9472    0.3328    
#siteID:FertStatus 0.002614 0.0006535     4   101  0.2197    0.9269  



emmeans(root_betadisp_MLE_fung_mod,pairwise~siteID)
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
GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis_soil=subset(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis,Root_soil=="Soil")

soil_betadisp_MLE_fung_mod=lmer(betadisp~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis_soil)
plot(soil_betadisp_MLE_fung_mod)
hist(resid(soil_betadisp_MLE_fung_mod))
qqPlot(resid(soil_betadisp_MLE_fung_mod))
shapiro.test(resid(soil_betadisp_MLE_fung_mod))
#W = 0.97726, p-value = 0.05284
simulateResiduals(fittedModel = soil_betadisp_MLE_fung_mod, plot = T)

anova(soil_betadisp_MLE_fung_mod)
#siteID            0.270100 0.067525     4 100.062 17.5100 6.387e-11 ***
#FertStatus        0.005119 0.005119     1  99.196  1.3274    0.2520    
#siteID:FertStatus 0.019024 0.004756     4  99.164  1.2333    0.3017    



emmeans(soil_betadisp_MLE_fung_mod,pairwise~siteID)
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


#MLE Beta-disp figures####
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
bact_betadisp_max=GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(betadisp))

bact_betadisp_max_disp_letters=merge(root_soil_bact_disp_letters,bact_betadisp_max,by=c("siteID","Root_soil"))




root_fung_disp_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("a","b","bc","ac","c"),
                                  Root_soil=rep("Root"))


soil_fung_disp_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("AB","A","C","B","A"),
                                  Root_soil=rep("Soil"))
root_soil_fung_disp_letters=rbind(root_fung_disp_letters,soil_fung_disp_letters)

fung_betadisp_max=GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(betadisp))

fung_betadisp_max_disp_letters=merge(root_soil_fung_disp_letters,fung_betadisp_max,by=c("siteID","Root_soil"))







(bact_betadisp_p3=ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),
                     fill=FertStatus,color=FertStatus))+
    geom_text(data = bact_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_fill_manual(values = c("black","grey"), name=NULL)+
    scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    scale_y_continuous(name = "Betadispersion",limits = c(0.22,0.57))+theme(axis.text.x = element_blank()))





(bact_betadisp_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_MLE_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = bact_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    scale_y_continuous(name = "Bacterial\nBetadispersion",limits = c(0.22,0.57))+
    theme(axis.text.x = element_blank(),legend.position = "none",
          strip.background = element_rect(fill="white",color = "black")))

(fung_betadisp_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_MLE_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = fung_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Fungal\nBetadispersion",limits = c(0.2,0.77))+ 
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none"))


(betadisp_2panel=plot_grid(bact_betadisp_p2,fung_betadisp_p2,ncol = 1,labels = c('a)', 'b)'), label_size = 30))

#NOT INCLUDED IN REPOSITORY
ggsave(plot_grid(betadisp_2panel,get_legend(bact_betadisp_p3),ncol = 2,rel_widths = c(4,1)),
       filename = "Betadisp_boxplot_All_Sites_p.png",path = here::here("Manuscript","MLE_comm_figs"),width = 20,height =10)

ggsave(plot_grid(betadisp_2panel,get_legend(bact_betadisp_p3),ncol = 2,rel_widths = c(4,1)),
       filename = "Betadisp_boxplot_All_Sites_p.svg",path = here::here("Manuscript","MLE_comm_figs"),width = 20,height =10)

#NOT INCLUDED IN REPOSITORY


