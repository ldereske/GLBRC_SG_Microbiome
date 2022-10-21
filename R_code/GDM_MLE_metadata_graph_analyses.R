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
library(gdm)
library(GGally)
library(missForest)
library(stringr)
library(otuSummary)
library(grid)
library(eulerr)

#Load code for calculating importance and p-values for GDMs with 3 or less factors
source(here::here("R_functions","gdm.varImp_MOD.R"))

#Metadata for community predictions


MMPRNT_2018.metadata_mapp_coverage= read.csv(here::here("Publish_data","MMPRNT_2018.metadata_mapp_coverage.csv"),
                                             header = T)

dim(MMPRNT_2018.metadata_mapp_coverage)
#870  73

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


#####All Site Abiotic####



MMPRNT_2018.metadata_mapp_coverage_all_sites=subset(MMPRNT_2018.metadata_mapp_coverage,siteID!="LUX"|collectionDate=="7/30/2018")
nrow(MMPRNT_2018.metadata_mapp_coverage_all_sites)
#486


MMPRNT_2018.metadata_mapp_coverage_all_sites=subset(MMPRNT_2018.metadata_mapp_coverage_all_sites,siteID!="LC"|collectionDate=="7/10/2018")
nrow(MMPRNT_2018.metadata_mapp_coverage_all_sites)
#294

#matching to the actual sequence samples

MMPRNT_2018.metadata_mapp_coverage_all_sites_G5=subset(MMPRNT_2018.metadata_mapp_coverage_all_sites,plotType=="G5")
nrow(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5)
#114
colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5))


MMPRNT_2018.metadata_mapp_all_sites_G5=MMPRNT_2018.metadata_mapp_coverage_all_sites_G5
MMPRNT_2018.metadata_mapp_all_sites_G5[,c("event_age","pH_raw_MMPRNT","pH_MMPRNT_subplot","rain_mm_tot_sum")]=NULL

colnames(MMPRNT_2018.metadata_mapp_all_sites_G5)
colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5))

cor(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5[,c("percent_soil_moisture_dry_weight","vwc_avg_mean","dry_matter_yield_mg_ha_mean",
                                                       "total.shoot.dry.weight","plant.height","specific.leaf.area","core_root_mass_subplot",
                                                       "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na","soil_temp_1_avg_mean","airtc_avg_mean","past_7d_rain",
                                                       "ph","p_ppm","k_ppm","ca_ppm","mg_ppm","ugN_fixed_gdrysoil_day","TOC","TON","TOC_TON_ratio"
)],use = "na.or.complete")

write.csv(cor(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5[,c("percent_soil_moisture_dry_weight","dry_matter_yield_mg_ha_mean",
                                                                 "total.shoot.dry.weight","plant.height","specific.leaf.area","core_root_mass_subplot",
                                                                 "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na","soil_temp_1_avg_mean","airtc_avg_mean","past_7d_rain",
                                                                 "ph","p_ppm","k_ppm","ca_ppm","mg_ppm","TOC","TON","TOC_TON_ratio"
)],use = "na.or.complete",method = "spearman"),file = here::here("R_files","ignorable_files","corr_mat_MLE_Abiotic_Biotic.csv"))






ggpairs(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5, columns = c("percent_soil_moisture_dry_weight","dry_matter_yield_mg_ha_mean","total.shoot.dry.weight","core_root_mass_subplot","ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na","past_7d_rain","ph","k_ppm","ca_ppm","TOC","TON","TOC_TON_ratio"), 
        upper = list(continuous = wrap("cor",size = 3)),
        lower = list(continuous = wrap("smooth",
                                       alpha = 0.3,
                                       size = 0.1)))




#Metadata Extraction 





MMPRNT_2018.metadata_mapp_all_sites_G5_lim=
  data.frame(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("SampleID","siteID","percent_soil_moisture_dry_weight","dry_matter_yield_mg_ha_mean",
                                                       "total.shoot.dry.weight","core_root_mass_subplot","ugN_NH4_g_dry_soil_na",
                                                       "ugN_NO3_g_dry_soil_na","past_7d_rain","ph","k_ppm","ca_ppm","TOC","TON","TOC_TON_ratio")])

colSums(is.na(MMPRNT_2018.metadata_mapp_all_sites_G5_lim))

MMPRNT_2018.metadata_mapp_all_sites_G5_lim|>group_by(siteID)|>summarise(biomass=sum(is.na(total.shoot.dry.weight)),
                                                                        Soil_N=sum(is.na(ugN_NO3_g_dry_soil_na)),
                                                                        root=sum(is.na(core_root_mass_subplot)),
                                                                        TOC_N=sum(is.na(TOC_TON_ratio)))


#Impute the missing data
summary(MMPRNT_2018.metadata_mapp_all_sites_G5_lim)
row.names(MMPRNT_2018.metadata_mapp_all_sites_G5_lim)=MMPRNT_2018.metadata_mapp_all_sites_G5_lim$SampleID

MMPRNT_2018.metadata_mapp_all_sites_G5_lim$siteID=as.factor(MMPRNT_2018.metadata_mapp_all_sites_G5_lim$siteID)

set.seed(2022)
MMPRNT_2018.metadata_mapp_all_sites_G5_lim.imp <- missForest(MMPRNT_2018.metadata_mapp_all_sites_G5_lim[,-1])


MMPRNT_2018.metadata_mapp_all_sites_G5_lim.imp$ximp

#check imputation error
MMPRNT_2018.metadata_mapp_all_sites_G5_lim.imp$OOBerror
#     NRMSE        PFC 
#0.05325305 0.00000000  

MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp=MMPRNT_2018.metadata_mapp_all_sites_G5_lim.imp$ximp
head(MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp)

MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp$sampleID_long=
  paste("MMPRNT",str_sub(row.names(MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp),start = 7,end = 11),sep = "-")
#GDM All Sites Bacterial Roots####



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data=sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)[,c("sampleID_long","SampleID","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data)
#113
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data)

colnames(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data)[colnames(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata=
  data.frame(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic=merge(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)
#113  16



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root,method = "bray")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF),GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)
nrow(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   61.339
#gdm.varImp(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_TabForm, geo = T, predSelect = T,nPerm = 100)



#Remove factors with less than 0.01
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","past_7d_rain","k_ppm","ca_ppm",
                                                                       "ugN_NO3_g_dry_soil_na","core_root_mass_subplot",
                                                                       "ugN_NH4_g_dry_soil_na","TOC_TON_ratio",
                                                                       "dry_matter_yield_mg_ha_mean","TON")]
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub1_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub1_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub1_gdmTabMod)
#plot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   61.339

gdm.varImp(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub1_TabForm, geo = T, predSelect = T,nPerm = 100)

#Final set of predictors returned: 
#Geographic
#ph
#past_7d_rain
#k_ppm


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "past_7d_rain","ph","k_ppm")]
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   59.682

gdm.varImp(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm, geo = T, nPerm = 100)
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
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm,varSets = bact_R_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART)

sum(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART[9:15,2])


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_x_long=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_x_long,
        GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_y_long)






VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")





bact_All_sites_root_meta_color=c("past_7d_rain"="#4575B4",
                                 "ph"="#7C001D",
                                 "k_ppm"="#F46D43",
                                 "Geographic"="black")



bact_All_sites_root_meta_names=c("past_7d_rain"="Seven day rain",
                                 "ph"="Soil pH",
                                 "k_ppm"="Soil K",
                                 "Geographic"="Spatial factor")


dev.off()
ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = bact_All_sites_root_meta_color,name=NULL,labels=bact_All_sites_root_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")





GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_x_long,
        GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long)

row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$site

head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)
All_sites_bact_root_dist=matrixConvert(dist(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub[,c("UTM_Lon_Cord", "UTM_Lat_Cord")]), 
                                       colname = c("sample1", "sample2", "distance"))#total distance 
dim(All_sites_bact_root_dist)
#6328    3
All_sites_bact_root_dist_site=merge(All_sites_bact_root_dist,
                                    sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)[,c("siteID","plotRep")],
                                    by.x = "sample1",by.y = "row.names")
dim(All_sites_bact_root_dist_site)
#6328    5
colnames(All_sites_bact_root_dist_site)[4:5]=c("s1_siteID","s1_plotRep")

All_sites_bact_root_dist_site=merge(All_sites_bact_root_dist_site,
                                    sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)[,c("siteID","plotRep")],
                                    by.x = "sample2",by.y = "row.names")
dim(All_sites_bact_root_dist_site)
#6328    7
colnames(All_sites_bact_root_dist_site)[6:7]=c("s2_siteID","s2_plotRep")


All_sites_bact_root_dist_site$s1_s2_siteID=with(All_sites_bact_root_dist_site,
                                                interaction(s1_siteID,s2_siteID))
unique(All_sites_bact_root_dist_site$s1_s2_siteID)
All_sites_bact_root_dist_site$site_comp=with(All_sites_bact_root_dist_site,
                                             ifelse(s1_siteID==s2_siteID,"within",
                                                    as.character(All_sites_bact_root_dist_site$s1_s2_siteID)))
unique(All_sites_bact_root_dist_site$site_comp)

All_sites_sq=All_sites_bact_root_dist_site%>%group_by(site_comp)%>%
  summarise(max_dis=max(distance),min_dis=min(distance))

head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)

MLE_pH_points= data.frame(pH_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$ph))
summary(MLE_pH_points)

(MLE_bact_root_pH_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="ph"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_root_meta_color,name=NULL,labels=bact_All_sites_root_meta_names)+ylim(c(-0.005,0.43))+
    geom_segment(data=MLE_pH_points,aes(x=pH_p,xend = pH_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil pH",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))



(MLE_bact_root_spat_gdm=ggplot(data=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="Geographic"))+
    geom_smooth(aes(x=x_predictors,y=y_values,color=metadata),linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_root_meta_color,name=NULL,labels=bact_All_sites_root_meta_names)+
    geom_rect(data = All_sites_sq,aes(xmin=min_dis,xmax = max_dis,group=site_comp),
              size=2,color="black",fill="black",ymin=-0.005,ymax=-Inf)+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)",limits = c(-0.005,0.43))+
    labs(x="Distance (m)")+
    theme(legend.position = "none"))

MLE_past_7d_rain_points= data.frame(past_7d_rain_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$past_7d_rain))
summary(MLE_past_7d_rain_points)

(MLE_bact_root_past_7d_rain_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="past_7d_rain"),
                                       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_root_meta_color,name=NULL,labels=bact_All_sites_root_meta_names)+ylim(c(-0.005,0.43))+
    geom_segment(data=MLE_past_7d_rain_points,aes(x=past_7d_rain_p,
                                                  xend = past_7d_rain_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Seven day rain accumulation (mm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))

MLE_k_ppm_points= data.frame(k_ppm_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$k_ppm))
summary(MLE_k_ppm_points)

(MLE_bact_root_k_ppm_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="k_ppm"),
                                aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_root_meta_color,name=NULL,labels=bact_All_sites_root_meta_names)+ylim(c(-0.005,0.43))+
    geom_segment(data=MLE_k_ppm_points,aes(x=k_ppm_p,
                                       xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))

(LEGEND_MLE_bact_root_k_ppm_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="k_ppm"),
                                       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_root_meta_color,name=NULL,labels=bact_All_sites_root_meta_names)+ylim(c(-0.005,0.43))+
    geom_segment(data=MLE_k_ppm_points,aes(x=k_ppm_p,
                                       xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))

plot_grid(plot_grid(MLE_bact_root_pH_gdm,MLE_bact_root_past_7d_rain_gdm,
                    MLE_bact_root_spat_gdm,MLE_bact_root_k_ppm_gdm,rel_widths = c(1.3,1),
                    labels = c("a)","b)","c)","d)"),label_size = 36,label_x = c(0,-0.1)),
          get_legend(LEGEND_MLE_bact_root_k_ppm_gdm),rel_widths = c(1,0.3))



ggplot()+ggtitle("Root Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Root Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(plot_grid(MLE_bact_root_pH_gdm,MLE_bact_root_past_7d_rain_gdm,
                              MLE_bact_root_spat_gdm,MLE_bact_root_k_ppm_gdm,rel_widths = c(1.3,1),
                              labels = c("a)","b)","c)","d)"),label_size = 36,label_x = c(0,-0.1)),
                    get_legend(LEGEND_MLE_bact_root_k_ppm_gdm),rel_widths = c(1,0.3)),nrow = 2,
          rel_heights = c(0.075,1))







#GDM All Sites Bacterial Soils####




GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data=sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)[,c("sampleID_long","SampleID","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data)
#114
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data)

colnames(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data)[colnames(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata=
  data.frame(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic=merge(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)
#114  16



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil,method = "bray")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF),GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)
nrow(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   82.883



#Remove factors with less than 0.01
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord","ph","ca_ppm",
                                                                       "TON","past_7d_rain","percent_soil_moisture_dry_weight",
                                                                       "k_ppm","ugN_NO3_g_dry_soil_na","dry_matter_yield_mg_ha_mean",
                                                                       "TOC_TON_ratio")]
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub1_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub1_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub1_gdmTabMod)
#82.781
gdm.varImp(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub1_TabForm, geo = T, predSelect = T, nPerm = 100)
#Final set of predictors returned: 
#Geographic
#ph
#ca_ppm
#TON
#past_7d_rain
#k_ppm






GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "past_7d_rain","ph","k_ppm","ca_ppm","TON")]
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   82.465

gdm.varImp(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm, geo = T, nPerm = 100)
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
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm,
                                                                                   varSet = bact_S_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART)

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_x_long=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_x_long,
        GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_long)





VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")


bact_All_sites_soil_meta_color=c("past_7d_rain"="#4575B4",
                                 "ph"="#7C001D",
                                 "k_ppm"="#F46D43",
                                 "ca_ppm"="#FEE090",
                                 "Geographic"="black",
                                 "TON"="#BC2A22")




bact_All_sites_soil_meta_names=c("past_7d_rain"="Seven day rain",
                                 "ph"="Soil pH",
                                 "k_ppm"="Soil K",
                                 "ca_ppm"="Soil Ca",
                                 "Geographic"="Spatial factor",
                                 "TON"="Soil TON")



#dev.off()
ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")



GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_x_long,
        GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long)

row.names(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$site





head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)

MLE_pH_points= data.frame(pH_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$ph))
summary(MLE_pH_points)

(MLE_bact_soil_pH_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="ph"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+ylim(c(-0.005,0.60))+
    geom_segment(data=MLE_pH_points,aes(x=pH_p,xend = pH_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil pH",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))



(MLE_bact_soil_spat_gdm=ggplot(data=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="Geographic"))+
    geom_smooth(aes(x=x_predictors,y=y_values,color=metadata),linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+
    geom_rect(data = All_sites_sq,aes(xmin=min_dis,xmax = max_dis,group=site_comp),
              size=2,color="black",fill="black",ymin=-0.005,ymax=-Inf)+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)",limits = c(-0.005,0.60))+
    labs(x="Distance (m)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))

MLE_past_7d_rain_points= data.frame(past_7d_rain_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$past_7d_rain))
summary(MLE_past_7d_rain_points)

(MLE_bact_soil_past_7d_rain_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="past_7d_rain"),
                                       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+ylim(c(-0.005,0.60))+
    geom_segment(data=MLE_past_7d_rain_points,aes(x=past_7d_rain_p,
                                                  xend = past_7d_rain_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Seven day rain accumulation (mm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))



MLE_TON_points= data.frame(TON_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$TON))
summary(MLE_TON_points)

(MLE_bact_soil_TON_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="TON"),
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+ylim(c(-0.005,0.60))+
    geom_segment(data=MLE_TON_points,aes(x=TON_p,
                                         xend = TON_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil TON",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))


MLE_ca_ppm_points= data.frame(ca_ppm_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$ca_ppm))
summary(MLE_ca_ppm_points)

(MLE_bact_soil_ca_ppm_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="ca_ppm"),
                                 aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+ylim(c(-0.005,0.60))+
    geom_segment(data=MLE_ca_ppm_points,aes(x=ca_ppm_p,
                                            xend = ca_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil Ca (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))


MLE_k_ppm_points= data.frame(k_ppm_p=unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$k_ppm))
summary(MLE_k_ppm_points)

(MLE_bact_soil_k_ppm_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="k_ppm"),
                                aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+ylim(c(-0.005,0.60))+
    geom_segment(data=MLE_k_ppm_points,aes(x=k_ppm_p,
                                           xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))

(LEGEND_MLE_bact_soil_k_ppm_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="k_ppm"),
                                       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_All_sites_soil_meta_color,name=NULL,labels=bact_All_sites_soil_meta_names)+ylim(c(-0.005,0.60))+
    geom_segment(data=MLE_k_ppm_points,aes(x=k_ppm_p,
                                           xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))

plot_grid(plot_grid(MLE_bact_soil_pH_gdm,MLE_bact_soil_ca_ppm_gdm,
                    MLE_bact_soil_spat_gdm,MLE_bact_soil_TON_gdm,
                    MLE_bact_soil_past_7d_rain_gdm,MLE_bact_soil_k_ppm_gdm,rel_widths = c(1.3,1,1),
                    labels = c("a)","b)","c)","d)","e)","f)"),label_size = 36,label_x = c(0,-0.1,-0.1)),
          get_legend(LEGEND_MLE_bact_soil_k_ppm_gdm),rel_widths = c(1,0.3))



ggplot()+ggtitle("Soil Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Soil Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(plot_grid(MLE_bact_soil_pH_gdm,MLE_bact_soil_ca_ppm_gdm,
                              MLE_bact_soil_spat_gdm,MLE_bact_soil_TON_gdm,
                              MLE_bact_soil_past_7d_rain_gdm,MLE_bact_soil_k_ppm_gdm,rel_widths = c(1.2,1,1),
                              labels = c("a)","b)","c)","d)","e)","f)"),label_size = 36,label_x = c(0,-0.1,-0.1)),
                    get_legend(LEGEND_MLE_bact_soil_k_ppm_gdm),rel_widths = c(1,0.3)),nrow = 2,
          rel_heights = c(0.075,1))





#GDM All Sites Fungal Roots####




GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data=sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)[,c("sampleID_long","sampleID_seq","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data)
#111
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data)

colnames(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data)[colnames(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata=
  data.frame(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic=merge(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)
#111  16



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root,method = "bray")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF),GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)
nrow(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   47.849




GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","core_root_mass_subplot","k_ppm","ca_ppm",
                                                                       "TOC","TOC_TON_ratio","percent_soil_moisture_dry_weight",
                                                                       "past_7d_rain","ugN_NO3_g_dry_soil_na")]
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub1_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub1_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub1_gdmTabMod)
#47.849

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub1_TabForm, geo = T, predSelect = T, nPerm = 100)



#Final set of predictors returned: 
#Geographic
#ph
#core_root_mass_subplot
#k_ppm


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","k_ppm","core_root_mass_subplot")]
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   46.709

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm, geo = T, nPerm = 100)
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
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm,varSets = fung_R_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART)

sum(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART[9:15,2])



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_x_long=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_x_long,
        GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_long)





VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")

fung_All_sites_root_meta_color=c("core_root_mass_subplot"="#006837",
                                 "ph"="#7C001D",
                                 "k_ppm"="#F46D43",
                                 "Geographic"="black")




fung_All_sites_root_meta_names=c("core_root_mass_subplot"="Root biomass",
                                 "ph"="Soil pH",
                                 "k_ppm"="Soil K",
                                 "Geographic"="Spatial factor")




ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = fung_All_sites_root_meta_color,name=NULL,labels=fung_All_sites_root_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_x_long,
        GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long)



head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub)

MLE_pH_points= data.frame(pH_p=unique(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$ph))
summary(MLE_pH_points)

(MLE_fung_root_pH_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="ph"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_root_meta_color,name=NULL,labels=fung_All_sites_root_meta_names)+ylim(c(-0.005,1.1))+
    geom_segment(data=MLE_pH_points,aes(x=pH_p,xend = pH_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil pH",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))



(MLE_fung_root_spat_gdm=ggplot(data=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="Geographic"))+
    geom_smooth(aes(x=x_predictors,y=y_values,color=metadata),linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_root_meta_color,name=NULL,labels=fung_All_sites_root_meta_names)+
    geom_rect(data = All_sites_sq,aes(xmin=min_dis,xmax = max_dis,group=site_comp),
              size=2,color="black",fill="black",ymin=-0.005,ymax=-Inf)+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)",limits = c(-0.005,1.1))+
    labs(x="Distance (m)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))

MLE_core_root_mass_subplot_points= data.frame(core_root_mass_subplot_p=unique(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$core_root_mass_subplot))
summary(MLE_core_root_mass_subplot_points)

(MLE_fung_root_core_root_mass_subplot_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="core_root_mass_subplot"),
                                                 aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_root_meta_color,name=NULL,labels=fung_All_sites_root_meta_names)+ylim(c(-0.005,1.1))+
    geom_segment(data=MLE_core_root_mass_subplot_points,aes(x=core_root_mass_subplot_p,
                                                            xend = core_root_mass_subplot_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Root biomass (g)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))

MLE_k_ppm_points= data.frame(k_ppm_p=unique(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_Space_metadata_abiotic_sub$k_ppm))
summary(MLE_k_ppm_points)

(MLE_fung_root_k_ppm_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="k_ppm"),
                                aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_root_meta_color,name=NULL,labels=fung_All_sites_root_meta_names)+ylim(c(-0.005,1.1))+
    geom_segment(data=MLE_k_ppm_points,aes(x=k_ppm_p,
                                           xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))

(LEGEND_MLE_fung_root_k_ppm_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_NS_long,metadata=="k_ppm"),
                                       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_root_meta_color,name=NULL,labels=fung_All_sites_root_meta_names)+ylim(c(-0.005,1.1))+
    geom_segment(data=MLE_k_ppm_points,aes(x=k_ppm_p,
                                           xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))

plot_grid(plot_grid(MLE_fung_root_pH_gdm,MLE_fung_root_spat_gdm,
                    MLE_fung_root_core_root_mass_subplot_gdm,MLE_fung_root_k_ppm_gdm,rel_widths = c(1.2,1),
                    labels = c("a)","b)","c)","d)"),label_size = 36,label_x = c(0,-0.1)),
          get_legend(LEGEND_MLE_fung_root_k_ppm_gdm),rel_widths = c(1,0.3))



ggplot()+ggtitle("Root Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Root Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(plot_grid(MLE_fung_root_pH_gdm,MLE_fung_root_spat_gdm,
                              MLE_fung_root_core_root_mass_subplot_gdm,MLE_fung_root_k_ppm_gdm,rel_widths = c(1.2,1),
                              labels = c("a)","b)","c)","d)"),label_size = 36,label_x = c(0,-0.1)),
                    get_legend(LEGEND_MLE_fung_root_k_ppm_gdm),rel_widths = c(1,0.3)),nrow = 2,
          rel_heights = c(0.075,1))



#GDM All Sites Fungal Soils####

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data=sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)[,c("sampleID_long","sampleID_seq","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data)
#112
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data)

colnames(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data)[colnames(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata=
  data.frame(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic=merge(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata,
                                                                        MMPRNT_2018.metadata_mapp_all_sites_G5_lim_imp, by = "sampleID_long")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic[,c("sampleID_long","siteID")]=NULL
dim(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)
#112  16



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil,method = "bray")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF),GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)
nrow(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                    XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                    predData=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   70.146



#Remove factors with less than 0.01
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "ph","dry_matter_yield_mg_ha_mean",
                                                                       "percent_soil_moisture_dry_weight",
                                                                       "TON","total.shoot.dry.weight","ca_ppm",
                                                                       "ugN_NH4_g_dry_soil_na","k_ppm","TOC",
                                                                       "past_7d_rain")]
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub1_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                         XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                         predData=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub1_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub1_gdmTabMod)
#70.146

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub1_TabForm, geo = T, predSelect = T, nPerm = 100)

#Final set of predictors returned: 
#Geographic
#ph
#dry_matter_yield_mg_ha_mean







GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                       "dry_matter_yield_mg_ha_mean","ph")]
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_bray_DF_site, 3, siteColumn="site", 
                                                                        XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                        predData=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   67.55

gdm.varImp_MOD(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm, geo = T, nPerm = 100)
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
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART=gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm,
                                                                                   varSet = fung_S_meta_varSet,partSpace= TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART)

sum(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART[9:15,2])


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$x,
                                                                          method = "range"))

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_x_long=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_x_long,
        GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long)


)

VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")

fung_All_sites_soil_meta_color=c("ph"="#7C001D",
                                 "Geographic"="black",
                                 "dry_matter_yield_mg_ha_mean"="#66BD63")




fung_All_sites_soil_meta_names=c("ph"="Soil pH",
                                 "Geographic"="Spatial factor",
                                 "dry_matter_yield_mg_ha_mean"="Subplot yield")




ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = fung_All_sites_soil_meta_color,name=NULL,labels=fung_All_sites_soil_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_x_long,
        GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long)



head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub)

MLE_pH_points= data.frame(pH_p=unique(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$ph))
summary(MLE_pH_points)

(MLE_fung_soil_pH_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="ph"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_soil_meta_color,name=NULL,labels=fung_All_sites_soil_meta_names)+ylim(c(-0.005,0.55))+
    geom_segment(data=MLE_pH_points,aes(x=pH_p,xend = pH_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil pH",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_blank()))



(MLE_fung_soil_spat_gdm=ggplot(data=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="Geographic"))+
    geom_smooth(aes(x=x_predictors,y=y_values,color=metadata),linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_soil_meta_color,name=NULL,labels=fung_All_sites_soil_meta_names)+
    geom_rect(data = All_sites_sq,aes(xmin=min_dis,xmax = max_dis,group=site_comp),
              size=2,color="black",fill="black",ymin=-0.005,ymax=-Inf)+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)",limits = c(-0.005,0.55))+
    labs(x="Distance (m)")+
    theme(legend.position = "none"))

MLE_dry_matter_yield_mg_ha_mean_points= data.frame(dry_matter_yield_mg_ha_mean_p=unique(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_Space_metadata_abiotic_sub$dry_matter_yield_mg_ha_mean))
summary(MLE_dry_matter_yield_mg_ha_mean_points)

(MLE_fung_soil_dry_matter_yield_mg_ha_mean_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="dry_matter_yield_mg_ha_mean"),
                                                      aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_soil_meta_color,name=NULL,labels=fung_All_sites_soil_meta_names)+ylim(c(-0.005,0.55))+
    geom_segment(data=MLE_dry_matter_yield_mg_ha_mean_points,aes(x=dry_matter_yield_mg_ha_mean_p,
                                                                 xend = dry_matter_yield_mg_ha_mean_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Subplot yield (Mg/ha)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))

(LEGEND_MLE_fung_soil_dry_matter_yield_mg_ha_mean_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_NS_long,metadata=="dry_matter_yield_mg_ha_mean"),
                                                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_All_sites_soil_meta_color,name=NULL,labels=fung_All_sites_soil_meta_names)+ylim(c(-0.005,0.55))+
    geom_segment(data=MLE_dry_matter_yield_mg_ha_mean_points,aes(x=dry_matter_yield_mg_ha_mean_p,
                                                                 xend = dry_matter_yield_mg_ha_mean_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Subplot yield (Mg/ha)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme())

plot_grid(MLE_fung_soil_spat_gdm,MLE_fung_soil_pH_gdm,
          MLE_fung_soil_dry_matter_yield_mg_ha_mean_gdm,
          get_legend(LEGEND_MLE_fung_soil_dry_matter_yield_mg_ha_mean_gdm),rel_widths = c(1.2,1),
          labels = c("a)","b)","c)",""),label_size = 36,label_x = c(0,-0.1))




ggplot()+ggtitle("Soil Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Soil Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(MLE_fung_soil_spat_gdm,MLE_fung_soil_pH_gdm,
                    MLE_fung_soil_dry_matter_yield_mg_ha_mean_gdm,
                    get_legend(LEGEND_MLE_fung_soil_dry_matter_yield_mg_ha_mean_gdm),rel_widths = c(1.2,1),
                    labels = c("a)","b)","c)",""),label_size = 36,label_x = c(0,-0.1)),nrow = 2,
          rel_heights = c(0.075,1))


#All Sites: Graphing the combined GDM graphs####

All_sample_GDM_color=c("Geographic"="black",#11111
                       # "MET measurement"="#313695"
                       "past_7d_rain"="#4575B4",#11
                       "percent_soil_moisture_dry_weight"="#313695",#1
                       "soil_temp_1_avg_mean"="#74ADD1",#1
                       #"Soil nutrients"="#A50026"
                       "ph"="#7C001D",#11111
                       "k_ppm"="#F46D43",#11111
                       "ca_ppm"="#FEE090",#111
                       "TON"="#BC2A22",#1
                       "p_ppm"="#D73027",#11 
                       #"Plant traits"="#006837"
                       "core_root_mass_subplot"="#006837",#1
                       "dry_matter_yield_mg_ha_mean"="#66BD63",#1
                       "collectionDate_N"="darkgrey"#1
)
brewer.pal(11,"RdYlGn")


All_sample_GDM_names=c("Geographic"="Spatial factor",#11111
                       "past_7d_rain"="Seven day\nrain accumulation",#11
                       "percent_soil_moisture_dry_weight"="Gravimetric\nsoil moisture",#1
                       "soil_temp_1_avg_mean"="Soil temp (24h avg)",#1
                       "ph"="Soil pH",#11111
                       "k_ppm"="Soil K",#11111
                       "ca_ppm"="Soil Ca",#111
                       "TON"="Soil TON",#1
                       "p_ppm"="Soil P",#11 
                       "core_root_mass_subplot"="Root biomass",#1
                       "dry_matter_yield_mg_ha_mean"="Subplot yield",#1
                       "collectionDate_N"="Temporal factor"#1
)




All_MLE_sample_GDM_color=All_sample_GDM_color[names(All_sample_GDM_color)%in%
                                                unique(c(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long$metadata))]

All_MLE_sample_GDM_names=All_sample_GDM_names[names(All_sample_GDM_names)%in%
                                                unique(c(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long$metadata))]



(MLE_bact_root_sub_GDM=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    theme_cowplot(font_size = 24)+labs(x=NULL,y="Root\nPartial ecological distance\n(Bray-Curtis)")+
    ggtitle("Bacteria")+annotation_custom(grobTree(textGrob(paste("p < 0.001\n59.7%"), x=0.1,  y=0.90, hjust=0,
                                                            gp=gpar(col="black", fontsize=26))))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 36)))

(MLE_bact_soil_sub_GDM=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n82.5%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Soil\nPartial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))

(MLE_fung_root_sub_GDM=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n46.7%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x=NULL,y=NULL)+ggtitle("Fungi")+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 36)))


(MLE_fung_soil_sub_GDM=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n67.6%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.position = "none"))


All_MLE_sample_GDM_order=c("Geographic","ph","TON","k_ppm","ca_ppm",
                           "past_7d_rain",
                           "core_root_mass_subplot","dry_matter_yield_mg_ha_mean")

(MLE_fung_soil_sub_GDM_legend=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long,
                                     aes(x=x_predictors,y=y_values,color=metadata))+
    geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names,limits=All_MLE_sample_GDM_order)+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.text.align = 0,legend.key.size = unit(2, 'cm')))

(four_panel_GDM_MLE=plot_grid(plot_grid(MLE_bact_root_sub_GDM,MLE_fung_root_sub_GDM,
                                        MLE_bact_soil_sub_GDM,MLE_fung_soil_sub_GDM,labels = c("a)","b)","c)","d)"),
                                        ncol = 2, rel_widths = c(1,1), label_size = 32,
                                        label_x = c(0,0.08,0,0.08),label_y = c(1,1,1.05,1.05),
                                        align = "v"),get_legend(MLE_fung_soil_sub_GDM_legend), rel_widths = c(1,0.3)))


ggsave(four_panel_GDM_MLE,
       filename = "GDM_loess_metadata_MLE_four_panel.png",path = here::here("Manuscript","GDM_metadata_figs"),width = 20,height =12)

ggsave(four_panel_GDM_MLE,
       filename = "GDM_loess_metadata_MLE_four_panel.svg",path = here::here("Manuscript","GDM_metadata_figs"),width = 20,height =12)



#All Sites: Graphing the Euler diagrams of partitioning of deviance explained####



GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART$community=rep("Fungi")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART$niche=rep("Soil")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART2=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART[9:15,]
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART$community=rep("Bacteria")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART$niche=rep("Soil")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART2=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART[9:15,]
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART$community=rep("Fungi")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART$niche=rep("Root")
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART2=
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART[9:15,]
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART$experiment=rep("MLE")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART$community=rep("Bacteria")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART$niche=rep("Root")
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART2=
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART[9:15,]


GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART=
  rbind(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_TabForm_PART2,
        GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_sub_TabForm_PART2,
        GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_sub_TabForm_PART2,
        GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_sub_TabForm_PART2)

unique(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART$VARIABLE_SET)
length(unique(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART$VARIABLE_SET))

GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART$VARIABLE_SET_fix=
  GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART$VARIABLE_SET

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

GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART$VARIABLE_SET_euler=
  GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART$VARIABLE_SET

GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART <- GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART %>% 
  mutate(VARIABLE_SET_euler=mapvalues(VARIABLE_SET_euler,
                                      from = VARIABLE_SET_order_f,
                                      to = euler_VARIABLE_SET_names))


MET_SOIL_SPAT_alone_color=c("MET"="#313695","Soil"="#A50026","Spatial"="darkgrey")
PLANT_SOIL_SPAT_alone_color=c("Plant"="#006837","Soil"="#A50026","Spatial"="darkgrey")
VARIABLE_SET_euler_alone_color=c("Spatial"="darkgrey",
                                 "Soil"="#A50026",
                                 "MET"="#313695",
                                 "Plant"="#006837")

MLE_Root_bact_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Root")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Root")$VARIABLE_SET_euler)

MLE_Root_bact_euler_fit=euler(MLE_Root_bact_euler_list)



error_plot(MLE_Root_bact_euler_fit)




(MLE_Root_bact_euler_p=plot(MLE_Root_bact_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = MET_SOIL_SPAT_alone_color))







MLE_Soil_bact_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Soil")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Soil")$VARIABLE_SET_euler)

MLE_Soil_bact_euler_fit=euler(MLE_Soil_bact_euler_list)


(MLE_Soil_bact_euler_p=plot(MLE_Soil_bact_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = MET_SOIL_SPAT_alone_color))

plot_grid(MLE_Root_bact_euler_p,MLE_Soil_bact_euler_p,ncol = 1)


MLE_Root_fung_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Fungi" &niche=="Root")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Fungi" &niche=="Root")$VARIABLE_SET_euler)

MLE_Root_fung_euler_fit=euler(MLE_Root_fung_euler_list)



error_plot(MLE_Root_fung_euler_fit)




(MLE_Root_fung_euler_p=plot(MLE_Root_fung_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = PLANT_SOIL_SPAT_alone_color))







MLE_Soil_fung_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Fungi" &niche=="Soil")$DEVIANCE), 
           subset(GLBRC018_OTU_MMPRNT_All_sites_G5_sub_TabForm_PART,community=="Fungi" &niche=="Soil")$VARIABLE_SET_euler)

MLE_Soil_fung_euler_fit=euler(MLE_Soil_fung_euler_list)

error_plot(MLE_Soil_fung_euler_fit)

(MLE_Soil_fung_euler_p=plot(MLE_Soil_fung_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = PLANT_SOIL_SPAT_alone_color))

(MLE_fung_soil_sub_GDM_legend2=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_sub_gdm_long,
                                      aes(x=x_predictors,y=y_values,color=metadata))+
    geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_MLE_sample_GDM_color,name=NULL,labels=All_MLE_sample_GDM_names,limits=All_MLE_sample_GDM_order)+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.key.size = unit(2, 'cm'),legend.position = "bottom"))

#legend.text.align = 0,
(six_panel_Euler_GDM_MLE=plot_grid(plot_grid(plot_grid(MLE_Root_bact_euler_p,MLE_Soil_bact_euler_p,ncol = 1,
                                                       rel_heights = c(1.4,0.6)),
                                             plot_grid(MLE_bact_root_sub_GDM,MLE_fung_root_sub_GDM,
                                                       MLE_bact_soil_sub_GDM,MLE_fung_soil_sub_GDM,
                                                       labels = c("a)","b)","c)","d)"),
                                                       ncol = 2, label_size = 40,
                                                       label_x = c(0,0.08,0,0.08),label_y = c(1.03,1.03,1.05,1.05),
                                                       align = "v"),plot_grid(MLE_Root_fung_euler_p,MLE_Soil_fung_euler_p,ncol = 1),
                                             rel_widths = c(0.2,1,0.2),ncol = 3),
                                   plot_grid(ggplot+geom_blank(),get_legend(MLE_fung_soil_sub_GDM_legend2),rel_widths = c(0.5,1)), rel_heights = c(1,0.2),ncol = 1))


ggsave(six_panel_Euler_GDM_MLE,
       filename = "MLE_Euler_GDM_Percentage_deviance_plot.png",path = here::here("Manuscript","GDM_metadata_figs"),width = 28,height =14)
ggsave(six_panel_Euler_GDM_MLE,
       filename = "MLE_Euler_GDM_Percentage_deviance_plot.svg",path = here::here("Manuscript","GDM_metadata_figs"),width = 28,height =14)


