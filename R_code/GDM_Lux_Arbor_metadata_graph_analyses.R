
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

"%w/o%" <- function(x,y)!('%in%'(x,y))
Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "English")

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

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)



bacterial_root_sampling_date=data.frame(sampleDate=c(unique(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,Root_soil=="Root")$collectionDate)),
                                        Root_soil=rep("Root"))
bacterial_soil_sampling_date=data.frame(sampleDate=c(unique(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,Root_soil=="Soil")$collectionDate)),
                                        Root_soil=rep("Soil"))

sampling_date=rbind(bacterial_root_sampling_date,
                    bacterial_soil_sampling_date[bacterial_soil_sampling_date$sampleDate%w/o%bacterial_root_sampling_date$sampleDate,])

#####GDM: LUX Abiotic and Biotic predictors####



MMPRNT_2018.metadata_mapp_coverage_LUX=subset(MMPRNT_2018.metadata_mapp_coverage,siteID=="LUX"&plotType=="G5")
nrow(MMPRNT_2018.metadata_mapp_coverage_LUX)
#360


colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_LUX))


MMPRNT_2018.metadata_mapp_LUX=MMPRNT_2018.metadata_mapp_coverage_LUX
MMPRNT_2018.metadata_mapp_LUX[,c("event_age","pH_raw_MMPRNT","pH_MMPRNT_subplot")]=NULL

colnames(MMPRNT_2018.metadata_mapp_LUX)

#We are missing MET from 3/19/2018, let's remove it

MMPRNT_2018.metadata_mapp_LUX_lim_date=subset(MMPRNT_2018.metadata_mapp_LUX,collectionDate!="3/19/2018")

colSums(is.na(MMPRNT_2018.metadata_mapp_LUX_lim_date))


#Metadata Extraction 



write.csv(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","percent_soil_moisture_dry_weight","plant.height",
                                           "core_root_mass_subplot","ugN_NH4_g_dry_soil_na",
                                           "ugN_NO3_g_dry_soil_na","vwc_avg_mean","soil_temp_1_avg_mean","airtc_avg_mean","past_7d_rain")]%>%
            group_by(collectionDate)%>%summarise_all(~mean(.,na.rm =T)),
          file = here::here("R_files","ignorable_files","LUX_ARBOR_Abiotic_Biotic_summary.csv"))






print(MMPRNT_2018.metadata_mapp_LUX_lim_date|>group_by(plotRep,FertStatus,collectionDate)|>summarise(biomass=sum(is.na(total.shoot.dry.weight)),
                                                                                                     leaf_trait=sum(is.na(specific.leaf.area)),
                                                                                                     met=sum(is.na(soil_temp_1_avg_mean)),
                                                                                                     root=sum(is.na(core_root_mass_subplot)),
                                                                                                     TOC_N=sum(is.na(TOC_TON_ratio))),n =112)


MMPRNT_2018.metadata_mapp_LUX_lim=data.frame(MMPRNT_2018.metadata_mapp_LUX_lim_date[,c("SampleID","collectionDate",
                                                                                       "percent_soil_moisture_dry_weight",
                                                                                       "dry_matter_yield_mg_ha_mean",
                                                                                       "total.shoot.dry.weight",
                                                                                       "specific.leaf.area","core_root_mass_subplot",
                                                                                       "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na",
                                                                                       "soil_temp_1_avg_mean","past_7d_rain",
                                                                                       "ph","p_ppm","k_ppm","ca_ppm","TOC","TON","TOC_TON_ratio")])

colSums(is.na(MMPRNT_2018.metadata_mapp_LUX_lim))

MMPRNT_2018.metadata_mapp_LUX_lim|>group_by(collectionDate)|>summarise(biomass=sum(is.na(total.shoot.dry.weight)),
                                                                       leaf_trait=sum(is.na(specific.leaf.area)),
                                                                       met=sum(is.na(soil_temp_1_avg_mean)),
                                                                       root=sum(is.na(core_root_mass_subplot)),
                                                                       TOC_N=sum(is.na(TOC_TON_ratio)))



# All missing shoot weight are due to a lack of emerged switchgrass, let's add zeros

MMPRNT_2018.metadata_mapp_LUX_lim$total.shoot.dry.weight=ifelse(MMPRNT_2018.metadata_mapp_LUX_lim$collectionDate=="11/5/2018"|
                                                                  MMPRNT_2018.metadata_mapp_LUX_lim$collectionDate=="4/30/2018",0,
                                                                MMPRNT_2018.metadata_mapp_LUX_lim$total.shoot.dry.weight)


MMPRNT_2018.metadata_mapp_LUX_lim$specific.leaf.area=ifelse(MMPRNT_2018.metadata_mapp_LUX_lim$collectionDate=="11/5/2018"|
                                                              MMPRNT_2018.metadata_mapp_LUX_lim$collectionDate=="4/30/2018",0,
                                                            MMPRNT_2018.metadata_mapp_LUX_lim$specific.leaf.area)

colSums(is.na(MMPRNT_2018.metadata_mapp_LUX_lim))
ncol(MMPRNT_2018.metadata_mapp_LUX_lim)


write.csv(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","percent_soil_moisture_dry_weight",
                                           "total.shoot.dry.weight",
                                           "specific.leaf.area","core_root_mass_subplot",
                                           "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na",
                                           "soil_temp_1_avg_mean","past_7d_rain",
                                           "TOC","TON","TOC_TON_ratio")]%>%
            group_by(collectionDate)%>%summarise_all(~mean(.,na.rm =T),~sd(.,na.rm =T)/sqrt(n())),
          file = here::here("R_files","ignorable_files","LUX_ARBOR_Abiotic_Biotic_summary_full.csv"))


#Impute the missing data
summary(MMPRNT_2018.metadata_mapp_LUX_lim)
row.names(MMPRNT_2018.metadata_mapp_LUX_lim)=MMPRNT_2018.metadata_mapp_LUX_lim$SampleID

MMPRNT_2018.metadata_mapp_LUX_lim$collectionDate=as.numeric(as.Date(MMPRNT_2018.metadata_mapp_LUX_lim$collectionDate,format = "%m/%d/%Y"))
set.seed(2022)
MMPRNT_2018.metadata_mapp_LUX_lim.imp <- missForest(MMPRNT_2018.metadata_mapp_LUX_lim[,-1])


MMPRNT_2018.metadata_mapp_LUX_lim.imp$ximp

#check imputation error
MMPRNT_2018.metadata_mapp_LUX_lim.imp$OOBerror
#     NRMSE 
#0.002044304 

MMPRNT_2018.metadata_mapp_LUX_lim_imp=MMPRNT_2018.metadata_mapp_LUX_lim.imp$ximp
head(MMPRNT_2018.metadata_mapp_LUX_lim_imp)

MMPRNT_2018.metadata_mapp_LUX_lim_imp$sampleID_long=
  paste("MMPRNT",str_sub(row.names(MMPRNT_2018.metadata_mapp_LUX_lim_imp),start = 7,end = 11),sep = "-")
head(MMPRNT_2018.metadata_mapp_LUX_lim_imp)

#GDM Lux Bacterial Roots####



GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data=sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)[,c("sampleID_long","SampleID","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data)
#142
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data)

colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data)[colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata=
  data.frame(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic=merge(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata,
                                                                  MMPRNT_2018.metadata_mapp_LUX_lim_imp, by = "sampleID_long")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic$collectionDate_N=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic$collectionDate
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic[,c("sampleID_long","collectionDate")]=NULL
dim(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic)
#142  20
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root,method = "bray")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF),GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site[1:10,1:10]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site)
#142
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic)
nrow(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                              XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                              predData=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   20.599


#Filter out factors with less than 0.01 Sum of coefficients

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord","percent_soil_moisture_dry_weight",
                                                                 "collectionDate_N","soil_temp_1_avg_mean","ca_ppm",
                                                                 "dry_matter_yield_mg_ha_mean","total.shoot.dry.weight", 
                                                                 "ugN_NO3_g_dry_soil_na","k_ppm","core_root_mass_subplot",
                                                                 "p_ppm","ph","past_7d_rain","specific.leaf.area")]
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub1_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                   XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                   predData=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub1_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub1_TabForm, geo=FALSE)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub1_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub1_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   20.558
gdm.varImp(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub1_TabForm, geo = F, predSelect = T,nPerm = 100)




#Important factors
#Final set of predictors returned: 
#collectionDate_N
#ca_ppm

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                 "collectionDate_N","ca_ppm")]
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                  XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                  predData=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm, geo=FALSE)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   15.118

gdm.varImp_MOD(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm, geo = F, nPerm = 100)
#                           All predictors
#Model deviance                    127.070
#Percent deviance explained         15.118
#Model p-value                       0.000
#Fitted permutations                95.000


LUX_bact_R_meta_varSet <- vector("list", 2)

names(LUX_bact_R_meta_varSet)=c("Time","soil")
LUX_bact_R_meta_varSet$Time=c("collectionDate_N")
LUX_bact_R_meta_varSet$soil=c("ca_ppm")
#LUX_bact_R_meta_varSet$MET=c("soil_temp_1_avg_mean")
summary(LUX_bact_R_meta_varSet)
#Variance partitioning
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART=
  gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm,varSets = LUX_bact_R_meta_varSet,partSpace= F)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART)



GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod)$x,
                                                                    method = "range"))

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_x_long=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_x_long,
        GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_long)




VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")

bact_LUX_root_meta_color=c("collectionDate_N"="darkgrey",
                           "ca_ppm"="#FEE090")




bact_LUX_root_meta_names=c("collectionDate_N"="Temporal factor",
                           "ca_ppm"="Soil Ca")




ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = bact_LUX_root_meta_color,name=NULL,labels=bact_LUX_root_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")






GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_NS_x_long=data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_NS_x_long,
        GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_NS_long)




VARIABLE_SET_alone_color=c("Spatial factor"="black",
                           "Soil nutrients"="#A50026",
                           "MET measurement"="#313695",
                           "Plant traits"="#006837")

bact_LUX_root_meta_color=c("collectionDate_N"="darkgrey",
                           "ca_ppm"="#FEE090")




bact_LUX_root_meta_names=c("collectionDate_N"="Temporal factor",
                           "ca_ppm"="Soil Ca")

ca_ppm_points= data.frame(ca_ppm_p=unique(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub$ca_ppm))
summary(ca_ppm_points)

(Lux_bact_root_ca_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="ca_ppm"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_root_meta_color,name=NULL,labels=bact_LUX_root_meta_names)+ylim(c(-0.005,0.15))+
    geom_segment(data=ca_ppm_points,aes(x=ca_ppm_p,xend = ca_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil Ca (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))
#



(Lux_bact_root_temporal_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="collectionDate_N"),
                                   aes(x=as.Date(x_predictors,origin="1970-01-01"),y=y_values,color=metadata))+
    geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_root_meta_color,name=NULL,labels=bact_LUX_root_meta_names)+
    geom_segment(data = subset(sampling_date,Root_soil=="Root"),aes(x=as.Date(sampleDate,format = "%m/%d/%Y"),
                                                                    xend = as.Date(sampleDate,format = "%m/%d/%Y"),y=-0.005,yend=-Inf),size=2,color="black")+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)", limits=c(-0.005,0.15))+
    scale_x_date(name="Collection date")+
    theme(legend.position = "none"))
ggplot()+ggtitle("Root Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Root Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(Lux_bact_root_temporal_gdm,Lux_bact_root_ca_gdm,rel_widths = c(1,1.3),
                    labels = c("a)","b)"),label_size = 36,label_x = c(0,-0.075)),nrow = 2,
          rel_heights = c(0.075,1))

#GDM Lux Bacterial Soil matching root####


#Lux Bacteria

unique(sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)$collectionDate)
#"9/17/2018" "5/29/2018" "10/3/2018" "8/20/2018" "7/30/2018" "6/25/2018"

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR=subset_samples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil,
                                                       collectionDate=="9/17/2018"|collectionDate=="5/29/2018"|
                                                         collectionDate=="10/3/2018"|collectionDate=="8/20/2018"|
                                                         collectionDate=="7/30/2018"|collectionDate=="6/25/2018")

nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)
#140
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)>0,GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)
ntaxa(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)
#18992

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data=sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)[,c("sampleID_long","SampleID","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data)
#140
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data)

colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data)[colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata=
  data.frame(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic=merge(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata,
                                                                     MMPRNT_2018.metadata_mapp_LUX_lim_imp, by = "sampleID_long")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic$collectionDate_N=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic$collectionDate
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic[,c("sampleID_long","collectionDate")]=NULL
dim(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)
#140  20
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR,method = "bray")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF),GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)
nrow(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF_site, 3, siteColumn="site", 
                                                                 XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                 predData=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   29.76

#Filter factors below 0.01

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                    "core_root_mass_subplot","percent_soil_moisture_dry_weight",
                                                                    "ugN_NO3_g_dry_soil_na","ca_ppm","dry_matter_yield_mg_ha_mean","ugN_NH4_g_dry_soil_na",
                                                                    "ph","p_ppm","k_ppm","specific.leaf.area","total.shoot.dry.weight",
                                                                    "collectionDate_N","TOC")]
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub1)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub1_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF_site, 3, siteColumn="site", 
                                                                      XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                      predData=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub1)

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub1_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub1_gdmTabMod)
#Percent Deviance Explained:   29.685
gdm.varImp(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub1_TabForm, geo = T, predSelect = T,nPerm = 100)

#Final set of predictors returned: with all factors
#Geographic
#dry_matter_yield_mg_ha_mean


#Final set of predictors returned: 
#Geographic
#percent_soil_moisture_dry_weight
#ca_ppm
#dry_matter_yield_mg_ha_mean


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                    "percent_soil_moisture_dry_weight","ca_ppm",
                                                                    "dry_matter_yield_mg_ha_mean")]
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_bray_DF_site, 3, siteColumn="site", 
                                                                     XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                     predData=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub)


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm, geo=T)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   25.397

#gdm.varImp(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm, geo = T, nPerm = 100)
gdm.varImp_MOD(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm, geo = T, nPerm = 100)
#                           All predictors
#Model deviance                     59.053
#Percent deviance explained         25.397
#Model p-value                       0.000
#Fitted permutations                99.000


LUX_bact_S_MR_meta_varSet <- vector("list", 2)

names(LUX_bact_S_MR_meta_varSet)=c("MET","soil_plant")
LUX_bact_S_MR_meta_varSet$MET=c("percent_soil_moisture_dry_weight")
LUX_bact_S_MR_meta_varSet$soil_plant=c("ca_ppm", "dry_matter_yield_mg_ha_mean")

summary(LUX_bact_S_MR_meta_varSet)
#Variance partitioning
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART=
  gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm,varSets = LUX_bact_S_MR_meta_varSet,partSpace= T)
summary(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART)

LUX_bact_S_MR_meta_varSet2 <- vector("list", 3)
names(LUX_bact_S_MR_meta_varSet2)=c("MET","soil","plant")
LUX_bact_S_MR_meta_varSet2$MET=c("percent_soil_moisture_dry_weight")
LUX_bact_S_MR_meta_varSet2$soil=c("ca_ppm")
LUX_bact_S_MR_meta_varSet2$plant=c("dry_matter_yield_mg_ha_mean")

summary(LUX_bact_S_MR_meta_varSet2)
#Variance partitioning
gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm,varSets = LUX_bact_S_MR_meta_varSet2,partSpace= F)



GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$x,
                                                                       method = "range"))

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_x_long=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_x_long,
        GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_long)

bact_LUX_soil_MR_meta_color=c("Geographic"="black",
                              "percent_soil_moisture_dry_weight"="#313695",
                              "ca_ppm"="#FEE090",
                              "dry_matter_yield_mg_ha_mean"="#66BD63")

bact_LUX_soil_MR_meta_names=c("Geographic"="Spatial factor",
                              "percent_soil_moisture_dry_weight"="Gravimetric\nsoil moisture",
                              "ca_ppm"="Soil Ca",
                              "dry_matter_yield_mg_ha_mean"="Subplot yield")





ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = bact_LUX_soil_MR_meta_color,name=NULL,labels=bact_LUX_soil_MR_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")





GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_x_long,
        GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long)
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long)

bact_LUX_soil_MR_meta_color=c("Geographic"="black",
                              "percent_soil_moisture_dry_weight"="#313695",
                              "ca_ppm"="#FEE090",
                              "dry_matter_yield_mg_ha_mean"="#66BD63")



bact_LUX_soil_MR_meta_names=c("Geographic"="Spatial factor",
                              "percent_soil_moisture_dry_weight"="Gravimetric\nsoil moisture",
                              "ca_ppm"="Soil Ca",
                              "dry_matter_yield_mg_ha_mean"="Subplot yield")


head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub)


percent_soil_moisture_dry_weight_points= data.frame(percent_soil_moisture_dry_weight_p=unique(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub$percent_soil_moisture_dry_weight))
summary(ca_ppm_points)

(Lux_bact_soil_moist_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="percent_soil_moisture_dry_weight"),
                                aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_soil_MR_meta_color,name=NULL,labels=bact_LUX_soil_MR_meta_names)+ylim(c(-0.005,0.08))+
    geom_segment(data=percent_soil_moisture_dry_weight_points,aes(x=percent_soil_moisture_dry_weight_p,
                                                                  xend = percent_soil_moisture_dry_weight_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Gravimetric soil moisture",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis <- 
  matrixConvert((dist(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata[,c("UTM_Lat_Cord", "UTM_Lon_Cord")], method = "euclidean")), 
                colname = c("sample1", "sample2", "distance"))#total distance 

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2=merge(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis,
                                                  sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)[,c("plotRep")],
                                                  by.x = "sample1",by.y="row.names")
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2)
colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2)[colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2)=="plotRep"]="s1_plotRep"


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m=merge(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2,
                                                  sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR)[,c("plotRep")],
                                                  by.x = "sample2",by.y="row.names")

head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m)
colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m)[colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m)=="plotRep"]="s2_plotRep"


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m$s1_s2_plotRep=with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m,
                                                                  interaction(s1_plotRep,s2_plotRep))
unique(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m$s1_s2_plotRep)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m$plot_comp=with(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m,
                                                              ifelse(s1_plotRep==s2_plotRep|s1_s2_plotRep=="4.3"|
                                                                       s1_s2_plotRep=="3.4","within",
                                                                     ifelse(s1_s2_plotRep=="1.4"|s1_s2_plotRep=="4.1"|
                                                                              s1_s2_plotRep=="1.3"|s1_s2_plotRep=="3.1",
                                                                            "long","short")))
unique(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m$plot_comp)

bacterial_sampling_dist=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m%>%group_by(plot_comp)%>%
  summarise(max_dis=max(distance),min_dis=min(distance))


(Lux_bact_soil_spat_gdm=ggplot(data=subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="Geographic"))+
    geom_smooth(aes(x=x_predictors,y=y_values,color=metadata),linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_soil_MR_meta_color,name=NULL,labels=bact_LUX_soil_MR_meta_names)+
    geom_rect(data = bacterial_sampling_dist,aes(xmin=min_dis,xmax = max_dis,group=plot_comp),
              size=2,color="black",fill="black",ymin=-0.005,ymax=-Inf)+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)", limits=c(-0.005,0.08))+
    labs(x="Distance (m)")+
    theme(legend.position = "none"))

dry_matter_yield_mg_ha_mean_points= data.frame(dry_matter_yield_mg_ha_mean_p=unique(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub$dry_matter_yield_mg_ha_mean))
summary(ca_ppm_points)

(Lux_bact_yield_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="dry_matter_yield_mg_ha_mean"),
                           aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_soil_MR_meta_color,name=NULL,labels=bact_LUX_soil_MR_meta_names)+ylim(c(-0.005,0.08))+
    geom_segment(data=dry_matter_yield_mg_ha_mean_points,aes(x=dry_matter_yield_mg_ha_mean_p,
                                                             xend = dry_matter_yield_mg_ha_mean_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Subplot yield (Mg/ha)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))

ca_ppm_points= data.frame(ca_ppm_p=unique(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub$ca_ppm))
summary(ca_ppm_points)

(Lux_bact_soil_ca_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="ca_ppm"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_soil_MR_meta_color,name=NULL,labels=bact_LUX_soil_MR_meta_names)+ylim(c(-0.005,0.08))+
    geom_segment(data=ca_ppm_points,aes(x=ca_ppm_p,xend = ca_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil Ca (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))


(LEGEND_Lux_bact_soil_ca_gdm=ggplot(subset(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="ca_ppm"),
                                    aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = bact_LUX_soil_MR_meta_color,name=NULL,labels=bact_LUX_soil_MR_meta_names)+ylim(c(-0.005,0.08))+
    geom_segment(data=ca_ppm_points,aes(x=ca_ppm_p,xend = ca_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil Ca (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))

plot_grid(plot_grid(Lux_bact_soil_spat_gdm,Lux_bact_soil_moist_gdm,
                    Lux_bact_yield_gdm,Lux_bact_soil_ca_gdm,rel_widths = c(1.3,1),
                    labels = c("a)","b)","c)","d)"),label_size = 36,label_x = c(0,-0.075)),
          get_legend(LEGEND_Lux_bact_soil_ca_gdm),rel_widths = c(1,0.3))



ggplot()+ggtitle("Soil Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Soil Bacteria")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(plot_grid(Lux_bact_soil_spat_gdm,Lux_bact_soil_moist_gdm,
                              Lux_bact_yield_gdm,Lux_bact_soil_ca_gdm,rel_widths = c(1.3,1),
                              labels = c("a)","b)","c)","d)"),label_size = 36,label_x = c(0,-0.1)),
                    get_legend(LEGEND_Lux_bact_soil_ca_gdm),rel_widths = c(1,0.3)),nrow = 2,
          rel_heights = c(0.075,1))


#GDM Lux Fungal Roots####



GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data=sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)[,c("sampleID_long","sampleID_seq","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data)
#141
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data)

colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data)[colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata=
  data.frame(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic=merge(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata,
                                                                  MMPRNT_2018.metadata_mapp_LUX_lim_imp, by = "sampleID_long")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic$collectionDate_N=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic$collectionDate
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic[,c("sampleID_long","collectionDate")]=NULL
dim(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic)
#141  20
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root,method = "bray")
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF),GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic)
nrow(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                              XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                              predData=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   17.161
#gdm.varImp(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_TabForm, geo = T, predSelect = T,nPerm = 100)



#Filter the factor with less than 0.01
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                 "percent_soil_moisture_dry_weight","soil_temp_1_avg_mean",
                                                                 "specific.leaf.area",
                                                                 "k_ppm","core_root_mass_subplot",
                                                                 "ca_ppm","TOC","collectionDate_N","TOC_TON_ratio","ph",
                                                                 "total.shoot.dry.weight")]
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub1_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                   XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                   predData=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub1)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub1_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub1_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub1_gdmTabMod)
#17.161

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub1_TabForm, geo = T, predSelect = T,nPerm = 100)

#Final set of predictors returned: 
#Geographic
#percent_soil_moisture_dry_weight
#soil_temp_1_avg_mean
#k_ppm
#ph




GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                 "soil_temp_1_avg_mean", "k_ppm","ph",
                                                                 "percent_soil_moisture_dry_weight")]
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", 
                                                                  XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                  predData=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm, geo=T)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   15.502

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm, geo = T, nPerm = 100)
#Geographic                                8.173
#soil_temp_1_avg_mean                     15.434
#k_ppm                                    10.724
#ph                                        7.769
#percent_soil_moisture_dry_weight         29.914

LUX_fung_R_meta_varSet <- vector("list", 2)

names(LUX_fung_R_meta_varSet)=c("MET","soil")
LUX_fung_R_meta_varSet$MET=c("percent_soil_moisture_dry_weight","soil_temp_1_avg_mean")
LUX_fung_R_meta_varSet$soil=c("k_ppm","ph")
summary(LUX_fung_R_meta_varSet)
#Variance partitioning
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART=
  gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm,varSets = LUX_fung_R_meta_varSet,partSpace= T)



GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod)$x,
                                                                    method = "range"))

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_x_long=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_x_long,
        GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_long)

fung_LUX_root_meta_color=c("percent_soil_moisture_dry_weight"="#313695",
                           "soil_temp_1_avg_mean"="#74ADD1",
                           "Geographic"="black",
                           "k_ppm"="#F46D43",
                           "ph"="#7C001D")




fung_LUX_root_meta_names=c("soil_temp_1_avg_mean"="Soil temp (24h avg)",
                           "percent_soil_moisture_dry_weight"="Gravimetric\nsoil moisture",
                           "Geographic"  = "Spatial factor",
                           "k_ppm"="Soil K",
                           "ph"="Soil pH")




ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")



GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_x_long,
        GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long)


fung_LUX_root_meta_color=c("percent_soil_moisture_dry_weight"="#313695",
                           "soil_temp_1_avg_mean"="#74ADD1",
                           "Geographic"="black",
                           "k_ppm"="#F46D43",
                           "ph"="#7C001D")




fung_LUX_root_meta_names=c("soil_temp_1_avg_mean"="Soil temp (24h avg)",
                           "percent_soil_moisture_dry_weight"="Gravimetric\nsoil moisture",
                           "Geographic"  = "Spatial factor",
                           "k_ppm"="Soil K",
                           "ph"="Soil pH")





head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub)


F_percent_soil_moisture_dry_weight_points= data.frame(percent_soil_moisture_dry_weight_p=unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub$percent_soil_moisture_dry_weight))
summary(F_percent_soil_moisture_dry_weight_points)

(Lux_fung_root_moist_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="percent_soil_moisture_dry_weight"),
                                aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+ylim(c(-0.005,0.48))+
    geom_segment(data=F_percent_soil_moisture_dry_weight_points,aes(x=percent_soil_moisture_dry_weight_p,
                                                                    xend = percent_soil_moisture_dry_weight_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Gravimetric soil moisture",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))

(Lux_fung_root_spat_gdm=ggplot(data=subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="Geographic"))+
    geom_smooth(aes(x=x_predictors,y=y_values,color=metadata),linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+
    geom_rect(data = bacterial_sampling_dist,aes(xmin=min_dis,xmax = max_dis,group=plot_comp),
              size=2,color="black",fill="black",ymin=-0.005,ymax=-Inf)+
    theme_cowplot(font_size = 24)+scale_y_continuous(name="Partial ecological distance\n(Bray-Curtis)", limits=c(-0.005,0.48))+
    labs(x="Distance (m)")+
    theme(legend.position = "none"))

soil_temp_1_avg_mean_points= data.frame(soil_temp_1_avg_mean_p=unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub$soil_temp_1_avg_mean))
summary(soil_temp_1_avg_mean_points)

(Lux_fung_root_soil_temp_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="soil_temp_1_avg_mean"),
                                    aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+ylim(c(-0.005,0.48))+
    geom_segment(data=soil_temp_1_avg_mean_points,aes(x=soil_temp_1_avg_mean_p,
                                                      xend = soil_temp_1_avg_mean_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil temp (24h avg. C)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))

k_ppm_points= data.frame(k_ppm_p=unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub$k_ppm))
summary(k_ppm_points)

(Lux_fung_root_K_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="k_ppm"),
                            aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+ylim(c(-0.005,0.48))+
    geom_segment(data=k_ppm_points,aes(x=k_ppm_p,xend = k_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil K (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))


ph_points= data.frame(ph_p=unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_Space_metadata_abiotic_sub$ph))
summary(ph_points)

(Lux_fung_root_pH_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="ph"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+ylim(c(-0.005,0.48))+
    geom_segment(data=ph_points,aes(x=ph_p,xend = ph_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil pH",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))


(LEGEND_Lux_fung_root_pH_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_NS_long,metadata=="ph"),
                                    aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_root_meta_color,name=NULL,labels=fung_LUX_root_meta_names)+ylim(c(-0.005,0.48))+
    geom_segment(data=ph_points,aes(x=ph_p,xend = ph_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil pH",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))

plot_grid(plot_grid(Lux_fung_root_moist_gdm,Lux_fung_root_soil_temp_gdm,
                    get_legend(LEGEND_Lux_fung_root_pH_gdm),
                    Lux_fung_root_spat_gdm,Lux_fung_root_K_gdm,Lux_fung_root_pH_gdm,rel_widths = c(1.3,1,1),
                    labels = c("a)","b)","","c)","d)","e)"),label_size = 36,label_x = c(0,-0.1,-0.1)))



ggplot()+ggtitle("Root Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Root Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(plot_grid(Lux_fung_root_moist_gdm,Lux_fung_root_soil_temp_gdm,
                              get_legend(LEGEND_Lux_fung_root_pH_gdm),
                              Lux_fung_root_spat_gdm,Lux_fung_root_K_gdm,Lux_fung_root_pH_gdm,rel_widths = c(1.3,1,1),
                              labels = c("a)","b)","","c)","d)","e)"),label_size = 36,label_x = c(0,-0.1,-0.1))),
          nrow = 2,rel_heights = c(0.075,1))


#GDM Lux Fungal Soil matching root####

unique(sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)$collectionDate)
#"9/17/2018" "5/29/2018" "10/3/2018" "8/20/2018" "7/30/2018" "6/25/2018"

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR=subset_samples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil,
                                                       collectionDate=="9/17/2018"|collectionDate=="5/29/2018"|
                                                         collectionDate=="10/3/2018"|collectionDate=="8/20/2018"|
                                                         collectionDate=="7/30/2018"|collectionDate=="6/25/2018")

nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)
#140
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)>0,GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)
ntaxa(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)
#2259

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data=sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)[,c("sampleID_long","sampleID_seq","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data)
#140
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data)

colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data)[colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata=
  data.frame(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_data[,c("sampleID_long","site","UTM_Lat_Cord","UTM_Lon_Cord")])

#Add in the abiotic data 



GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic=merge(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata,
                                                                     MMPRNT_2018.metadata_mapp_LUX_lim_imp, by = "sampleID_long")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic$collectionDate_N=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic$collectionDate
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic[,c("sampleID_long","collectionDate")]=NULL
dim(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)
#140  20
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR,method = "bray")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF)
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF),GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF)
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)
nrow(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF_site, 3, siteColumn="site", 
                                                                 XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                 predData=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   29.268

#Factors with less than 0.01
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub1=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic[,c(
    "site", "UTM_Lat_Cord", "UTM_Lon_Cord",
    "p_ppm", "core_root_mass_subplot","ca_ppm",
    "percent_soil_moisture_dry_weight",
    "dry_matter_yield_mg_ha_mean","TOC","k_ppm",
    "ugN_NO3_g_dry_soil_na","soil_temp_1_avg_mean",
    "TON","ugN_NH4_g_dry_soil_na","collectionDate_N")]

summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub1)

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_TabForm_sub1<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF_site, 3, siteColumn="site", 
                                                                      XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                      predData=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub1)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_gdmTabMod_sub1<- gdm(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_TabForm_sub1, geo=F)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_gdmTabMod_sub1)
#29.266
gdm.varImp(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_TabForm_sub1, geo = F, predSelect = T,nPerm = 100)

#Final set of predictors returned: 
#p_ppm
#ca_ppm



GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic[,c("site","UTM_Lon_Cord","UTM_Lat_Cord",
                                                                    "p_ppm","ca_ppm")]
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_bray_DF_site, 3, siteColumn="site", 
                                                                     XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                     predData=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub)


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm, geo=FALSE)
summary(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:   25.979

gdm.varImp(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm, geo = F, nPerm = 100)
#Error in matrix(NA, 4, nVars, dimnames = list(c("Model deviance", "Percent deviance explained",  : 
#length of 'dimnames' [2] not equal to array extent

LUX_fung_S_MR_meta_varSet <- vector("list", 1)

names(LUX_fung_S_MR_meta_varSet)=c("soil")
LUX_fung_S_MR_meta_varSet$soil=c("p_ppm","ca_ppm")
summary(LUX_fung_S_MR_meta_varSet)
#Variance partitioning
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART=
  gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm,varSets = LUX_fung_S_MR_meta_varSet,partSpace= F)



gdm.varImp_MOD(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm, geo = F, nPerm = 100)
#All predictors
#Model deviance                    273.484
#Percent deviance explained         25.979
#Model p-value                       0.000
#Fitted permutations                99.000


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_x=data.frame(decostand(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$x,
                                                                       method = "range"))

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_x_long=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_x|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_x_long,
        GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long)

fung_LUX_soil_MR_meta_color=c("p_ppm"="#D73027",
                              "ca_ppm"="#FEE090")




fung_LUX_soil_MR_meta_names=c("p_ppm"="Soil P",
                              "ca_ppm"  = "Soil Ca")




ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long,
       aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
  scale_color_manual(values = fung_LUX_soil_MR_meta_color,name=NULL,labels=fung_LUX_soil_MR_meta_names)+
  theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Partial ecological distance\n(Bray-Curtis)")



GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_x_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$x)|>
  pivot_longer(cols = everything(), names_to = "metadata",values_to = "x_predictors")

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long=
  data.frame(isplineExtract(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdmTabMod)$y)|>
  pivot_longer(cols = everything(), names_to = "metadata2",values_to = "y_values")


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long=
  cbind(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_x_long,
        GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_y_long)
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long)


head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub)

ca_ppm_points= data.frame(ca_ppm_p=unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub$ca_ppm))
summary(ca_ppm_points)

(Lux_fung_soil_ca_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="ca_ppm"),
                             aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_soil_MR_meta_color,name=NULL,labels=fung_LUX_soil_MR_meta_names)+ylim(c(-0.005,0.47))+
    geom_segment(data=ca_ppm_points,aes(x=ca_ppm_p,xend = ca_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil Ca (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank()))
#



p_ppm_points= data.frame(p_ppm_p=unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_Space_metadata_abiotic_sub$p_ppm))
summary(p_ppm_points)

(Lux_fung_soil_p_ppm_gdm=ggplot(subset(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_NS_long,metadata=="p_ppm"),
                                aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = fung_LUX_soil_MR_meta_color,name=NULL,labels=fung_LUX_soil_MR_meta_names)+ylim(c(-0.005,0.47))+
    geom_segment(data=p_ppm_points,aes(x=p_ppm_p,xend = p_ppm_p),y=-0.005,yend=-Inf,size=2,color="black")+
    theme_cowplot(font_size = 24)+labs(x="Soil P (ppm)",y="Partial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))


ggplot()+ggtitle("Soil Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))
plot_grid(get_title(ggplot()+ggtitle("Soil Fungi")+theme(plot.title = element_text(hjust = 0.5,size = 36))),
          plot_grid(Lux_fung_soil_p_ppm_gdm,Lux_fung_soil_ca_gdm,rel_widths = c(1,1.3),
                    labels = c("a)","b)"),label_size = 36,label_x = c(0,-0.075)),nrow = 2,
          rel_heights = c(0.075,1))

#Graphing the combined GDM graphs####

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


All_LUX_sample_GDM_color=All_sample_GDM_color[names(All_sample_GDM_color)%in%
                                                unique(c(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long$metadata))]

All_LUX_sample_GDM_names=All_sample_GDM_names[names(All_sample_GDM_names)%in%
                                                unique(c(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_long$metadata,
                                                         GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long$metadata))]



(Lux_bact_root_sub_GDM=ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_LUX_sample_GDM_color,name=NULL,labels=All_LUX_sample_GDM_names)+
    theme_cowplot(font_size = 24)+labs(x=NULL,y="Root\nPartial ecological distance\n(Bray-Curtis)")+ggtitle("Bacteria")+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n15.1%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 36)))


(Lux_bact_soil_MR_sub_GDM=ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_gdm_long,
                                 aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_LUX_sample_GDM_color,name=NULL,labels=All_LUX_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n25.4%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y="Soil\nPartial ecological distance\n(Bray-Curtis)")+
    theme(legend.position = "none"))


(Lux_fung_root_sub_GDM=ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_gdm_long,
                              aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_LUX_sample_GDM_color,name=NULL,labels=All_LUX_sample_GDM_names)+
    theme_cowplot(font_size = 24)+labs(x=NULL,y=NULL)+ggtitle("Fungi")+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n15.5%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 36)))

(Lux_fung_soil_MR_sub_GDM=ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long,
                                 aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_LUX_sample_GDM_color,name=NULL,labels=All_LUX_sample_GDM_names)+
    annotation_custom(grobTree(textGrob(paste("p < 0.001\n26.0%"), x=0.1,  y=0.90, hjust=0,
                                        gp=gpar(col="black", fontsize=26))))+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.position = "none"))


All_LUX_sample_GDM_order=c("Geographic","collectionDate_N","ph","p_ppm","k_ppm","ca_ppm",
                           "percent_soil_moisture_dry_weight","soil_temp_1_avg_mean","dry_matter_yield_mg_ha_mean")

(Lux_fung_soil_MR_sub_GDM_legend=ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long,
                                        aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_LUX_sample_GDM_color,name=NULL,labels=All_LUX_sample_GDM_names,limits=All_LUX_sample_GDM_order)+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.text.align = 0,legend.key.size = unit(2, 'cm')))


plot_grid(Lux_bact_root_sub_GDM,Lux_fung_root_sub_GDM,
          Lux_bact_soil_MR_sub_GDM,Lux_fung_soil_MR_sub_GDM,
          ncol = 2, rel_widths = c(1,1),
          align = "v")


(four_panel_GDM_LUX=plot_grid(plot_grid(Lux_bact_root_sub_GDM,Lux_fung_root_sub_GDM,
                                        Lux_bact_soil_MR_sub_GDM,Lux_fung_soil_MR_sub_GDM,labels = c("a)","b)","c)","d)"),
                                        ncol = 2, rel_widths = c(1,1), label_size = 32,
                                        label_x = c(0,0.08,0,0.08),label_y = c(1,1,1.05,1.05),
                                        align = "v"),get_legend(Lux_fung_soil_MR_sub_GDM_legend), rel_widths = c(1,0.3)))




#Lux Arbor: Graphing the Euler diagrams od partitioning of deviance explained####



GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART$experiment=rep("Lux Arbor")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART$community=rep("Fungi")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART$niche=rep("Soil")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART[1,]
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2[1,]=
  c("soil alone",GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2[1,2:5])


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART$experiment=rep("Lux Arbor")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART$community=rep("Bacteria")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART$niche=rep("Soil")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART[9:15,]
colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2)=
  c("variableSet", "deviance", "experiment", "community", "niche")


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART$experiment=rep("Lux Arbor")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART$community=rep("Fungi")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART$niche=rep("Root")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART2=
  GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART[9:15,]
colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART2)=
  c("variableSet", "deviance", "experiment", "community", "niche")


GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART$experiment=rep("Lux Arbor")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART$community=rep("Bacteria")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART$niche=rep("Root")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART2=
  GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART[5:7,]


GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART=
  rbind(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2,
        GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_MR_sub_TabForm_PART2,
        GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_sub_TabForm_PART2,
        GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_sub_TabForm_PART2)




unique(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART$variableSet)
length(unique(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART$variableSet))
variableSet_order_f=c("soil alone",
                      "MET alone",
                      "soil_plant alone",
                      "geo alone",
                      "MET intersect soil_plant, exclude geo" ,
                      "MET intersect geo, exclude soil_plant",
                      "geo intersect soil_plant, exclude MET",
                      "soil_plant intersect MET intersect geo",
                      "MET intersect soil, exclude geo",
                      "MET intersect geo, exclude soil",
                      "geo intersect soil, exclude MET",
                      "soil intersect MET intersect geo",
                      "Time alone","Time intersect soil"   
)

variableSet_euler_names=c("soil alone"="Soil",
                          "MET alone"="MET",
                          "soil_plant alone"="Soil-Plant",
                          "geo alone"="Spatial",
                          "MET intersect soil_plant, exclude geo"="Soil-Plant&MET",
                          "MET intersect geo, exclude soil_plant"="Spatial&MET",
                          "geo intersect soil_plant, exclude MET"="Spatial&Soil-Plant",
                          "soil_plant intersect MET intersect geo"="Spatial&Soil-Plant&MET",
                          "MET intersect soil, exclude geo"="Soil&MET",
                          "MET intersect geo, exclude soil"="Spatial&MET",
                          "geo intersect soil, exclude MET"="Spatial&Soil",
                          "soil intersect MET intersect geo"="Spatial&Soil&MET",
                          "Time alone"="Temporal",
                          "Time intersect soil"="Temporal&Soil")


length(variableSet_euler_names)

GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART$variableSet_euler=
  GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART$variableSet


GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART <- GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART %>% 
  mutate(variableSet_euler=mapvalues(variableSet_euler,
                                     from = variableSet_order_f,
                                     to = variableSet_euler_names))


MET_SOIL_SPAT_alone_color=c("MET"="#313695","Soil"="#A50026","Spatial"="darkgrey")
MET_SOILPLANT_SPAT_alone_color=c("MET"="#313695","Soil-Plant"="#53342F","Spatial"="darkgrey")
PLANT_SOIL_SPAT_alone_color=c("Plant"="#006837","Soil"="#A50026","Spatial"="darkgrey")
TEMPO_SPAT_alone_color=c("Temporal"="lightgrey","Spatial"="darkgrey")
TEMPO_SOIL_alone_color=c("Temporal"="lightgrey","Soil"="#A50026")
variableSet_euler_alone_color=c("Spatial"="darkgrey",
                                "Soil"="#A50026",
                                "MET"="#313695",
                                "Plant"="#006837")

LUX_Root_bact_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Root")$deviance), 
           subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Root")$variableSet_euler)

LUX_Root_bact_euler_fit=euler(LUX_Root_bact_euler_list)



error_plot(LUX_Root_bact_euler_fit)




(LUX_Root_bact_euler_p=plot(LUX_Root_bact_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = TEMPO_SOIL_alone_color))







LUX_Soil_bact_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Soil")$deviance), 
           subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Bacteria" &niche=="Soil")$variableSet_euler)

LUX_Soil_bact_euler_fit=euler(LUX_Soil_bact_euler_list)
error_plot(LUX_Soil_bact_euler_fit)

(LUX_Soil_bact_euler_p=plot(LUX_Soil_bact_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = MET_SOILPLANT_SPAT_alone_color))

plot_grid(LUX_Root_bact_euler_p,LUX_Soil_bact_euler_p,ncol = 1)


LUX_Root_fung_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Fungi" &niche=="Root")$deviance), 
           subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Fungi" &niche=="Root")$variableSet_euler)

LUX_Root_fung_euler_fit=euler(LUX_Root_fung_euler_list)



error_plot(LUX_Root_fung_euler_fit)




(LUX_Root_fung_euler_p=plot(LUX_Root_fung_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = MET_SOIL_SPAT_alone_color))







LUX_Soil_fung_euler_list <- 
  setNames(as.numeric(subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Fungi" &niche=="Soil")$deviance), 
           subset(GLBRC018_OTU_MMPRNT_LUX_G5_sub_TabForm_PART,community=="Fungi" &niche=="Soil")$variableSet_euler)

LUX_Soil_fung_euler_fit=euler(LUX_Soil_fung_euler_list)

error_plot(LUX_Soil_fung_euler_fit)

(LUX_Soil_fung_euler_p=plot(LUX_Soil_fung_euler_fit,quantities = list(fontsize = 20),
                            labels = list(fontsize = 24,font=4),
                            fill = "#A50026"))


(Lux_fung_soil_MR_sub_GDM_legend2=ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR_sub_gdm_long,
                                         aes(x=x_predictors,y=y_values,color=metadata))+geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_color_manual(values = All_LUX_sample_GDM_color,name=NULL,labels=All_LUX_sample_GDM_names,limits=All_LUX_sample_GDM_order)+
    theme_cowplot(font_size = 24)+labs(x="Scaled predictors",y=NULL)+
    theme(legend.position = "bottom",legend.key.size = unit(2, 'cm')))


(six_panel_Euler_GDM_LUX=plot_grid(plot_grid(plot_grid(LUX_Root_bact_euler_p,LUX_Soil_bact_euler_p,ncol = 1),
                                             plot_grid(Lux_bact_root_sub_GDM,Lux_fung_root_sub_GDM,
                                                       Lux_bact_soil_MR_sub_GDM,Lux_fung_soil_MR_sub_GDM,
                                                       labels = c("a)","b)","c)","d)"),
                                                       ncol = 2, label_size = 40,
                                                       label_x = c(0,0.08,0,0.08),label_y = c(1.03,1.03,1.05,1.05),
                                                       align = "v"),plot_grid(LUX_Root_fung_euler_p,LUX_Soil_fung_euler_p,ncol = 1),
                                             rel_widths = c(0.2,1,0.2),ncol = 3),
                                   plot_grid(ggplot+geom_blank(),get_legend(Lux_fung_soil_MR_sub_GDM_legend2),rel_widths = c(0.5,1)), rel_heights = c(1,0.2),ncol = 1))


ggsave(six_panel_Euler_GDM_LUX,
       filename = "LUX_Euler_GDM_Percentage_deviance_plot.png",path = here::here("Manuscript","GDM_metadata_figs"),width = 28,height =14)
ggsave(six_panel_Euler_GDM_LUX,
       filename = "LUX_Euler_GDM_Percentage_deviance_plot.svg",path = here::here("Manuscript","GDM_metadata_figs"),width = 28,height =14)


