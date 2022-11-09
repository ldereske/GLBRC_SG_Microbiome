
## ---------------------------
##
## Script name: Community characterization of the switchgrass microbiomes across 
## the Lux Arbor Growing Season
##
## Purpose of script: Graphing and characterization of the bacteria and fungi
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
library(gdm)
library(ggnewscale)
library(otuSummary)
"%w/o%" <- function(x,y)!('%in%'(x,y))
Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "English")

#Load code for calculating confidence intervals for GDMs 
source(here::here("R_functions","dataUncertainty.R"))
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

#Let's ordinate
set.seed(2021)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#20: scale factor of the gradient < sfgrmin
#0.06852062  


plot_ordination(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_ord)+
  geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=Root_soil),size=4)+geom_text(aes(label=SampleID))+
  theme_bw()



#Amp129 is an outlier so I am going to remove it 

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,SampleID!="Amp129")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
#495
rm(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)

#Let's ordinate
set.seed(2021)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#20: scale factor of the gradient < sfgrmin
#0.06852062  


plot_ordination(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_ord)+
  geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=Root_soil),size=4)+geom_text(aes(label=sampleID_seq))+
  theme_bw()


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

summary(as.numeric(as.Date(sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)$collectionDate, format="%m/%d/%Y")))




#Let's generate NMDS
set.seed(2022)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root,method = "NMDS")
#*** Solution reached
#0.2373194     

set.seed(2022)
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#1: no. of iterations >= maxit
#19: stress ratio > sratmax
#0.2194553    


set.seed(2022)
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#1: no. of iterations >= maxit
#18: stress ratio > sratmax
#1: scale factor of the gradient < sfgrmin
#0.237239  

set.seed(2022)
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#1: no. of iterations >= maxit
#16: stress ratio > sratmax
#3: scale factor of the gradient < sfgrmin
#0.1960742     



GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_NMDS_points=merge(sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root),GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_ord$points,
                                                       by="row.names")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_NMDS_points_sum=GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_NMDS_points%>%group_by(collectionDate,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_NMDS_points=merge(sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root),GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_ord$points,
                                                       by="row.names")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_NMDS_points_sum=GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_NMDS_points%>%group_by(collectionDate,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_NMDS_points=merge(sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil),GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_ord$points,
                                                       by="row.names")
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_NMDS_points_sum=GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_NMDS_points%>%group_by(collectionDate,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_NMDS_points=merge(sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil),GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_ord$points,
                                                       by="row.names")
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_NMDS_points_sum=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_NMDS_points%>%group_by(collectionDate,FertStatus,plotRep)%>%
  summarise(MDS1_mean=mean(MDS1),MDS1_se=sd(MDS1)/sqrt(n()),MDS2_mean=mean(MDS2),MDS2_se=sd(MDS2)/sqrt(n()))


#####GDM of Space and Time####


#Bacteria
#Roots 





GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data=sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)[,c("SampleID","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)
#142
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)[colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_metadata=data.frame(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root,method = "bray")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF),GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site[1:10,1:10]

summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                  predData=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_metadata)


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:  12.901



date_varSet <- vector("list", 1)
names(date_varSet)=c("date")
date_varSet$date=c("collectionDate_N")
summary(date_varSet)
#Variance partitioning
gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm,varSets = date_varSet,partSpace= TRUE)
#        VARIABLE_SET DEVIANCE
#1                geo     1.31
#2               date    11.46
#3         geo & date    12.90
#4         date alone    11.59
#5          geo alone     1.45
#6        UNEXPLAINED    87.10



gdm.varImp_MOD(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm, geo = T, nPerm = 100)
#                           All predictors
#Model deviance                    130.389
#Percent deviance explained         12.901
#Model p-value                       0.000
#Fitted permutations                95.000

#$`Predictor p-values`
#All predictors
#Geographic                    0
#collectionDate_N              0

#Root uncertainy 

root_bact_comb_splin_uncert=dataUncertainty(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm,0.7,100, geo = T)
root_bact_comb_splin_uncert[,"Root_soil"]="Root"


#Soils

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data=sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)[,c("SampleID","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)
#353
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)[colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_metadata=data.frame(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil,method = "bray")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF),GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF_site[1:10,1:10]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray_DF_site, 3, siteColumn="site", XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                  predData=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_metadata)


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:  15.54

gdm.varImp_MOD(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_TabForm, geo = T, nPerm = 100)
#                           All predictors
#Model deviance                    447.246
#Percent deviance explained         15.540
#Model p-value                       0.000

date_varSet <- vector("list", 1)
names(date_varSet)=c("date")
date_varSet$date=c("collectionDate_N")
summary(date_varSet)
#Variance partitioning
gdm.partition.deviance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_TabForm,varSets = date_varSet,partSpace= TRUE)
#        VARIABLE_SET DEVIANCE
#1                geo    15.08
#2               date     0.43
#3         geo & date    15.54
#4         date alone     0.46
#5          geo alone    15.11
#6        UNEXPLAINED    84.46
#7 ALL VARIABLES (NA)       NA






#Soil uncertainly 

soil_bact_comb_splin_uncert=dataUncertainty(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_TabForm,0.7,100, geo = T)
soil_bact_comb_splin_uncert[,"Root_soil"]="Soil"




root_soil_bact_comb_splin=rbind(root_bact_comb_splin,soil_bact_comb_splin)



colnames(root_bact_comb_splin_uncert)
colnames(soil_bact_comb_splin_uncert)
root_soil_bact_comb_splin_uncert=rbind(root_bact_comb_splin_uncert,soil_bact_comb_splin_uncert)



#Fungi


#Roots
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data=sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)[,c("sampleID_seq","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)
#141
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)[colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_metadata=data.frame(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root,method = "bray")
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF),GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site[1:10,1:10]
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                  predData=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_metadata)


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:  7.206

gdm.varImp_MOD(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_TabForm, geo = T, nPerm = 100)
#                           All predictors
#Model deviance                    582.902
#Percent deviance explained          7.206
#Model p-value                       0.000
#Fitted permutations                92.000

date_varSet <- vector("list", 1)
names(date_varSet)=c("date")
date_varSet$date=c("collectionDate_N")
summary(date_varSet)
#Variance partitioning
gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_TabForm,varSets = date_varSet,partSpace= TRUE)
#        VARIABLE_SET DEVIANCE
#1                geo     4.09
#2               date     3.00
#3         geo & date     7.21
#4         date alone     3.12
#5          geo alone     4.21
#6        UNEXPLAINED    92.79
#7 ALL VARIABLES (NA)       NA





root_fung_comb_splin_uncert=dataUncertainty(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_TabForm,0.7,100, geo = T)
root_fung_comb_splin_uncert[,"Root_soil"]="Root"

#Soils


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data=sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil)[,c("sampleID_seq","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)
#349
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)[colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_metadata=data.frame(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil,method = "bray")
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF=data.frame(as.matrix(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray))
row.names(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF),GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF_site[1:10,1:10]
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF_site)
summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_TabForm<- formatsitepair(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray_DF_site, 3, siteColumn="site", XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                                  predData=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_metadata)


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod<- gdm(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_TabForm, geo=TRUE)
summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod)
plot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:  15.06

gdm.varImp_MOD(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_TabForm, geo = T, nPerm = 100)

date_varSet <- vector("list", 1)
names(date_varSet)=c("date")
date_varSet$date=c("collectionDate_N")
summary(date_varSet)
#Variance partitioning
gdm.partition.deviance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_TabForm,varSets = date_varSet,partSpace= TRUE)
#        VARIABLE_SET DEVIANCE
#1                geo    14.22
#2               date     0.79
#3         geo & date    15.06
#4         date alone     0.84
#5          geo alone    14.27
#6        UNEXPLAINED    84.94
#7 ALL VARIABLES (NA)       NA



#Soil uncertainty

soil_fung_comb_splin_uncert=dataUncertainty(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_TabForm,0.7,100, geo = T)
soil_fung_comb_splin_uncert[,"Root_soil"]="Soil"


root_soil_fung_comb_splin=rbind(root_fung_comb_splin,soil_fung_comb_splin)




colnames(root_fung_comb_splin_uncert)
colnames(soil_fung_comb_splin_uncert)
root_soil_fung_comb_splin_uncert=rbind(root_fung_comb_splin_uncert,soil_fung_comb_splin_uncert)


#Mark the actual sampling dates

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)



bacterial_root_sampling_date=data.frame(sampleDate=c(unique(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,Root_soil=="Root")$collectionDate)),
                                        Root_soil=rep("Root"))
bacterial_soil_sampling_date=data.frame(sampleDate=c(unique(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,Root_soil=="Soil")$collectionDate)),
                                        Root_soil=rep("Soil"))

sampling_date=rbind(bacterial_root_sampling_date,
                    bacterial_soil_sampling_date[bacterial_soil_sampling_date$sampleDate%w/o%bacterial_root_sampling_date$sampleDate,])



#Mark the actual sampling distances



GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis <- 
  matrixConvert((dist(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data[,c("UTM_Lat_Cord", "UTM_Lon_Cord")], method = "euclidean")), 
                colname = c("sample1", "sample2", "distance"))#total distance 

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2=merge(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis,
                                                  sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)[,c("plotRep")],
                                                  by.x = "sample1",by.y="row.names")
head(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2)
colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2)[colnames(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2)=="plotRep"]="s1_plotRep"


GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m=merge(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_dis_m2,
                                                 sample_data(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)[,c("plotRep")],
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

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis <- 
  matrixConvert((dist(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data[,c("UTM_Lat_Cord", "UTM_Lon_Cord")], method = "euclidean")), 
                colname = c("sample1", "sample2", "distance"))#total distance 

GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m2=merge(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis,
                                                  sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)[,c("plotRep")],
                                                  by.x = "sample1",by.y="row.names")
head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m2)
colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m2)[colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m2)=="plotRep"]="s1_plotRep"


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m=merge(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m2,
                                                 sample_data(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_MR)[,c("plotRep")],
                                                 by.x = "sample2",by.y="row.names")

head(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m)
colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m)[colnames(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m)=="plotRep"]="s2_plotRep"


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m$s1_s2_plotRep=with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m,
                                                              interaction(s1_plotRep,s2_plotRep))
unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m$s1_s2_plotRep)
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m$plot_comp=with(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m,
                                                          ifelse(s1_plotRep==s2_plotRep|s1_s2_plotRep=="4.3"|
                                                                   s1_s2_plotRep=="3.4","within",
                                                                 ifelse(s1_s2_plotRep=="1.4"|s1_s2_plotRep=="4.1"|
                                                                          s1_s2_plotRep=="1.3"|s1_s2_plotRep=="3.1",
                                                                        "long","short")))
unique(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m$plot_comp)

fungal_sampling_dist=GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_dis_m%>%group_by(plot_comp)%>%
  summarise(max_dis=max(distance),min_dis=min(distance))






######NMDS and GDM graph: Compartment colored NMDS####

Date_F <- function(x){
  format(as.Date(x, origin = '1970-01-01',"%Y-%m-%d"),"%b")
}

(LUX_compart_bact_root_fill_color_sum_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                         y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    stat_ellipse(aes(group=plotRep, fill=factor(plotRep)),geom = "polygon",alpha=0.5)+
    scale_fill_manual(values = c("#DEDDDD","#C2C1C1","#A6A6A6","#8B8A8A"),name="Plot rep")+new_scale("fill")+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=1)+theme_bw(base_size = 24)+
    geom_errorbar(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),width=0.01)+geom_errorbarh(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),height=0.01)+
    xlab("NMDS1")+ylab("NMDS2")+
    annotate("text", x = -0.363, y = -0.4, label = "Stress = 0.237", size=6)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_shape_manual(values=c(19,24),name="N fert")+scale_fill_manual(values = c("black","white"),name="N fert")+
    theme(panel.grid = element_blank(),legend.position = "none"))


(LUX_compart_fung_root_fill_color_sum_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                         y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    stat_ellipse(aes(group=plotRep, fill=factor(plotRep)),geom = "polygon",alpha=0.5)+
    scale_fill_manual(values = c("#DEDDDD","#C2C1C1","#A6A6A6","#8B8A8A"),name="Plot rep")+new_scale("fill")+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=1)+theme_bw(base_size = 24)+
    geom_errorbar(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),width=0.01)+geom_errorbarh(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),height=0.01)+
    xlab("NMDS1")+ylab("NMDS2")+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_shape_manual(values=c(19,24),name="N fert")+scale_fill_manual(values = c("black","white"),name="N fert")+
    annotate("text", x = -0.333, y = -0.615, label = "Stress = 0.219", size=6)+
    theme(panel.grid = element_blank(),legend.position = "none"))

(LUX_compart_bact_soil_fill_color_sum_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                         y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    stat_ellipse(aes(group=plotRep, fill=factor(plotRep)),geom = "polygon",alpha=0.5)+
    scale_fill_manual(values = c("#DEDDDD","#C2C1C1","#A6A6A6","#8B8A8A"),name="Plot rep")+new_scale("fill")+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=1)+theme_bw(base_size = 24)+
    xlab("NMDS1")+ylab("NMDS2")+
    geom_errorbar(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),width=0.01)+geom_errorbarh(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),height=0.01)+
    annotate("text", x = -0.255, y = -0.317, label = "Stress = 0.237", size=6)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+scale_shape_manual(values=c(19,24),name="N fert")+scale_fill_manual(values = c("black","white"),name="N fert")+
    theme(panel.grid = element_blank(),legend.position = "none",axis.title.y = element_blank()))


(LUX_compart_fung_soil_fill_color_sum_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                         y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    stat_ellipse(aes(group=plotRep, fill=factor(plotRep)),geom = "polygon",alpha=0.5)+
    scale_fill_manual(values = c("#DEDDDD","#C2C1C1","#A6A6A6","#8B8A8A"),name="Plot rep")+new_scale("fill")+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=1)+theme_bw(base_size = 24)+
    xlab("NMDS1")+ylab("NMDS2")+
    geom_errorbar(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),width=0.01)+
    geom_errorbarh(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),height=0.01)+
    annotate("text", x = -0.51, y = -0.776, label = "Stress = 0.196", size=6)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+scale_shape_manual(values=c(19,24),name="N fert")+scale_fill_manual(values = c("black","white"),name="N fert")+
    theme(panel.grid = element_blank(),legend.position = "none",axis.title.y = element_blank()))


(Time_root_soil_bact_GDM_uncert_marks_color_p <- ggplot(subset(root_soil_bact_comb_splin_uncert, factor=="collectionDate_N"),aes(x=as.Date(fullPlotX,origin="1970-01-01"),y=fullPlotY,color=Root_soil)) +
    geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+
    scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)",limits = c(-0.005,0.174))+
    scale_colour_manual(values = c("#33A02C","#B15928"),name=NULL)+scale_x_date(name= NULL)+
    geom_smooth(se=FALSE,method = 'loess',size=2,aes(x=as.Date(highBoundX,origin="1970-01-01"),y=as.Date(highBoundY,origin="1970-01-01")),linetype="dashed")+
    geom_smooth(se=FALSE,method = 'loess',size=2,aes(x=as.Date(lowBoundX,origin="1970-01-01"),y=as.Date(lowBoundY,origin="1970-01-01")),linetype="dashed")+
    geom_segment(data = sampling_date,aes(x=as.Date(sampleDate,format = "%m/%d/%Y"),xend = as.Date(sampleDate,format = "%m/%d/%Y"),y=-0.005,yend=-Inf),size=2,color="black")+
    theme_bw(base_size=24)+theme(legend.position = c(0.1,0.9),panel.grid = element_blank()))

(Time_root_soil_fung_GDM_uncert_marks_color_p <- ggplot(subset(root_soil_fung_comb_splin_uncert, factor=="collectionDate_N"),aes(x=as.Date(fullPlotX,origin="1970-01-01"),y=fullPlotY,color=Root_soil)) +
    geom_smooth(linetype="solid",se=FALSE,method = 'loess',size=2)+scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)",limits = c(-0.005,0.355))+
    geom_smooth(se=FALSE,method = 'loess',size=2,aes(x=as.Date(highBoundX,origin="1970-01-01"),y=as.Date(highBoundY,origin="1970-01-01")),linetype="dashed")+
    geom_smooth(se=FALSE,method = 'loess',size=2,aes(x=as.Date(lowBoundX,origin="1970-01-01"),y=as.Date(lowBoundY,origin="1970-01-01")),linetype="dashed")+
    scale_colour_manual(values = c("#33A02C","#B15928"),name=NULL)+scale_x_date(name= NULL)+
    geom_segment(data = sampling_date,aes(x=as.Date(sampleDate,format = "%m/%d/%Y"),xend = as.Date(sampleDate,format = "%m/%d/%Y"),y=-0.005,yend=-Inf),size=2,color="black")+
    theme_bw(base_size=24)+theme(legend.position = c(0.1,0.9),panel.grid = element_blank()))


(Dist_root_soil_bact_GDM_uncert_marks_color_p <- ggplot(subset(root_soil_bact_comb_splin_uncert, factor=="Geographic")) +
    geom_smooth(aes(x=fullPlotX,y=fullPlotY,color=Root_soil,linetype=Root_soil),se=FALSE,method = 'loess',size=2)+
    scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)",limits = c(-0.005,0.174))+
    scale_linetype_manual(values=c("solid","longdash"),name=NULL)+
    geom_smooth(aes(x=highBoundX,y=highBoundY,color=Root_soil),se=FALSE,method = 'loess',size=2,linetype="dashed")+
    geom_smooth(aes(x=lowBoundX,y=lowBoundY,color=Root_soil),se=FALSE,method = 'loess',size=2,linetype="dashed")+
    scale_colour_manual(values = c("#33A02C","#B15928"),name=NULL)+scale_x_continuous(name= "Distance (m)")+
    geom_rect(data = bacterial_sampling_dist,aes(xmin=min_dis,xmax = max_dis,ymin=-0.005,ymax=-Inf,group=plot_comp),
              size=2,color="black",fill="black")+
    theme_bw(base_size=24)+theme(legend.position = c(0.1,0.9),panel.grid = element_blank(),axis.title.y = element_blank()))

(Dist_root_soil_fung_GDM_uncert_marks_color_p <- ggplot(subset(root_soil_fung_comb_splin_uncert, factor=="Geographic")) +
    geom_smooth(aes(x=fullPlotX,y=fullPlotY,color=Root_soil,linetype=Root_soil),se=FALSE,method = 'loess',size=2)+
    scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)",limits = c(-0.005,0.355))+
    scale_linetype_manual(values=c("solid","longdash"),name=NULL)+
    geom_smooth(aes(x=highBoundX,y=highBoundY,color=Root_soil),se=FALSE,method = 'loess',size=2,linetype="dashed")+
    geom_smooth(aes(x=lowBoundX,y=lowBoundY,color=Root_soil),se=FALSE,method = 'loess',size=2,linetype="dashed")+
    scale_colour_manual(values = c("#33A02C","#B15928"),name=NULL, )+scale_x_continuous(name= "Distance (m)")+
    geom_rect(data = fungal_sampling_dist,aes(xmin=min_dis,xmax = max_dis,ymin=-0.005,ymax=-Inf,group=plot_comp),
              size=2,color="black",fill="black")+
    theme_bw(base_size=24)+theme(legend.position = c(0.1,0.9),panel.grid = element_blank(),axis.title.y = element_blank()))


(Lux_arbor_BACT_4_sum_GDM_color_panel=plot_grid(LUX_compart_bact_root_fill_color_sum_p2,
                                                LUX_compart_bact_soil_fill_color_sum_p2,
                                                Time_root_soil_bact_GDM_uncert_marks_color_p,
                                                Dist_root_soil_bact_GDM_uncert_marks_color_p,
                                                ncol = 2,align = "hv",
                                                label_size = 26,
                                                rel_heights = c(1,1)))


(Lux_arbor_FUNG_4_sum_GDM_color_panel=plot_grid(LUX_compart_fung_root_fill_color_sum_p2,
                                                LUX_compart_fung_soil_fill_color_sum_p2,
                                                Time_root_soil_fung_GDM_uncert_marks_color_p,
                                                Dist_root_soil_fung_GDM_uncert_marks_color_p,
                                                ncol = 2,align = "hv",
                                                label_size = 26,
                                                rel_heights = c(1,1)))

ggsave(Lux_arbor_BACT_4_sum_GDM_color_panel, filename = "Bact_NMDS_GDM_Lux_Arbor_sum_fert_compart_plot_color_p.png",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height = 15)


ggsave(Lux_arbor_BACT_4_sum_GDM_color_panel, filename = "Bact_NMDS_GDM_Lux_Arbor_sum_fert_compart_plot_color_p.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height = 15)


ggsave(Lux_arbor_FUNG_4_sum_GDM_color_panel, filename = "Fung_NMDS_GDM_Lux_Arbor_sum_fert_compart_plot_color_p.png",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height = 15)


ggsave(Lux_arbor_FUNG_4_sum_GDM_color_panel, filename = "Fung_NMDS_GDM_Lux_Arbor_sum_fert_compart_plot_color_p.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height = 15)

(LEGEND_LUX_compart_bact_root_fill_color_sum_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                                y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    stat_ellipse(aes(group=plotRep, fill=factor(plotRep)),geom = "polygon",alpha=0.5)+
    scale_fill_manual(values = c("#DEDDDD","#C2C1C1","#A6A6A6","#8B8A8A"),name="Plot rep")+new_scale("fill")+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 24)+
    geom_errorbar(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),width=0.01)+geom_errorbarh(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),height=0.01)+
    xlab("NMDS1")+ylab("NMDS2")+
    annotate("text", x = -0.363, y = -0.4, label = "Stress = 0.237", size=6)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_shape_manual(values=c(19,24),name="N fert")+scale_fill_manual(values = c("black","white"),name="N fert")+
    theme(panel.grid = element_blank(),axis.title.x = element_blank()))


(LEGEND_LUX_compart_bact_soil_fill_color_sum_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_NMDS_points_sum,aes(x=MDS1_mean,xmax=MDS1_mean+MDS1_se,xmin=MDS1_mean-MDS1_se,
                                                                                                                y=MDS2_mean,ymax=MDS2_mean+MDS2_se,ymin=MDS2_mean-MDS2_se))+
    stat_ellipse(aes(group=plotRep, fill=factor(plotRep)),geom = "polygon",alpha=0.5)+
    scale_fill_manual(values = c("#DEDDDD","#C2C1C1","#A6A6A6","#8B8A8A"),name="Plot rep")+new_scale("fill")+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw(base_size = 24)+
    xlab("NMDS1")+ylab("NMDS2")+
    geom_errorbar(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),width=0.01)+geom_errorbarh(aes(color=as.Date(collectionDate, format="%m/%d/%Y")),height=0.01)+
    annotate("text", x = -0.255, y = -0.317, label = "Stress = 0.238", size=6)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+scale_shape_manual(values=c(19,24),name="N fert")+scale_fill_manual(values = c("black","white"),name="N fert")+
    theme(panel.grid = element_blank()))



(LEGEND_Lux_arbor_GDM_color_panel=plot_grid(get_legend(LEGEND_LUX_compart_bact_root_fill_color_sum_p2),
                                            get_legend(LEGEND_LUX_compart_bact_soil_fill_color_sum_p2)))


ggsave(LEGEND_Lux_arbor_GDM_color_panel, filename = "LEGEND_NMDS_GDM_Lux_Arbor_sum_fert_compart_plot_color_p.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 17,height = 15)



