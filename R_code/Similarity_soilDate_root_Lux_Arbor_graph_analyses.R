## ---------------------------
##
## Script name: Temporal charaterization og the similarity between root and soil communities 
## across the Lux Arbor Growing Season
##
## Purpose of script: Modelling and graphing of the similarity of roots and soil communities of 
## switchgrass (Panicum virgatum) to determine the timing of peak similarity
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
library(otuSummary)

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



#LUX: Bacteria Pairwise distance####



#Let's graph the pairwise distance between communities

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,method = "bray")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M <- matrixConvert(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis, 
                                                           colname = c("sample1", "sample2", "bray"))#total distance 
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M,
                                                  GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate","plotRep","FertStatus",
                                                                                             "UTM_Lat_Cord", "UTM_Lon_Cord","sampleID_long")],
                                                  by.x = "sample1",by.y = "row.names")
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2)
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2)[4:10]=
  c("s1_Root_soil","s1_collectionDate","s1_plotRep","s1_FertStatus","s1_UTM_Lat_Cord", "s1_UTM_Lon_Cord","s1_samp_MMPRNT")

nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2)
#122265
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2,
                                                     GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate","plotRep","FertStatus",
                                                                                                "UTM_Lat_Cord", "UTM_Lon_Cord","sampleID_long")],
                                                     by.x = "sample2",by.y = "row.names")

colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)[11:17]=
  c("s2_Root_soil","s2_collectionDate","s2_plotRep","s2_FertStatus","s2_UTM_Lat_Cord", "s2_UTM_Lon_Cord","s2_samp_MMPRNT")
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)
#122265

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s1_s2_Root_soil=with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta, 
                                                                    interaction(s1_Root_soil,s2_Root_soil))

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s1_s2_collectionDate=with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta, 
                                                                         interaction(s1_collectionDate,s2_collectionDate))

#Pairwise distance calculation

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$distance=
  with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta, sqrt((s1_UTM_Lat_Cord-s2_UTM_Lat_Cord)^2+(s1_UTM_Lon_Cord-s2_UTM_Lon_Cord)^2))
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$similarity=
  with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta,1-bray)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta,
                                                         s1_Root_soil!=s2_Root_soil&s1_plotRep==s2_plotRep&s1_FertStatus==s2_FertStatus)



unique(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS$s1_s2_Root_soil)
#Root.Soil
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS)
#6265




#Previous sampling date 
unique(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS$s1_collectionDate)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS,
                                                              s1_s2_collectionDate=="5/29/2018.5/15/2018"|
                                                                s1_s2_collectionDate=="9/17/2018.9/4/2018" |
                                                                s1_s2_collectionDate=="8/20/2018.8/8/2018"|
                                                                s1_s2_collectionDate=="7/30/2018.7/9/2018"|
                                                                s1_s2_collectionDate=="6/25/2018.6/11/2018"|
                                                                s1_s2_collectionDate=="10/3/2018.9/17/2018")
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev)
#420
unique(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev$s1_s2_collectionDate)

#Same sample pairwise distance between organs
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta,
                                                              s2_samp_MMPRNT==s1_samp_MMPRNT)

#Combine sampling on same day and previous day
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev$comp_date=rep("prev")
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev)
dim(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev)
#420  22
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub$comp_date=rep("same")
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)
dim(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)
#138  22
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub=rbind(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev,
                                                                 GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)





similarity_date_dist_bact_mod=lmer(similarity~s1_collectionDate*s1_FertStatus*comp_date+(1|s1_plotRep)+(1|s1_plotRep:s1_samp_MMPRNT),
                                   data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub)
plot(similarity_date_dist_bact_mod)
hist(resid(similarity_date_dist_bact_mod))
qqPlot(resid(similarity_date_dist_bact_mod))
shapiro.test(resid(similarity_date_dist_bact_mod))
#W = 0.99486, p-value = 0.05908
simulateResiduals(fittedModel = similarity_date_dist_bact_mod, plot = T)



anova(similarity_date_dist_bact_mod)
#s1_collectionDate                         0.0199445 0.0039889     5 140.85  9.1209 1.552e-07 ***
#s1_FertStatus                             0.0000968 0.0000968     1 140.80  0.2213    0.6388    
#comp_date                                 0.0110175 0.0110175     1 405.43 25.1924 7.784e-07 ***
#s1_collectionDate:s1_FertStatus           0.0016365 0.0003273     5 140.83  0.7484    0.5886    
#s1_collectionDate:comp_date               0.0119142 0.0023828     5 405.41  5.4485 7.255e-05 ***
#s1_FertStatus:comp_date                   0.0000034 0.0000034     1 405.44  0.0078    0.9295    
#s1_collectionDate:s1_FertStatus:comp_date 0.0017232 0.0003446     5 405.41  0.7881    0.5587   


emmip(similarity_date_dist_bact_mod,  comp_date~ s1_collectionDate| s1_FertStatus)

emmeans(similarity_date_dist_bact_mod,pairwise~comp_date|s1_collectionDate,type = "response", adjust="fdr")


#$contrasts
#s1_collectionDate = 10/3/2018:
#contrast    estimate      SE  df t.ratio p.value
#prev - same  0.01295 0.00514 406   2.517  0.0122

#s1_collectionDate = 5/29/2018:
#  contrast    estimate      SE  df t.ratio p.value
#prev - same  0.02501 0.00504 406   4.967  <.0001

#s1_collectionDate = 6/25/2018:
#  contrast    estimate      SE  df t.ratio p.value
#prev - same  0.01410 0.00493 404   2.860  0.0045

#s1_collectionDate = 7/30/2018:
# contrast    estimate      SE  df t.ratio p.value
#prev - same  0.00748 0.00504 404   1.484  0.1385

#s1_collectionDate = 8/20/2018:
#  contrast    estimate      SE  df t.ratio p.value
#prev - same  0.01299 0.00514 407   2.527  0.0119

#s1_collectionDate = 9/17/2018:
#  contrast    estimate      SE  df t.ratio p.value
#prev - same -0.01046 0.00500 405  -2.094  0.0369

#Results are averaged over the levels of: s1_FertStatus 
#Degrees-of-freedom method: kenward-roger 






#####LUX: Fungi Pairwise distance####

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5, method = "bray")
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M <- matrixConvert(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis, 
                                                           colname = c("sample1", "sample2", "bray"))#total distance 
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M)
#119805
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M,
                                                  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate","plotRep","FertStatus",
                                                                                             "UTM_Lat_Cord", "UTM_Lon_Cord","sampleID_long")],
                                                  by.x = "sample1",by.y = "row.names")
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2)
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2)[4:10]=
  c("s1_Root_soil","s1_collectionDate","s1_plotRep","s1_FertStatus","s1_UTM_Lat_Cord", "s1_UTM_Lon_Cord","s1_samp_MMPRNT")

nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2)
#119805
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2,
                                                     GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate","plotRep","FertStatus",
                                                                                                "UTM_Lat_Cord", "UTM_Lon_Cord","sampleID_long")],
                                                     by.x = "sample2",by.y = "row.names")

colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)[11:17]=
  c("s2_Root_soil","s2_collectionDate","s2_plotRep","s2_FertStatus","s2_UTM_Lat_Cord", "s2_UTM_Lon_Cord","s2_samp_MMPRNT")
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)
#119805

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s1_s2_Root_soil=with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta, 
                                                                    interaction(s1_Root_soil,s2_Root_soil))

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s1_s2_collectionDate=with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta, 
                                                                         interaction(s1_collectionDate,s2_collectionDate))

#Pairwise distance calculation

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$distance=
  with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta, sqrt((s1_UTM_Lat_Cord-s2_UTM_Lat_Cord)^2+(s1_UTM_Lon_Cord-s2_UTM_Lon_Cord)^2))

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$similarity=
  with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta,1-bray)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta,
                                                         s1_Root_soil!=s2_Root_soil&s1_plotRep==s2_plotRep&s1_FertStatus==s2_FertStatus)



unique(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS$s1_s2_Root_soil)
#Soil.Root
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS)
#6152


#Previous sampling date 
unique(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS$s1_collectionDate)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS,
                                                              s1_s2_collectionDate=="5/15/2018.5/29/2018"|
                                                                s1_s2_collectionDate=="9/4/2018.9/17/2018" |
                                                                s1_s2_collectionDate=="8/8/2018.8/20/2018"|
                                                                s1_s2_collectionDate=="7/9/2018.7/30/2018"|
                                                                s1_s2_collectionDate=="6/11/2018.6/25/2018"|
                                                                s1_s2_collectionDate=="9/17/2018.10/3/2018")
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev)
#411


#Same sample pairwise distance between organs
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta,
                                                              s2_samp_MMPRNT==s1_samp_MMPRNT)
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)
#137
unique(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub$s1_s2_Root_soil)
#Soil.Root
unique(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub$s2_collectionDate)




GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev$comp_date=rep("prev")
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev)
dim(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev)
#411  22
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub$comp_date=rep("same")
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)
dim(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)
#137  22
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub=rbind(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev,
                                                                 GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_org_sub)





similarity_date_dist_fung_mod=lmer(logit(similarity)~s2_collectionDate*s2_FertStatus*comp_date+(1|s2_plotRep)+(1|s2_plotRep:s2_samp_MMPRNT),
                                   data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub)
plot(similarity_date_dist_fung_mod)
hist(resid(similarity_date_dist_fung_mod))
qqPlot(resid(similarity_date_dist_fung_mod))
shapiro.test(resid(similarity_date_dist_fung_mod))
#W = 0.99581, p-value = 0.1514
simulateResiduals(fittedModel = similarity_date_dist_fung_mod, plot = T)



anova(similarity_date_dist_fung_mod)
#s2_collectionDate                         3.8714  0.7743     5 140.12 10.7359 9.450e-09 ***
#s2_FertStatus                             0.1262  0.1262     1 140.15  1.7504    0.1880    
#comp_date                                 3.4914  3.4914     1 396.12 48.4102 1.439e-11 ***
#s2_collectionDate:s2_FertStatus           0.2618  0.0524     5 140.11  0.7260    0.6050    
#s2_collectionDate:comp_date               0.5290  0.1058     5 396.08  1.4671    0.1996    
#s2_FertStatus:comp_date                   0.0745  0.0745     1 396.07  1.0333    0.3100    
#s2_collectionDate:s2_FertStatus:comp_date 0.2639  0.0528     5 396.10  0.7319    0.5998  



emmip(similarity_date_dist_fung_mod,  comp_date~ s2_collectionDate| s2_FertStatus)

emmeans(similarity_date_dist_fung_mod,pairwise~comp_date|s2_collectionDate,type = "response", adjust="fdr")



#Similarity Plots####


comp_date_names=c("prev"="2-week", "same"="Same day")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean=
  GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub|>group_by(s1_collectionDate,comp_date)|>
  summarise(response=mean(similarity),SE=sd(similarity)/sqrt(n()))

pairwise_comp_bact_sim_lux=data.frame(subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean,
                                             s1_collectionDate!="7/30/2018")|>group_by(s1_collectionDate)|>
                                        summarise(mean_max=max(response+SE)))|>
  mutate(s1_collectionDate=as.Date(s1_collectionDate,format="%m/%d/%Y"))|>mutate(x_min=s1_collectionDate-5,x_max=s1_collectionDate+5)
pairwise_comp_bact_sim_lux$annot_text=c("*","***","**"," * ","  *  ")

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean|>ungroup()|>summarise(mean_max=max(response+SE), mean_min=min(response-SE))

(bact_pairwise_dist_soil_root_time=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean,
                                          aes(x=as.Date(s1_collectionDate,format="%m/%d/%Y"),y=response))+
    geom_point(aes(color=comp_date),size=10,position = position_dodge(15))+
    geom_errorbar(aes(color=comp_date,ymin=response-SE,ymax=response+SE),position = position_dodge(15))+
    geom_signif(data=pairwise_comp_bact_sim_lux,mapping = aes(y_position = mean_max+0.002, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+scale_y_continuous(name="Bacteria\nSimilarity between\nroots and soils",limits = c(0.18,0.26))+
    scale_fill_manual(values = c("black","grey"), name=NULL)+scale_color_manual(values = c("darkgrey","black"),label=comp_date_names, name=NULL)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+
    theme(axis.text.x = element_blank(),legend.position = "none"))






GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean=
  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub|>group_by(s2_collectionDate,comp_date)|>
  summarise(response=mean(similarity),SE=sd(similarity)/sqrt(n()))


pairwise_comp_fung_sim_lux=data.frame(subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean,
                                             s2_collectionDate!="9/17/2018")|>group_by(s2_collectionDate)|>
                                        summarise(mean_max=max(response+SE)))|>
  mutate(s2_collectionDate=as.Date(s2_collectionDate,format="%m/%d/%Y"))|>mutate(x_min=s2_collectionDate-5,x_max=s2_collectionDate+5)
pairwise_comp_fung_sim_lux$annot_text=c("**","***","#"," * "," *** ")


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean|>ungroup()|>summarise(mean_max=max(response+SE), mean_min=min(response-SE))


(fung_pairwise_dist_soil_root_time=ggplot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_RS_prev_sub_mean,
                                          aes(x=as.Date(s2_collectionDate,format="%m/%d/%Y"),y=response))+
    geom_point(aes(color=comp_date),size=10,position = position_dodge(15))+
    geom_errorbar(aes(color=comp_date,ymin=response-SE,ymax=response+SE),position = position_dodge(15))+
    geom_signif(data=pairwise_comp_fung_sim_lux,mapping = aes(y_position = mean_max+0.002, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+scale_y_continuous(name = "Fungi\nSimilarity between\nroots and soils", limits = c(0.045,0.12))+
    scale_fill_manual(values = c("black","grey"), name=NULL)+scale_color_manual(values = c("darkgrey","black"),label=comp_date_names, name=NULL)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+
    theme(legend.position = "none"))

LENGEND_pairwise_dist_soil_root_time=get_legend(
  fung_pairwise_dist_soil_root_time + theme(legend.position = "left")
)

plot_grid(plot_grid(bact_pairwise_dist_soil_root_time,fung_pairwise_dist_soil_root_time,
                    labels = c("a)","b)"),label_size = 36,rel_heights = c(1,1.1),ncol = 1),
          LENGEND_pairwise_dist_soil_root_time,rel_widths = c(2,0.5))

ggsave(plot_grid(plot_grid(bact_pairwise_dist_soil_root_time,fung_pairwise_dist_soil_root_time,
                           labels = c("a)","b)"),label_size = 36,rel_heights = c(1,1.1),ncol = 1),
                 LENGEND_pairwise_dist_soil_root_time,rel_widths = c(2,0.5)),
       filename = "LUX_pairwise_dist_soil_root_time_p.png",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 13,height =10)

ggsave(plot_grid(plot_grid(bact_pairwise_dist_soil_root_time,fung_pairwise_dist_soil_root_time,
                           labels = c("a)","b)"),label_size = 36,rel_heights = c(1,1.1),ncol = 1),
                 LENGEND_pairwise_dist_soil_root_time,rel_widths = c(2,0.5)),
       filename = "LUX_pairwise_dist_soil_root_time_p.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 13,height =10)

