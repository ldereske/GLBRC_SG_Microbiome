## ---------------------------
##
## Script name: Community correlation between bacteria and fungi in the switchgrass microbiomes of 
## the Marginal Lands Experiment (MLE)
##
## Purpose of script: Analysis and graphing of the correspondence between bacteria and fungi
## in roots and soils of switchgrass (Panicum virgatum)
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
library(lmerTest)
library(emmeans)
library(DHARMa)
library(ggsignif)
library(RColorBrewer)

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


# Protestes: MATCH SAMPLES -------------------------------------------------------------------------------
# Match samples to the same ones for fungi and bacteria ---------------------------------------
#Code Selection Written by Gian Maria Niccolò Benucci
#Modified by Lukas Bell-Dereske
# MLE
# root



intersect(
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_root@sam_data$sampleID_long,
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_root@sam_data$sampleID_long
) -> MLE_good_root
length(MLE_good_root)

MLE_physeq_fungi_root <-
  subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE_G5_root,
                 sampleID_long %in% MLE_good_root)
MLE_physeq_fungi_root@otu_table <-
  MLE_physeq_fungi_root@otu_table[which(rowSums(MLE_physeq_fungi_root@otu_table) > 0), ]
MLE_physeq_fungi_root

sample_names(MLE_physeq_fungi_root) <-
  MLE_physeq_fungi_root@sam_data$sampleID_long
MLE_physeq_fungi_root@sam_data
head(MLE_physeq_fungi_root@otu_table)

MLE_physeq_bact_root <-
  subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE_G5_root,
                 sampleID_long %in% MLE_good_root)
MLE_physeq_bact_root@otu_table <-
  MLE_physeq_bact_root@otu_table[which(rowSums(MLE_physeq_bact_root@otu_table) > 0), ]
MLE_physeq_bact_root

sample_names(MLE_physeq_bact_root) <-
  MLE_physeq_bact_root@sam_data$sampleID_long
MLE_physeq_bact_root@sam_data

# soil



intersect(
  GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil@sam_data$sampleID_long,
  GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil@sam_data$sampleID_long
) -> MLE_good_soil
length(MLE_good_soil)


MLE_physeq_fungi_soil <-
  subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil,
                 sampleID_long %in% MLE_good_soil)
MLE_physeq_fungi_soil@otu_table <-
  MLE_physeq_fungi_soil@otu_table[which(rowSums(MLE_physeq_fungi_soil@otu_table) > 0), ]
MLE_physeq_fungi_soil

sample_names(MLE_physeq_fungi_soil) <-
  MLE_physeq_fungi_soil@sam_data$sampleID_long
MLE_physeq_fungi_soil@sam_data

MLE_physeq_bact_soil <-
  subset_samples(GLBRC018_OTU_bact_MMPRNT_MLE_G5_soil,
                 sampleID_long %in% MLE_good_soil)
MLE_physeq_bact_soil@otu_table <-
  MLE_physeq_bact_soil@otu_table[which(rowSums(MLE_physeq_bact_soil@otu_table) > 0),]
MLE_physeq_bact_soil

sample_names(MLE_physeq_bact_soil) <-
  MLE_physeq_bact_soil@sam_data$sampleID_long
MLE_physeq_bact_soil@sam_data




# ********************************-------------------------------------------------------------
# PROCRUSTES ROTATION ANALYSIS ----------------------------------------------------------------


# calculating NMDS ----------------------------------------------------------------------------
set.seed(2021)
MLE_nmds_its_root = ordinate(MLE_physeq_fungi_root, method ="NMDS", distance="bray")
MLE_nmds_its_root
stressplot(MLE_nmds_its_root)

MLE_nmds_16s_root = ordinate(MLE_physeq_bact_root, method ="NMDS", distance="bray")
MLE_nmds_16s_root
stressplot(MLE_nmds_16s_root)

MLE_nmds_its_soil = ordinate(MLE_physeq_fungi_soil, method ="NMDS", distance="bray")
MLE_nmds_its_soil
stressplot(MLE_nmds_its_soil)

MLE_nmds_16s_soil = ordinate(MLE_physeq_bact_soil, method ="NMDS", distance="bray")
MLE_nmds_16s_soil
stressplot(MLE_nmds_16s_soil)


# Extract vectors and incorporate sample data -------------------------------------------------
MLE_df_its_root <- data.frame(MLE_nmds_its_root$points)
MLE_df_16s_root <- data.frame(MLE_nmds_16s_root$points)

MLE_df_its_root <-
  left_join(
    tibble::rownames_to_column(
      as(MLE_physeq_fungi_root@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(MLE_df_its_root),
    by = "rowname"
  )

MLE_df_16s_root <-
  left_join(
    tibble::rownames_to_column(
      as(MLE_physeq_bact_root@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(MLE_df_16s_root),
    by = "rowname"
  )

rownames(MLE_df_its_root) <- MLE_df_its_root$rowname
head(MLE_df_its_root)

rownames(MLE_df_16s_root) <- MLE_df_16s_root$rowname
head(MLE_df_16s_root)


MLE_df_its_soil <- data.frame(MLE_nmds_its_soil$points)
MLE_df_16s_soil <- data.frame(MLE_nmds_16s_soil$points)

MLE_df_its_soil <-
  left_join(
    tibble::rownames_to_column(
      as(MLE_physeq_fungi_soil@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(MLE_df_its_soil),
    by = "rowname"
  )

MLE_df_16s_soil <-
  left_join(
    tibble::rownames_to_column(
      as(MLE_physeq_bact_soil@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(MLE_df_16s_soil),
    by = "rowname"
  )

rownames(MLE_df_its_soil) <- MLE_df_its_soil$rowname
head(MLE_df_its_soil)

rownames(MLE_df_16s_soil) <- MLE_df_16s_soil$rowname
head(MLE_df_16s_soil)

# rearranging the sample order ----------------------------------------------------------------
identical(rownames(MLE_df_its_root), rownames(MLE_df_16s_root))
MLE_order_root <- match(rownames(MLE_df_16s_root), rownames(MLE_df_its_root))
MLE_df_its_root <- MLE_df_its_root[MLE_order_root, ]

identical(rownames(MLE_df_its_soil), rownames(MLE_df_16s_soil))
MLE_order_soil <- match(rownames(MLE_df_16s_soil), rownames(MLE_df_its_soil))
MLE_df_its_soil <- MLE_df_its_soil[MLE_order_soil, ]



# Run procrustes analysis ---------------------------------------------------------------------
MLE_pro_root <- protest(
  MLE_df_its_root[, c("MDS1", "MDS2")],
  MLE_df_16s_root[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
MLE_pro_root
MLE_pro_root$ss

MLE_pro_soil <- protest(
  MLE_df_its_soil[, c("MDS1", "MDS2")],
  MLE_df_16s_soil[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
MLE_pro_soil
MLE_pro_soil$ss



# Plotting ------------------------------------------------------------------------------------
# m2 statistic provides an overall measure of the concordance between two data sets 
# (i.e. how close the two data configuration match),
# Procrustes analysis allows one to determine how much variance in one matrix is 
# attributable to the variance in the other. Further, visual inspection of a 
# Procrustes plot, in which the residuals between points from each matrix are
# mapped, can allow the identification of individual objects that have (relatively)
# unusual concordance.

MLE_pro_root_vectors=cbind(MLE_pro_root$Yrot,MLE_pro_root$X)
colnames(MLE_pro_root_vectors)=c("yNMDS1","yNMDS2","xNMDS1","xNMDS2")
MLE_df_pro_its_root=merge(MLE_df_its_root,MLE_pro_root_vectors,by="row.names")

MLE_pro_soil_vectors=cbind(MLE_pro_soil$Yrot,MLE_pro_soil$X)
colnames(MLE_pro_soil_vectors)=c("yNMDS1","yNMDS2","xNMDS1","xNMDS2")
MLE_df_pro_its_soil=merge(MLE_df_its_soil,MLE_pro_soil_vectors,by="row.names")


#####MLE Procrutes Errors####

residuals(MLE_pro_root)
length(residuals(MLE_pro_root))
MLE_res_root_plot=merge(data.frame(residuals(MLE_pro_root)),
                              data.frame(sample_data(MLE_physeq_bact_root)),
                              by="row.names")
head(MLE_res_root_plot)
colnames(MLE_res_root_plot)[2]="Res_Err"
dim(MLE_res_root_plot)
#110   38

#Analyses of the residuals

#MLE

#Roots



hist(MLE_res_root_plot$Res_Err)
qqPlot(MLE_res_root_plot$Res_Err)


MLE_res_ERR_root_plot_mod=lmer(sqrt(Res_Err)~siteID*FertStatus+(1|plotRep),data = MLE_res_root_plot)
plot(MLE_res_ERR_root_plot_mod)
hist(resid(MLE_res_ERR_root_plot_mod))
qqPlot(resid(MLE_res_ERR_root_plot_mod))
shapiro.test(resid(MLE_res_ERR_root_plot_mod))
#W = 0.98249, p-value = 0.1587
simulateResiduals(fittedModel = MLE_res_ERR_root_plot_mod, plot = T)

anova(MLE_res_ERR_root_plot_mod)
#siteID            0.229877 0.057469     4 97.978 27.0276 4.043e-15 ***
#FertStatus        0.003595 0.003595     1 97.417  1.6908   0.19657    
#siteID:FertStatus 0.027958 0.006990     4 97.355  3.2872   0.01425 * 




emmeans(MLE_res_ERR_root_plot_mod,pairwise~siteID)

#contrast  estimate     SE   df t.ratio p.value
#ESC - HAN  0.00536 0.0147 99.1   0.364  0.9962
#ESC - LC  -0.02115 0.0138 97.3  -1.535  0.5424
#ESC - LUX  0.02445 0.0136 97.0   1.796  0.3816
#ESC - RHN  0.10715 0.0135 97.1   7.956  <.0001
#HAN - LC  -0.02651 0.0150 99.7  -1.767  0.3985
#HAN - LUX  0.01909 0.0147 99.1   1.296  0.6945
#HAN - RHN  0.10179 0.0147 99.5   6.937  <.0001
#LC - LUX   0.04560 0.0138 97.3   3.310  0.0112
#LC - RHN   0.12830 0.0136 97.1   9.422  <.0001
#LUX - RHN  0.08270 0.0135 97.1   6.140  <.0001


emmeans(MLE_res_ERR_root_plot_mod,pairwise~FertStatus|siteID)
#siteID = RHN:
#  contrast      estimate     SE   df t.ratio p.value
#Fert - Unfert  0.06485 0.0188 97.0   3.445  0.0008





#Soil

residuals(MLE_pro_soil)
length(residuals(MLE_pro_soil))
MLE_res_soil_plot=merge(data.frame(residuals(MLE_pro_soil)),
                              data.frame(sample_data(MLE_physeq_bact_soil)),
                              by="row.names")
head(MLE_res_soil_plot)
colnames(MLE_res_soil_plot)[2]="Res_Err"
dim(MLE_res_soil_plot)
#112   38

hist(MLE_res_soil_plot$Res_Err)
qqPlot(MLE_res_soil_plot$Res_Err)


MLE_res_ERR_soil_plot_mod=lmer((Res_Err)~siteID*FertStatus+(1|plotRep),data = MLE_res_soil_plot)
plot(MLE_res_ERR_soil_plot_mod)
hist(resid(MLE_res_ERR_soil_plot_mod))
qqPlot(resid(MLE_res_ERR_soil_plot_mod))
shapiro.test(resid(MLE_res_ERR_soil_plot_mod))
#W = 0.99035, p-value = 0.6163

simulateResiduals(fittedModel = MLE_res_ERR_soil_plot_mod, plot = T)

anova(MLE_res_ERR_soil_plot_mod)
#siteID            0.0274431 0.0068608     4 99.969 23.8419 7.303e-14 ***
#FertStatus        0.0003780 0.0003780     1 99.224  1.3136   0.25450    
#siteID:FertStatus 0.0023218 0.0005804     4 99.203  2.0171   0.09784 .  




emmeans(MLE_res_ERR_soil_plot_mod,pairwise~siteID)

#contrast              estimate      SE    df t.ratio p.value
# ESC - HAN -0.04632 0.00539 101.5  -8.587  <.0001
#ESC - LC  -0.00600 0.00490  99.0  -1.225  0.7369
#ESC - LUX -0.00870 0.00495  99.1  -1.757  0.4045
#ESC - RHN -0.00165 0.00495  99.1  -0.333  0.9973
#HAN - LC   0.04033 0.00539 101.5   7.475  <.0001
#HAN - LUX  0.03762 0.00546 101.6   6.895  <.0001
#HAN - RHN  0.04467 0.00542 101.0   8.244  <.0001
#LC - LUX  -0.00271 0.00495  99.1  -0.546  0.9822
#LC - RHN   0.00435 0.00495  99.1   0.877  0.9046
#LUX - RHN  0.00705 0.00501  99.2   1.407  0.6247


summary(MLE_res_root_plot)
summary(MLE_res_soil_plot)

#


#####MLE Procrutes Figure####
# Plotting ------------------------------------------------------------------------------------
# m2 statistic provides an overall measure of the concordance between two data sets 
# (i.e. how close the two data configuration match),
# Procrustes analysis allows one to determine how much variance in one matrix is 
# attributable to the variance in the other. Further, visual inspection of a 
# Procrustes plot, in which the residuals between points from each matrix are
# mapped, can allow the identification of individual objects that have (relatively)
# unusual concordance.

brewer.pal(5,"Set1")
site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")

site_colors=c("ESC"="#E41A1C", "HAN"="#377EB8", "LC"="#4DAF4A", 
              "LUX"="#984EA3", "RHN"="#FF7F00")

GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_MLE_G5,Root_soil=="Soil")
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil=
  prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)>0,GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)
nsamples(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil)
#112
set.seed(2021)
GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil,method = "NMDS")
#*** Solution reached
#0.1712595   

(MLE_compart_fung_soil_legend=plot_ordination(GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil,GLBRC018_OTU_fung_MMPRNT_MLE_G5_soil_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

pairwise_MLE_fert_root_procrustes=data.frame(y_bot=c(0.054),
                                                   x_min=c(4.8),
                                                   x_max=c(5.2),
                                                   annot_text=c("***"))

MLE_root_procrustes_letters=data.frame(siteID=c("LUX","LC",
                                                      "ESC","HAN","RHN"),
                                             sig_let=c("a","b","ab","ab","c"))

MLE_root_procrustes_max=MLE_res_root_plot%>%group_by(siteID)%>%summarise(max_value=max(Res_Err))

MLE_root_procrustes_letters_max=merge(MLE_root_procrustes_letters,MLE_root_procrustes_max,by="siteID")



(MLE_Root_PRO_nmds=ggplot(MLE_df_pro_its_root) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS1") + # to adjust decimals
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS2",
                       limits = c(-0.15,0.13)) +
    scale_shape_manual("", 
                       values = c(19, 2), 
                       labels=c("Fert","Unfert")) +
    theme(panel.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line = element_line()) +
    geom_point(aes(x=yNMDS1, y=yNMDS2, 
                   colour=siteID, shape=FertStatus), size=3, alpha=0.9) +
    geom_point(aes(x=xNMDS1, y=xNMDS2, 
                   colour=siteID, shape=FertStatus), size=3, alpha=0.9) +
    geom_segment(aes(x=xNMDS1,y=xNMDS2,xend=yNMDS1,yend=yNMDS2,
                     colour=siteID), size=0.3,
                 arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
    geom_hline(yintercept = 0,linetype="dashed") +
    geom_vline(xintercept = 0,linetype="dashed") +
    geom_abline(slope = -1/MLE_pro_root$rotation[1,2]) +
    geom_abline(slope = MLE_pro_root$rotation[1,2])+theme_bw(base_size = 20)+
    scale_color_brewer(palette = "Set1",name=NULL,
                       labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    annotate("text", -Inf, Inf, label = paste("italic(m) ^ 2 == ",signif(MLE_pro_root$ss,3),sep = ""), 
             parse = TRUE, size = 5, hjust = -1.37, vjust = 2) +
    annotate("text", -Inf, Inf, label = paste("Stress bacteria = ",signif(MLE_nmds_16s_root$stress,3),sep = ""), 
             parse = F, size = 5, hjust = -0.1, vjust = 5) +
    annotate("text", -Inf, Inf, label = paste("Stress fungi = ",signif(MLE_nmds_its_root$stress,3),sep = ""), 
             parse = F, size = 5, hjust = -0.26, vjust = 7) +ggtitle("Root")+
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 32)))




(MLE_Soil_PRO_nmds=ggplot(MLE_df_pro_its_soil) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS1") + # to adjust decimals
    scale_y_continuous(labels = NULL,name=NULL,
                       limits = c(-0.15,0.13)) +
    scale_shape_manual("", 
                       values = c(19, 2), 
                       labels=c("Fert","Unfert")) +
    theme(panel.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line = element_line()) +
    geom_point(aes(x=yNMDS1, y=yNMDS2, 
                   colour=siteID, shape=FertStatus), size=3, alpha=0.9) +
    geom_point(aes(x=xNMDS1, y=xNMDS2, 
                   colour=siteID, shape=FertStatus), size=3, alpha=0.9) +
    geom_segment(aes(x=xNMDS1,y=xNMDS2,xend=yNMDS1,yend=yNMDS2,
                     colour=siteID), size=0.3,
                 arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
    geom_hline(yintercept = 0,linetype="dashed") +
    geom_vline(xintercept = 0,linetype="dashed") +
    geom_abline(slope = -1/MLE_pro_soil$rotation[1,2]) +
    geom_abline(slope = MLE_pro_soil$rotation[1,2])+theme_bw(base_size = 20)+
    scale_color_brewer(palette = "Set1",name=NULL,
                       labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    annotate("text", -Inf, Inf, label = paste("italic(m) ^ 2 == ",signif(MLE_pro_soil$ss,3),sep = ""), 
             parse = TRUE, size = 5, hjust = -1.37, vjust = 2) +
    annotate("text", -Inf, Inf, label = paste("Stress bacteria = ",signif(MLE_nmds_16s_soil$stress,3),sep = ""), 
             parse = F, size = 5, hjust = -0.1, vjust = 5) +
    annotate("text", -Inf, Inf, label = paste("Stress fungi = ",signif(MLE_nmds_its_soil$stress,3),sep = ""), 
             parse = F, size = 5, hjust = -0.26, vjust = 7) +ggtitle("Soil")+
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 32)))


(procrustes_NMDS_MLE_2panel=plot_grid(MLE_Root_PRO_nmds,MLE_Soil_PRO_nmds,ncol = 2,rel_widths = c(1,0.8),
                                      labels = c("a)","b)"),label_size = 32,label_x = c(0,-0.05)))

plot_grid(procrustes_NMDS_MLE_2panel,get_legend(MLE_compart_fung_soil_legend),ncol = 2,rel_widths = c(4,0.5))


(MLE_procrustes_p3=ggplot(MLE_res_root_plot)+
    geom_boxplot(aes(y=Res_Err, x=factor(siteID,levels = site_order),
                     fill=FertStatus,color=FertStatus))+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_fill_manual(values = c("black","grey"), name=NULL)+
    scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    scale_y_continuous(name = "Residual error",limits = c(0,0.16))+theme(axis.text.x = element_blank()))


(MLE_procrustes_root_p3=ggplot(MLE_res_root_plot)+
    geom_boxplot(aes(y=Res_Err, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+theme_cowplot(font_size = 24)+
    geom_text(data = MLE_root_procrustes_letters_max, aes(x=siteID, y =max_value, label = sig_let), 
              vjust=0, size=10,color="black",alpha=1)+
    geom_signif(data=pairwise_MLE_fert_root_procrustes,aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 10,manual = T)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_x_discrete(labels=site_labels,name=NULL)+
    scale_alpha_manual(values = c(1,0.1))+
    scale_y_continuous(name = "Residual error",limits = c(0,0.16))+
    theme(legend.position = "none"))


MLE_soil_procrustes_letters=data.frame(siteID=c("LUX","LC",
                                                      "ESC","HAN","RHN"),
                                             sig_let=c("A","A","A","B","A"))
MLE_soil_procrustes_max=MLE_res_soil_plot%>%group_by(siteID)%>%summarise(max_value=max(Res_Err))

MLE_soil_procrustes_letters_max=merge(MLE_soil_procrustes_letters,MLE_soil_procrustes_max,by="siteID")


(MLE_procrustes_soil_p3=ggplot(MLE_res_soil_plot)+
    geom_boxplot(aes(y=Res_Err, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+theme_cowplot(font_size = 24)+
    geom_text(data = MLE_soil_procrustes_letters_max, aes(x=siteID, y =max_value, label = sig_let), 
              vjust=0, size=10,color="black",alpha=1)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_x_discrete(labels=site_labels,name=NULL)+
    scale_alpha_manual(values = c(1,0.1))+
    scale_y_continuous(name = "Residual error",limits = c(0,0.16))+
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "none"))


(procrustes_MLE_2panel2=plot_grid(MLE_procrustes_root_p3,MLE_procrustes_soil_p3,ncol = 2,rel_widths = c(1,0.8),
                                  labels = c("c)","d)"),label_size = 32,label_x = c(0,-0.2)))

plot_grid(procrustes_MLE_2panel2,get_legend(MLE_procrustes_p3),ncol = 2,rel_widths = c(4,0.5))

(procrustes_NMDS_MLE_2panel=plot_grid(MLE_Root_PRO_nmds,MLE_Soil_PRO_nmds,ncol = 2,rel_widths = c(1,0.8),
                                      labels = c("a)","b)"),label_size = 32,label_x = c(0,-0.05)))

plot_grid(procrustes_NMDS_MLE_2panel,get_legend(MLE_compart_fung_soil_legend),ncol = 2,rel_widths = c(4,0.5))


plot_grid(MLE_Root_PRO_nmds,MLE_Soil_PRO_nmds,get_legend(MLE_compart_fung_soil_legend),
          MLE_procrustes_root_p3,MLE_procrustes_soil_p3,get_legend(MLE_procrustes_p3),
          labels = c("a)","b)"," ","c)","d)"," "),label_size = 32,label_x = c(0,0.1,0,0,0.1,0),
          rel_widths = c(2,2,0.5),align = "v",axis = "l")


#NOT INCLUDED IN REPOSITORY
ggsave(plot_grid(MLE_Root_PRO_nmds,MLE_Soil_PRO_nmds,get_legend(MLE_compart_fung_soil_legend),
                 MLE_procrustes_root_p3,MLE_procrustes_soil_p3,get_legend(MLE_procrustes_p3),
                 labels = c("a)","b)"," ","c)","d)"," "),label_size = 32,label_x = c(0,0.05,0,0,0.05,0),
                 rel_widths = c(2,2,0.5),align = "v",axis = "l",rel_heights = c(1.2,1)),
       filename = "procrustes_NMDS_boxplot_All_Sites_date_color_p.png",path = here::here("Manuscript","MLE_comm_figs"),width = 22,height =15)

ggsave(plot_grid(MLE_Root_PRO_nmds,MLE_Soil_PRO_nmds,get_legend(MLE_compart_fung_soil_legend),
                 MLE_procrustes_root_p3,MLE_procrustes_soil_p3,get_legend(MLE_procrustes_p3),
                 labels = c("a)","b)"," ","c)","d)"," "),label_size = 32,label_x = c(0,0.05,0,0,0.05,0),
                 rel_widths = c(2,2,0.5),align = "v",axis = "l",rel_heights = c(1.2,1)),
       filename = "procrustes_NMDS_boxplot_All_Sites_date_color_p.svg",path = here::here("Manuscript","MLE_comm_figs"),width = 22,height =15)

#NOT INCLUDED IN REPOSITORY

