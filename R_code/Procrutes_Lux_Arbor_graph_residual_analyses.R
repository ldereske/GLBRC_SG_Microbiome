
## ---------------------------
##
## Script name: Community correlation between bacteria and fungi in the switchgrass microbiomes across 
## the Lux Arbor Growing Season
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
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(ggsignif)
library(vegan)
library(car)

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






#Procrustes Lux Arbor------------
# Lux Arbor MATCH SAMPLES -------------------------------------------------------------------------------
# Match samples to the same ones for fungi and bacteria ---------------------------------------

#Code Selection Written by Gian Maria Niccolò Benucci
#Modified by Lukas Bell-Dereske

# Lux Arbor time series
set.seed(2021)


physeq_bact_LUX_root <-
  GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil %>%
  subset_samples(sampType %in% "Root") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)

physeq_fungi_LUX_root <-
  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5 %>%
  subset_samples(Root_soil %in% "Root") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)


physeq_bact_LUX_soil <-
  GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil %>%
  subset_samples(sampType %in% "Soil") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)

physeq_fungi_LUX_soil <-
  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5 %>%
  subset_samples(Root_soil %in% "Soil") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)


intersect(
  physeq_bact_LUX_root@sam_data$sampleID_long,
  physeq_fungi_LUX_root@sam_data$sampleID_long
) -> good_time_root
length(good_time_root)
#139

intersect(
  physeq_bact_LUX_soil@sam_data$sampleID_long,
  physeq_fungi_LUX_soil@sam_data$sampleID_long
) -> good_time_soil
length(good_time_soil)
#342

physeq_bact_LUX_root_filt <-
  subset_samples(physeq_bact_LUX_root, sampleID_long %in% good_time_root) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_bact_LUX_root_filt

physeq_fungi_LUX_root_filt <-
  subset_samples(physeq_fungi_LUX_root, sampleID_long %in% good_time_root) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_fungi_LUX_root_filt

physeq_bact_LUX_soil_filt <-
  subset_samples(physeq_bact_LUX_soil, sampleID_long %in% good_time_soil) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_bact_LUX_soil_filt

physeq_fungi_LUX_soil_filt <-
  subset_samples(physeq_fungi_LUX_soil, sampleID_long %in% good_time_soil) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_fungi_LUX_soil_filt

sample_names(physeq_bact_LUX_root_filt) <-
  physeq_bact_LUX_root_filt@sam_data$sampleID_long
physeq_bact_LUX_root_filt@sam_data

sample_names(physeq_fungi_LUX_root_filt) <-
  physeq_fungi_LUX_root_filt@sam_data$sampleID_long
physeq_fungi_LUX_root_filt@sam_data

sample_names(physeq_bact_LUX_soil_filt) <-
  physeq_bact_LUX_soil_filt@sam_data$sampleID_long
physeq_bact_LUX_soil_filt@sam_data

sample_names(physeq_fungi_LUX_soil_filt) <-
  physeq_fungi_LUX_soil_filt@sam_data$sampleID_long
physeq_fungi_LUX_soil_filt@sam_data

#NMDS
# LUX arbor Bacteria and Fungi
set.seed(2021)
nmds_16s_soil_time <- ordinate(
  physeq_bact_LUX_soil_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_16s_soil_time
stressplot(nmds_16s_soil_time)
plot(nmds_16s_soil_time, "sites")

nmds_16s_root_time <- ordinate(
  physeq_bact_LUX_root_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_16s_root_time
stressplot(nmds_16s_root_time)

nmds_its_soil_time <- ordinate(
  physeq_fungi_LUX_soil_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_its_soil_time
plot(nmds_its_soil_time, "sites")
stressplot(nmds_its_soil_time)

nmds_its_root_time <- ordinate(
  physeq_fungi_LUX_root_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_its_root_time
stressplot(nmds_its_root_time)
plot(nmds_its_root_time, "sites")




# Lux Arbor 
df_its_root_time <- data.frame(nmds_its_root_time$points)
df_16s_root_time <- data.frame(nmds_16s_root_time$points)

df_its_root_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_fungi_LUX_root_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus,plotRep)
    ),
    tibble::rownames_to_column(df_its_root_time),
    by = "rowname"
  )

df_16s_root_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_bact_LUX_root_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus,plotRep)
    ),
    tibble::rownames_to_column(df_16s_root_time),
    by = "rowname"
  )

rownames(df_its_root_time) <- df_its_root_time$rowname
head(df_its_root_time)

rownames(df_16s_root_time) <- df_16s_root_time$rowname
head(df_16s_root_time)

df_its_soil_time <- data.frame(nmds_its_soil_time$points)
df_16s_soil_time <- data.frame(nmds_16s_soil_time$points)

df_its_soil_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_fungi_LUX_soil_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus,plotRep)
    ),
    tibble::rownames_to_column(df_its_soil_time),
    by = "rowname"
  )

df_16s_soil_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_bact_LUX_soil_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus,plotRep)
    ),
    tibble::rownames_to_column(df_16s_soil_time),
    by = "rowname"
  )

rownames(df_its_soil_time) <- df_its_soil_time$rowname
head(df_its_soil_time)

rownames(df_16s_soil_time) <- df_16s_soil_time$rowname
head(df_16s_soil_time)





identical(rownames(df_its_root_time), rownames(df_16s_root_time))
order_root_time <-
  match(rownames(df_16s_root_time), rownames(df_its_root_time))
df_its_root_time <- df_its_root_time[order_root_time, ]

identical(rownames(df_its_soil_time), rownames(df_16s_soil_time))
order_soil_time <-
  match(rownames(df_16s_soil_time), rownames(df_its_soil_time))
df_its_soil_time <- df_its_soil_time[order_soil_time, ]



# Lux Arbor 
pro_root_time <- protest(
  df_its_root_time[, c("MDS1", "MDS2")],
  df_16s_root_time[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
pro_root_time
pro_root_time$ss

pro_soil_time <- protest(
  df_its_soil_time[, c("MDS1", "MDS2")],
  df_16s_soil_time[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
pro_soil_time
pro_soil_time$ss

pro_root_time_vectors=cbind(pro_root_time$Yrot,pro_root_time$X)
colnames(pro_root_time_vectors)=c("yNMDS1","yNMDS2","xNMDS1","xNMDS2")
df_pro_its_root_time=merge(df_its_root_time,pro_root_time_vectors,by="row.names")


Date_F <- function(x){
  format(as.Date(x, origin = '1970-01-01',"%Y-%m-%d"),"%b")
}

head(df_pro_its_root_time)

ggplot(df_pro_its_root_time) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS1") + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS2") +
  scale_shape_manual("", 
                     values = c(19, 2), 
                     labels=c("Fert","Unfert")) +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line()) +
  geom_point(aes(x=yNMDS1, y=yNMDS2, 
                 colour=as.Date(collectionDate, format="%m/%d/%Y"), 
                 shape=FertStatus), size=1.5, alpha=0.9) +
  geom_point(aes(x=xNMDS1, y=xNMDS2, 
                 colour=as.Date(collectionDate, format="%m/%d/%Y"),
                 shape=FertStatus), size=1.5, alpha=0.9) +
  geom_segment(aes(x=xNMDS1,y=xNMDS2,xend=yNMDS1,yend=yNMDS2,
                   colour=as.Date(collectionDate, format="%m/%d/%Y")), size=0.3,
               arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_abline(slope = -1/pro_root_time$rotation[1,2]) +
  geom_abline(slope = pro_root_time$rotation[1,2])+theme_bw(base_size = 20)+
  scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name=NULL,label=Date_F,n.breaks=3)+
  annotate("text", -Inf, Inf, label = paste("italic(m) ^ 2 == ",signif(pro_root_time$ss,3),sep = ""), 
           parse = TRUE, size = 3, hjust = -2.37, vjust = 2) +
  annotate("text", -Inf, Inf, label = paste("Stress bacteria = ",signif(nmds_16s_root_time$stress,3),sep = ""), 
           parse = F, size = 3, hjust = -1.1, vjust = 5) +
  annotate("text", -Inf, Inf, label = paste("Stress fungi = ",signif(nmds_its_root_time$stress,3),sep = ""), 
           parse = F, size = 3, hjust = -1.26, vjust = 7) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())


pro_soil_time_vectors=cbind(pro_soil_time$Yrot,pro_soil_time$X)
colnames(pro_soil_time_vectors)=c("yNMDS1","yNMDS2","xNMDS1","xNMDS2")
df_pro_its_soil_time=merge(df_its_soil_time,pro_soil_time_vectors,by="row.names")

head(df_pro_its_soil_time)

ggplot(df_pro_its_soil_time) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS1") + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS2") +
  scale_shape_manual("", 
                     values = c(19, 2), 
                     labels=c("Fert","Unfert")) +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line()) +
  geom_point(aes(x=yNMDS1, y=yNMDS2, 
                 colour=as.Date(collectionDate, format="%m/%d/%Y"), 
                 shape=FertStatus), size=1.5, alpha=0.9) +
  geom_point(aes(x=xNMDS1, y=xNMDS2, 
                 colour=as.Date(collectionDate, format="%m/%d/%Y"),
                 shape=FertStatus), size=1.5, alpha=0.9) +
  geom_segment(aes(x=xNMDS1,y=xNMDS2,xend=yNMDS1,yend=yNMDS2,
                   colour=as.Date(collectionDate, format="%m/%d/%Y")), size=0.3,
               arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_abline(slope = -1/pro_soil_time$rotation[1,2]) +
  geom_abline(slope = pro_soil_time$rotation[1,2])+theme_bw(base_size = 20)+
  scale_color_gradient(low = "#FFFF99",high = "#B15928",name=NULL,label=Date_F,n.breaks=4)+ 
  annotate("text", -Inf, Inf, label = paste("italic(m) ^ 2 == ",signif(pro_soil_time$ss,3),sep = ""), 
           parse = TRUE, size = 3, hjust = -2.37, vjust = 2) +
  annotate("text", -Inf, Inf, label = paste("Stress bacteria = ",signif(nmds_16s_soil_time$stress,3),sep = ""), 
           parse = F, size = 3, hjust = -1.1, vjust = 5) +
  annotate("text", -Inf, Inf, label = paste("Stress fungi = ",signif(nmds_its_soil_time$stress,3),sep = ""), 
           parse = F, size = 3, hjust = -1.26, vjust = 7) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#Procrustes errors Lux Arbor-----

#Roots

length(residuals(pro_root_time))
#139
res_root_time_plot=merge(data.frame(residuals(pro_root_time)),
                         data.frame(sample_data(physeq_bact_LUX_root_filt)),
                         by="row.names")


head(res_root_time_plot)
colnames(res_root_time_plot)[colnames(res_root_time_plot)=="residuals.pro_root_time."]="Res_Err"
dim(res_root_time_plot)
#139   38

hist(res_root_time_plot$Res_Err)
qqPlot(res_root_time_plot$Res_Err)


res_ERR_root_time_plot_mod=lmer(sqrt(Res_Err)~collectionDate*FertStatus+(1|plotRep),data = res_root_time_plot)
plot(res_ERR_root_time_plot_mod)
hist(resid(res_ERR_root_time_plot_mod))
qqPlot(resid(res_ERR_root_time_plot_mod))
shapiro.test(resid(res_ERR_root_time_plot_mod))
#W = 0.98877, p-value = 0.324
simulateResiduals(fittedModel = res_ERR_root_time_plot_mod, plot = T)

anova(res_ERR_root_time_plot_mod)
#collectionDate            0.041138 0.008228     5   127  2.4315 0.0384418 *  
#FertStatus                0.040384 0.040384     1   127 11.9349 0.0007494 ***
#collectionDate:FertStatus 0.051017 0.010203     5   127  3.0154 0.0131870 * 




emmeans(res_ERR_root_time_plot_mod,pairwise~collectionDate)
#$contrasts
#contrast                  estimate     SE  df t.ratio p.value
#(10/3/2018) - (5/29/2018)  0.03653 0.0174 124   2.104  0.2920
#(10/3/2018) - (6/25/2018)  0.01189 0.0170 124   0.700  0.9817
#(10/3/2018) - (7/30/2018)  0.03290 0.0172 124   1.915  0.3983
#(10/3/2018) - (8/20/2018) -0.01261 0.0170 124  -0.742  0.9762
#(10/3/2018) - (9/17/2018)  0.01674 0.0172 124   0.974  0.9253
#(5/29/2018) - (6/25/2018) -0.02464 0.0172 124  -1.435  0.7060
#(5/29/2018) - (7/30/2018) -0.00363 0.0174 124  -0.209  0.9999
#(5/29/2018) - (8/20/2018) -0.04914 0.0172 124  -2.861  0.0547
#(5/29/2018) - (9/17/2018) -0.01979 0.0174 124  -1.139  0.8642
#(6/25/2018) - (7/30/2018)  0.02102 0.0170 124   1.237  0.8178
#(6/25/2018) - (8/20/2018) -0.02450 0.0168 124  -1.459  0.6909
#(6/25/2018) - (9/17/2018)  0.00485 0.0170 124   0.286  0.9997
#(7/30/2018) - (8/20/2018) -0.04551 0.0170 124  -2.679  0.0868
#(7/30/2018) - (9/17/2018) -0.01616 0.0172 124  -0.941  0.9349
#(8/20/2018) - (9/17/2018)  0.02935 0.0170 124   1.728  0.5163

emmeans(res_ERR_root_time_plot_mod,pairwise~FertStatus|collectionDate)
#$contrasts
#collectionDate = 10/3/2018:
#  contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert   0.0683 0.0243 124   2.811  0.0057

#collectionDate = 5/29/2018:
#  contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert  -0.0460 0.0248 125  -1.851  0.0666

#collectionDate = 7/30/2018:
#contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert   0.0515 0.0243 124   2.120  0.0360

#collectionDate = 8/20/2018:
#  contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert   0.0617 0.0237 124   2.600  0.0105

#collectionDate = 9/17/2018:
#  contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert   0.0498 0.0243 124   2.050  0.0425



#Soil
pro_soil_time

length(residuals(pro_soil_time))
#342
res_soil_time_plot=merge(data.frame(residuals(pro_soil_time)),
                         data.frame(sample_data(physeq_bact_LUX_soil_filt)),
                         by="row.names")


head(res_soil_time_plot)
colnames(res_soil_time_plot)[colnames(res_soil_time_plot)=="residuals.pro_soil_time."]="Res_Err"
dim(res_soil_time_plot)
#342   38

hist(res_soil_time_plot$Res_Err)
qqPlot(res_soil_time_plot$Res_Err)


res_ERR_soil_time_plot_mod=lmer(log(Res_Err)~collectionDate*FertStatus+(1|plotRep),data = res_soil_time_plot)
plot(res_ERR_soil_time_plot_mod)
hist(resid(res_ERR_soil_time_plot_mod))
qqPlot(resid(res_ERR_soil_time_plot_mod))
shapiro.test(resid(res_ERR_soil_time_plot_mod))
#W = 0.98909, p-value = 0.01171
simulateResiduals(fittedModel = res_ERR_soil_time_plot_mod, plot = T)

anova(res_ERR_soil_time_plot_mod)
#collectionDate            6.0054 0.42895    14 309.08  1.0719 0.38260  
#FertStatus                1.8043 1.80433     1 309.06  4.5086 0.03452 *
#collectionDate:FertStatus 4.1438 0.29599    14 309.06  0.7396 0.73383  



emmeans(res_ERR_soil_time_plot_mod,pairwise~FertStatus)
emmeans(res_ERR_soil_time_plot_mod,pairwise~FertStatus|collectionDate)


#Procrustes combine figures Lux Arbor----

Date_F <- function(x){
  format(as.Date(x, origin = '1970-01-01',"%Y-%m-%d"),"%b")
}



(lux_procruste_p3=ggplot(res_root_time_plot, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=Res_Err,
                                                 fill=FertStatus,color=FertStatus))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
    scale_fill_manual(values = c("black","lightgrey"), name=NULL)+scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    ylab("procrusteersion")+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(axis.text.x = element_blank()))




pairwise_LUX_fert_root_procrustes=data.frame(y_bot=c(0.105,0.08,0.11,0.091,0.095),
                                             x_min=as.Date(c("2018-05-26","2018-07-27","2018-08-17","2018-09-14","2018-09-30")),
                                             x_max=as.Date(c("2018-06-01","2018-08-02","2018-08-23","2018-09-20","2018-10-07")),
                                             annot_text=c("#","*"," * ","  *  ","**"))


(LUX_root_PRO_nmds=ggplot(df_pro_its_root_time) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS1") + # to adjust decimals
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS2") +
    scale_shape_manual(values = c(19, 2), 
                       labels=c("Fert","Unfert"),
                       name=NULL) +
    theme(panel.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line = element_line()) +
    geom_point(aes(x=yNMDS1, y=yNMDS2, 
                   colour=as.Date(collectionDate, format="%m/%d/%Y"), 
                   shape=FertStatus), size=3, alpha=0.9) +
    geom_point(aes(x=xNMDS1, y=xNMDS2, 
                   colour=as.Date(collectionDate, format="%m/%d/%Y"),
                   shape=FertStatus), size=3, alpha=0.9) +
    geom_segment(aes(x=xNMDS1,y=xNMDS2,xend=yNMDS1,yend=yNMDS2,
                     colour=as.Date(collectionDate, format="%m/%d/%Y")), size=0.3,
                 arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
    geom_hline(yintercept = 0,linetype="dashed") +
    geom_vline(xintercept = 0,linetype="dashed") +
    geom_abline(slope = -1/pro_root_time$rotation[1,2]) +
    geom_abline(slope = pro_root_time$rotation[1,2])+theme_bw(base_size = 20)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name=NULL,label=Date_F,n.breaks=3)+
    annotate("text", -Inf, Inf, label = paste("italic(m) ^ 2 == ",signif(pro_root_time$ss,3),sep = ""), 
             parse = TRUE, size = 3, hjust = -1.37, vjust = 2) +
    annotate("text", -Inf, Inf, label = paste("Stress bacteria = ",signif(nmds_16s_root_time$stress,3),sep = ""), 
             parse = F, size = 3, hjust = -0.1, vjust = 5) +
    annotate("text", -Inf, Inf, label = paste("Stress fungi = ",signif(nmds_its_root_time$stress,3),sep = ""), 
             parse = F, size = 3, hjust = -0.26, vjust = 7) +ggtitle("Root")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 32)))



(LUX_soil_PRO_nmds=ggplot(df_pro_its_soil_time) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS1") + # to adjust decimals
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05),name="NMDS2") +
    scale_shape_manual("", 
                       values = c(19, 2), 
                       labels=c("Fert","Unfert")) +
    theme(panel.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line = element_line()) +
    geom_point(aes(x=yNMDS1, y=yNMDS2, 
                   colour=as.Date(collectionDate, format="%m/%d/%Y"), 
                   shape=FertStatus), size=3, alpha=0.9) +
    geom_point(aes(x=xNMDS1, y=xNMDS2, 
                   colour=as.Date(collectionDate, format="%m/%d/%Y"),
                   shape=FertStatus), size=3, alpha=0.9) +
    geom_segment(aes(x=xNMDS1,y=xNMDS2,xend=yNMDS1,yend=yNMDS2,
                     colour=as.Date(collectionDate, format="%m/%d/%Y")), size=0.3,
                 arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
    geom_hline(yintercept = 0,linetype="dashed") +
    geom_vline(xintercept = 0,linetype="dashed") +
    geom_abline(slope = -1/pro_soil_time$rotation[1,2]) +
    geom_abline(slope = pro_soil_time$rotation[1,2])+theme_bw(base_size = 20)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name=NULL,label=Date_F,n.breaks=4)+ 
    annotate("text", -Inf, Inf, label = paste("italic(m) ^ 2 == ",signif(pro_soil_time$ss,3),sep = ""), 
             parse = TRUE, size = 3, hjust = -1.37, vjust = 2) +
    annotate("text", -Inf, Inf, label = paste("Stress bacteria = ",signif(nmds_16s_soil_time$stress,3),sep = ""), 
             parse = F, size = 3, hjust = -0.1, vjust = 5) +
    annotate("text", -Inf, Inf, label = paste("Stress fungi = ",signif(nmds_its_soil_time$stress,3),sep = ""), 
             parse = F, size = 3, hjust = -0.26, vjust = 7) +ggtitle("Soil")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 32),axis.title.y = element_blank(),legend.key.width = unit(2.5,"line")))

(procrustes_NMDS_LUX_2panel=plot_grid(LUX_root_PRO_nmds,LUX_soil_PRO_nmds,ncol = 2,rel_widths = c(1,0.8),
                                      labels = c("a)","b)"),label_size = 32,label_x = c(0,-0.05)))





(LUX_root_procruste_p3=ggplot(res_root_time_plot, 
                              aes(y=Res_Err))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(x=as.Date(collectionDate, format="%m/%d/%Y"),fill=as.Date(collectionDate, format="%m/%d/%Y"),
                     group=interaction(collectionDate,FertStatus),
                     alpha=FertStatus),position = position_dodge(12))+
    scale_fill_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_color_gradient(low = "#B2DF8A",high = "#33A02C",name="Date",label=Date_F,n.breaks=3)+
    scale_alpha_manual(values = c(1,0.4))+
    geom_signif(data=pairwise_LUX_fert_root_procrustes,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    scale_y_continuous(name="Residual error",limits = c(0,0.21))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="6 week",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32)))

(LUX_soil_procruste_p3=ggplot(res_soil_time_plot, 
                              aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=Res_Err,
                                  fill=as.Date(collectionDate, format="%m/%d/%Y"),
                                  alpha=FertStatus))+
    scale_shape_manual(values=c(19,24),name=NULL)+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(9))+
    scale_fill_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_color_gradient(low = "#FFFF99",high = "#B15928",name="Date",label=Date_F,n.breaks=4)+
    scale_y_continuous(name="Residual error",limits = c(0,0.21))+
    scale_alpha_manual(values = c(1,0.4))+
    theme_cowplot(font_size = 24)+scale_x_date(date_breaks ="2 month",date_labels = "%b",name = NULL)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 32),
          axis.title.y = element_blank(),axis.text.y = element_blank()))


(procrustes_LUX_2panel2=plot_grid(LUX_root_procruste_p3,LUX_soil_procruste_p3,ncol = 2,rel_widths = c(1,0.8),
                                  labels = c("c)","d)"),label_size = 32,label_x = c(0,-0.2)))

plot_grid(procrustes_LUX_2panel2,get_legend(lux_procruste_p3),ncol = 2,rel_widths = c(4,0.5))



plot_grid(LUX_root_PRO_nmds,LUX_soil_PRO_nmds,get_legend(LUX_soil_procruste_p3),
          LUX_root_procruste_p3,LUX_soil_procruste_p3,get_legend(lux_procruste_p3),
          labels = c("a)","b)"," ","c)","d)"," "),label_size = 32,label_x = c(0,0.1,0,0,0.1,0),
          rel_widths = c(2,2,0.5),align = "v",axis = "l")


#NOT INCLUDED IN REPOSITORY
ggsave(plot_grid(LUX_root_PRO_nmds,LUX_soil_PRO_nmds,get_legend(LUX_soil_procruste_p3),
                 LUX_root_procruste_p3,LUX_soil_procruste_p3,get_legend(lux_procruste_p3),
                 labels = c("a)","b)"," ","c)","d)"," "),label_size = 32,label_x = c(0,0.05,0,0,0.05,0),
                 rel_widths = c(2,2,0.5),align = "v",axis = "l",rel_heights = c(1.2,1)),
       filename = "procrustes_NMDS_boxplot_LUX_date_color_p.png",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 15,height =12)

ggsave(plot_grid(LUX_root_PRO_nmds,LUX_soil_PRO_nmds,get_legend(LUX_soil_procruste_p3),
                 LUX_root_procruste_p3,LUX_soil_procruste_p3,get_legend(lux_procruste_p3),
                 labels = c("a)","b)"," ","c)","d)"," "),label_size = 32,label_x = c(0,0.05,0,0,0.05,0),
                 rel_widths = c(2,2,0.5),align = "v",axis = "l",rel_heights = c(1.2,1)),
       filename = "procrustes_NMDS_boxplot_LUX_date_color_p.svg",path = here::here("Manuscript","Lux_Arbor_comm_figs"),width = 15,height =12)
#NOT INCLUDED IN REPOSITORY


