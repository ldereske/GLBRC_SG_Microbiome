###Code for exploring mock samples
here::i_am("R_code/OTU_bacterial_fungal_community_analyses_20210623.R")
library(here)
library(phyloseq)
library(plyr); library(dplyr)
library(reshape2)
library(tidyr)
library(seqinr)
library(lme4)
library(car)
library(ggplot2)
library(vegan)
library(otuSummary)
library(stringr)
library(microDecon)
library(cowplot)
library(lmerTest)
library(emmeans)
library(gdm)
pa=function(x)(ifelse(x>0,1,0))
p5=function(x)(ifelse(x>4,x,0))
"%w/o%" <- function(x,y)!('%in%'(x,y))

#####START: Pre-processing Bacterial community####

#Read in the fasta file

rep_set.GLBRC018_bact_OTU<- read.fasta(here::here("Bact_HPCC_out","rep_set_combined_soil_roots_merged_16S_GLBRC018_otus.fasta"), as.string = TRUE, set.attributes = FALSE)
head(rep_set.GLBRC018_bact_OTU)
length(rep_set.GLBRC018_bact_OTU)
#48842

#Load seq mapping
GLBRC018_bact_map=read.csv(here::here("Bact_HPCC_out","Mapping_bact_GLBRC_2018.csv"), header = T)
summary(GLBRC018_bact_map)
nrow(GLBRC018_bact_map)
#1723
head(GLBRC018_bact_map)
#Let's combine the map with the metadata

#Metadata not included in this GIT directory. Let me know if you need the raw files
MMPRNT_16s.metadata_mapp_raw=read.table(here::here("D:/MMPRNT_16S_016-018/ABIOTIC_data/MMPRNT.metadata_mapp_raw_comp_OUTLIER_filter.txt"))

MMPRNT.metadata_mapp_trunc=MMPRNT_16s.metadata_mapp_raw[,c("sampleID_long","collectionDate","siteID","plotType","FertStatus","plotRep","year","dis_Mid","UTM_Lat_Cord","UTM_Lon_Cord")]

MMPRNT.metadata_mapp_trunc[1:10,1:10]
nrow(MMPRNT.metadata_mapp_trunc)
#2700
#colnames(GLBRC018_bact_map)[colnames(GLBRC018_bact_map)=="field_sampleID"]="sampleID_long"

GLBRC018_bact_map_metadata=merge(GLBRC018_bact_map,MMPRNT.metadata_mapp_trunc, by="sampleID_long", all.x = T)
head(GLBRC018_bact_map_metadata)
summary(GLBRC018_bact_map_metadata)
colnames(GLBRC018_bact_map_metadata)
row.names(GLBRC018_bact_map_metadata)=GLBRC018_bact_map_metadata$SampleID 
GLBRC018_bact_map_metadata[1:10,1:10]

#OTU table
otu_GLBRC018_bact=otu_table(read.table(here::here("Bact_HPCC_out","combined_soil_roots_merged_16S_GLBRC018_OTU_table.txt"),sep = "\t", header = T, row.names = 1), taxa_are_rows = T)
otu_GLBRC018_bact[1:10,1:10]


#Load Taxa table
#####CONSTAX IS TROUBLING ME####
#I am going to use CONSTAX for taxonomy since it goes deeper and catches some crap sequences
GLBRC018_bact_raw_CONSTAX_UNITE8.2_all = read.delim(here::here("Bact_HPCC_out","constax_V2_SILVA132_1","combined_taxonomy.txt"),
                                                                     sep = "\t",header = T,fill=T, row.names = 1)
head(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all)
#48842

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all[GLBRC018_bact_raw_CONSTAX_UNITE8.2_all==""|
                                                          is.na(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all)]="Unknown"
head(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all)
unique(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_1_Consensus)

#let's look into the strength of classifications for domain

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Dom_unkn=ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_1_RDP=="Unknown"|
                                                           GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_1_BLAST=="Unknown"|
                                                           GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_1_SINTAX=="Unknown", "TRUE","FALSE")


nrow(subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all,Dom_unkn==TRUE))
#5444

#How many TAXA have a classification at Domain across the three measures but differ in their phyla classification

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all,Dom_unkn==FALSE)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW)
#43398

Rank_2_RDP
GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Dom_diff=ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Rank_1_RDP!=
                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Rank_1_BLAST|
                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Rank_1_RDP!=
                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Rank_1_SINTAX|
                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Rank_1_BLAST!=
                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW$Rank_1_SINTAX, "TRUE","FALSE")




nrow(subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW,Dom_diff==TRUE))
#219

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNW,Dom_diff==TRUE)



sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_RDP==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_BLAST,0,1))
#219
sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_BLAST,0,1))
#0
sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_RDP==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_SINTAX,0,1))
#219

sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_SINTAX!=
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_BLAST&
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif$Rank_1_RDP,1,0))
#0

unique(with(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif,interaction(Rank_1_RDP,Rank_1_SINTAX,Rank_1_BLAST)))
#Bacteria_1.Archaea_1.Archaea_1    Eukaryota_1.Bacteria_1.Bacteria_1

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif_df=data.frame(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif) %>% group_by(Rank_1_SINTAX,Rank_1_BLAST,Rank_1_RDP)%>%summarise_at(vars(Rank_2_RDP),~n())

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif_df[order(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif_df$Rank_2_RDP,decreasing = T),]
#Rank_1_SINTAX Rank_1_BLAST Rank_1_RDP  Rank_2_RDP
#1 Archaea_1     Archaea_1    Bacteria_1         207
#2 Bacteria_1    Bacteria_1   Eukaryota_1         12


write.csv(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif, here::here("Bact_HPCC_out","constax_V2_SILVA132_1","RDP_bug_combined_taxonomy.csv"))


#I am going to output the rep seqs for these weird OTUs
length(rep_set.GLBRC018_bact_OTU)
#8039

rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug=rep_set.GLBRC018_bact_OTU[names(rep_set.GLBRC018_bact_OTU) %in% 
                                                                      c(row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_DOM_KNw_dif))]
head(rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug)
length(rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug)
#219


write.fasta(sequences =rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug, names = names(rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug), 
            file.out =here::here("Bact_HPCC_out","constax_V2_SILVA132_1","rep_set.GLBRC018_bact_OTU_cons_RDP_bug.fna"))


#let's look into the strength of classifications

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phyla_unkn=ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_2_RDP=="Unknown"|
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_2_BLAST=="Unknown"|
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Rank_2_SINTAX=="Unknown", "TRUE","FALSE")


nrow(subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all,Phyla_unkn==TRUE))
#27789

#How many TAXA have a classification at Phyla across the three measures but differ in their phyla classification

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all,Phyla_unkn==FALSE)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW)
#21053

Rank_2_RDP
GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Phyla_diff=ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Rank_2_RDP!=
                                                                                   GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Rank_2_BLAST|
                                                                                   GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Rank_2_RDP!=
                                                                                   GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Rank_2_SINTAX|
                                                                                   GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Rank_2_BLAST!=
                                                                                   GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW$Rank_2_SINTAX, "TRUE","FALSE")




nrow(subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW,Phyla_diff==TRUE))
#12766

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW,Phyla_diff==TRUE)



sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_RDP==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_BLAST,0,1))
#12759
sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_BLAST,0,1))
#337
sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_RDP==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_SINTAX,0,1))
#12762

sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_SINTAX!=
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_BLAST&
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Rank_2_RDP,1,0))
#4

unique(with(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif,interaction(Rank_2_SINTAX,Rank_2_RDP,Rank_2_BLAST)))


GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_df=data.frame(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif) %>% group_by(Rank_2_RDP,Rank_2_SINTAX,Rank_2_BLAST)%>%summarise_at(vars(Rank_1_RDP),~n())

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_df[order(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_df$Rank_1_RDP,decreasing = T),]

#Rank_2_SINTAX                  Rank_2_BLAST                   Rank_2_RDP       Rank_1_RDP
#<chr>                          <chr>                          <chr>                 <int>
#1 Planctomycetota_1              Planctomycetota_1              Proteobacteria_1       3370
#2 Myxococcota_1                  Myxococcota_1                  Proteobacteria_1       2124
#3 Verrucomicrobiota_1            Verrucomicrobiota_1            Proteobacteria_1       1659
#4 Bdellovibrionota_1             Bdellovibrionota_1             Proteobacteria_1       1124
#5 Bacteroidota_1                 Bacteroidota_1                 Proteobacteria_1        997
#6 Acidobacteriota_1              Acidobacteriota_1              Proteobacteria_1        961
#7 Gemmatimonadota_1              Gemmatimonadota_1              Proteobacteria_1        524
#8 Patescibacteria_1              Patescibacteria_1              Proteobacteria_1        334
#9 Dependentiae_1                 Dependentiae_1                 Proteobacteria_1        256
#10 Sar324_clade(marine_group_b)_1 SAR324_clade(Marine_group_B)_1 Proteobacteria_1        234


#Let's output the sintax buggy results 
GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif, Phylum_SINTAX!=Phylum_RDP&Phylum_SINTAX!=Phylum_BLAST)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX)
#13
write.csv(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX, here::here("Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_combined_taxonomy.csv"))


#Let's load in the consensus taxonomy to out output the strange OTUs 

GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus = read.delim(here::here("Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","constax_taxonomy.txt"),
                                                    sep = "\t",header = T,fill=T, row.names = 1)
head(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus)
#8039


GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus_SINTAX=GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus[row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus)%in%row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX),]
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus_SINTAX)
#13


write.csv(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus_SINTAX, here::here("Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_constax_taxonomy.csv"))

#now let's load in the actual Syntax database

GLBRC018_bact_raw_SINTAX_UNITE8.2 = read.delim(here::here("Bact_HPCC_out","taxonomy_combined_merged_ITS1_2_GLBRC_otus.SINTAX"),
                                                          sep = "\t",header = F,fill=T, row.names = 1)

head(GLBRC018_bact_raw_SINTAX_UNITE8.2)
nrow(GLBRC018_bact_raw_SINTAX_UNITE8.2)
#8039


GLBRC018_bact_raw_SINTAX_UNITE8.2_SINTAX=GLBRC018_bact_raw_SINTAX_UNITE8.2[row.names(GLBRC018_bact_raw_SINTAX_UNITE8.2)%in%row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX),]
nrow(GLBRC018_bact_raw_SINTAX_UNITE8.2_SINTAX)
#13

write.csv(GLBRC018_bact_raw_SINTAX_UNITE8.2_SINTAX,here::here("SG_Microbiome/Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_taxonomy_combined_merged_ITS1_2_GLBRC_otus.csv"))



#total number of phyla classification that differ with unknowns included

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phyla_diff=ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phylum_RDP!=
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phylum_BLAST|
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phylum_RDP!=
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phylum_SINTAX|
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phylum_BLAST!=
                                                                               GLBRC018_bact_raw_CONSTAX_UNITE8.2_all$Phylum_SINTAX, "TRUE","FALSE")
nrow(subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all,Phyla_diff==TRUE))
#2276

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all,Phyla_diff==TRUE)



sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_RDP==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_BLAST,0,1))
#1922
sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_BLAST,0,1))
#2001
sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_RDP==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX,0,1))
#648

sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX!=
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_BLAST&
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_RDP,1,0))
#1628


#Let's look at the times when RDP is UNKNOWN

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif, 
                                                        Phylum_RDP=="Unknown"&Phylum_BLAST!="Unknown"&
                                                          Phylum_SINTAX!="Unknown")
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B)
#108

sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B$Phylum_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B$Phylum_BLAST,0,1))
#4

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B_SINTAX=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B,
                                                               Phylum_SINTAX!=Phylum_BLAST)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B_SINTAX)
#4


#Let's look at the times when BLAST is UNKNOWN

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif, 
                                                        Phylum_RDP!="Unknown"&Phylum_BLAST=="Unknown"&
                                                          Phylum_SINTAX!="Unknown")
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R)
#329

sum(ifelse(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R$Phylum_SINTAX==
             GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R$Phylum_RDP,0,1))
#10

GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R_SINTAX=subset(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R,
                                                               Phylum_SINTAX!=Phylum_RDP)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R_SINTAX)
#10
GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX=rbind(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_B_SINTAX,GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_S_R_SINTAX)
write.csv(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX, here::here("SG_Microbiome/Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_UNKNOWN_combined_taxonomy.csv"))
head(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX)

#Let's load in the consensus taxonomy to out output the strange OTUs 


head(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus)
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus)
#8039


GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus_SINTAX_UNK=GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus[row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus)%in%row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX),]
nrow(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus_SINTAX_UNK)
#14


write.csv(GLBRC018_bact_raw_CONSTAX_UNITE8.2_concensus_SINTAX_UNK, here::here("SG_Microbiome/Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_UNK_bug_constax_taxonomy.csv"))

#now let's load in the actual Syntax database



head(GLBRC018_bact_raw_SINTAX_UNITE8.2)
nrow(GLBRC018_bact_raw_SINTAX_UNITE8.2)
#8039


GLBRC018_bact_raw_SINTAX_UNITE8.2_SINTAX_UNK=GLBRC018_bact_raw_SINTAX_UNITE8.2[row.names(GLBRC018_bact_raw_SINTAX_UNITE8.2)%in%row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX),]
nrow(GLBRC018_bact_raw_SINTAX_UNITE8.2_SINTAX_UNK)
#14

write.csv(GLBRC018_bact_raw_SINTAX_UNITE8.2_SINTAX_UNK, here::here("SG_Microbiome/Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_UNK_bug_taxonomy_combined_merged_ITS1_2_GLBRC_otus.csv"))


#I am going to output the rep seqs for these weird OTUs
length(rep_set.GLBRC018_bact_OTU)
#8039

rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug=rep_set.GLBRC018_bact_OTU[names(rep_set.GLBRC018_bact_OTU) %in% 
                                                                      c(row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX),
                                                                        row.names(GLBRC018_bact_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX))]
head(rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug)
length(rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug)
#27


write.fasta(sequences =rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug, names = names(rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug), 
            file.out =here::here("SG_Microbiome/Bact_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","rep_set.GLBRC018_bact_OTU_cons_SINTAX_bug.fna"))

#####CONSTAX IS TROUBLING ME###

#I am going to use the sintax classification for now 


GLBRC018_bact_raw_SINTAX_UNITE8.2 = read.delim(here::here("Bact_HPCC_out","taxonomy_combined_merged_ITS1_2_GLBRC_otus.SINTAX"),
                                               sep = "\t",header = F,fill=T)


nrow(GLBRC018_bact_raw_SINTAX_SILVA123)
#48842
head(GLBRC018_bact_raw_SINTAX_SILVA123)
GLBRC018_bact_raw_SINTAX_SILVA123_80c=GLBRC018_bact_raw_SINTAX_SILVA123[,c(1,4)]
head(GLBRC018_bact_raw_SINTAX_SILVA123_80c)

GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep=GLBRC018_bact_raw_SINTAX_SILVA123_80c %>% separate(V4, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ",")
row.names(GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep)=GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep$V1
GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep$V1=NULL
GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep[is.na(GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep)] <- "Unknown"
GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep_mat=as.matrix(GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep)
head(GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep_mat)
unique(data.frame(GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep_mat)$Domain)
TAXA_80C_UNITE8.2_SILVA123=tax_table(GLBRC018_bact_raw_SINTAX_SILVA123_80c_sep_mat)
head(TAXA_80C_UNITE8.2_SILVA123)




GLBRC018_OTU_bact_phyl=phyloseq(otu_GLBRC018_bact,sample_data(GLBRC018_bact_map_metadata),TAXA_80C_UNITE8.2_SILVA123)
nsamples(GLBRC018_OTU_bact_phyl)
#1703
ntaxa(GLBRC018_OTU_bact_phyl)
#48842
sum(taxa_sums(GLBRC018_OTU_bact_phyl))
#55627577
mean(sample_sums(GLBRC018_OTU_bact_phyl))
#32664.46
min(sample_sums(GLBRC018_OTU_bact_phyl))
#1
max(sample_sums(GLBRC018_OTU_bact_phyl))
#179058
sort(sample_sums(GLBRC018_OTU_bact_phyl))[1:119]


#Amp396     Amp721     Amp722     Amp790      Amp16      Amp40      Amp52   CTRL12B5      Amp28      Amp76 MMPRNT1365  MERDSSG42     Amp481     Amp763 MMPRNT1932     Amp784     Amp541 
#1          1          1          3         21         22         33         44         45        161        162        236        240        378        433        438        499 
#Amp654     Amp610     Amp522     Amp554     Amp761     Amp762     Amp755     Amp756     Amp683   CTRL6F12     Amp759     Amp493     Amp667     Amp600     Amp760 MMPRNT1644     Amp565 
#509        536        582        630        644        673        731        813        862        881        948        960        972        975       1003       1017       1040 
#Amp438     Amp534     Amp648     Amp498     Amp402     Amp614     Amp563     Amp467     Amp367     Amp424      Amp95 MMPRNT1371     Amp397     Amp268     Amp234      Amp12     Amp718 
#1052       1053       1070       1072       1102       1133       1267       1341       1341       1466       1628       1633       1715       1717       1779       1860       1973 
#Amp267      Amp80      Amp46     Amp364     Amp458 MMPRNT2007     Amp743     Amp198     Amp258     Amp191      Amp24     Amp246     Amp210     Amp365     Amp287     Amp270     Amp222 
#2006       2033       2100       2249       2254       2319       2407       2553       2565       2577       2609       2645       2773       2873       2935       3105       3124 
#Amp277     Amp275     Amp288     Amp213     Amp176      Amp74     Amp170     Amp282     Amp142      Amp72 MMPRNT1724     Amp111     Amp161     Amp203     Amp345     Amp227      Amp68 
#3252       3267       3291       3371       3730       3738       3812       3851       4027       4102       4435       4474       4726       4726       4953       5144       5189 
#Amp294      Amp96      Amp92     Amp572    CTRL8C8     Amp368     Amp560     Amp524      Amp42     Amp298       Amp4     Amp548 MMPRNT1916     Amp291     Amp307 Marshall44 MMPRNT1456 
#5309       5364       5652       5726       6248       6424       6599       6608       7289       7519       7736       8090       8107       8577       8910       8923       8977 
#MMPRNT1919 MMPRNT1921 MMPRNT1277      Amp88 MMPRNT1920      Amp64  MERDSSG85 MMPRNT1913 MMPRNT1528 MMPRNT1918     Amp782 MMPRNT1741 MMPRNT1739     Amp778 MERDSSG103 MMPRNT1455      Amp50 
#9136       9282       9920       9951      10276      10386      10560      10880      11104      11238      11292      11318      11411      11418      11534      11754      11841 


sort(sample_sums(GLBRC018_OTU_bact_phyl),decreasing =T)[1:119]


#Amp289 Amp120 Amp214 Amp319 Amp201 Amp302 Amp325 Amp308 Amp226 Amp323 Amp329 Amp333 Amp297 Amp334 Amp309 Amp327 Amp320 Amp335 Amp332 Amp331 Amp314 Amp303 Amp290 Amp284 Amp311 Amp202 Amp299 Amp296 
#179058 153643 144787 142867 132676 131724 130303 130040 129948 129822 126375 126282 123955 123584 122367 118117 117697 117319 115060 114589 114180 114171 113891 112729 112476 112369 111126 110297 
#Amp301 Amp292 Amp326 Amp324 Amp306 Amp204 Amp336 Amp225 Amp720 Amp295 Amp127 Amp321 Amp108 Amp552 Amp250 Amp107 Amp328 Amp344 Amp313 Amp322 Amp310 Amp156 Amp356 Amp305 Amp716 Amp112 Amp300 Amp362 
#110045 109414 109283 107784 105969 105355 104899 104083 103384 103357 102877  99497  99016  98532  98512  98298  98001  97940  97440  96973  96779  96740  96200  96192  95604  93616  92656  92167 
#Amp216 Amp273 Amp251 Amp704 Amp315 Amp285 Amp330 Amp689 Amp286 Amp357 Amp168 Amp678 Amp228 Amp262 Amp209 Amp266 Amp237 Amp316 Amp318 Amp264 Amp304 Amp353 Amp238 Amp245 Amp252 Amp117 Amp337 Amp205 
#91181  91084  90955  89904  89684  89142  88830  86545  86278  85758  85747  85491  85275  85167  84933  84718  84502  84439  84321  84256  84216  84211  83025  81650  81496  81285  81219  80839 
#Amp261 Amp564 Amp317 Amp502 Amp348 Amp706 Amp694 Amp272 Amp248 Amp121 Amp358 Amp206 Amp113 Amp276 Amp215 Amp528 Amp346 Amp745 Amp274 Amp491 Amp241 Amp526 Amp349 Amp682 Amp566 Amp673 Amp746 Amp192 
#80721  80635  80614  80145  79585  79507  79402  79349  79207  78897  78283  78275  78129  78000  77875  77736  76632  76597  76505  76351  76192  75867  75607  74491  74487  74091  73933  73721 
#Amp162 Amp254 Amp347 Amp568 Amp278 Amp575 Amp114 
#72984  72801  72729  72478  72181  72005  71380 


#There is a difference i the depth of sequencing between labs

#####Mock community before removing bact######

GLBRC018_OTU_bact.mock=subset_samples(GLBRC018_OTU_bact_phyl,seq_samp_type=="mock")
nsamples(GLBRC018_OTU_bact.mock)
#17
GLBRC018_OTU_bact.mock=prune_taxa(taxa_sums(GLBRC018_OTU_bact.mock) > 0, GLBRC018_OTU_bact.mock)
ntaxa(GLBRC018_OTU_bact.mock)
#202


sum(otu_table(GLBRC018_OTU_bact.mock))
#1145163

max(sample_sums(GLBRC018_OTU_bact.mock))
#153643
min(sample_sums(GLBRC018_OTU_bact.mock))
#19353

sort(sample_sums(GLBRC018_OTU_bact.mock))
#ZymoMockC ZymoMockA ZymoMockB    Amp661    Amp525    Amp589     Amp31    Amp429     Amp23    Amp280    Amp395    Amp491    Amp689    Amp127    Amp334    Amp201    Amp120 
#####19353     19516     20124     34188     39111     41661     41757     43351     46198     60380     67092     76351     86545    102877    123584    132676    153643 


sort(taxa_sums(GLBRC018_OTU_bact.mock),decreasing = T)

rep_set.GLBRC018_bact_OTU_mock=rep_set.GLBRC018_bact_OTU[names(rep_set.GLBRC018_bact_OTU) %in% taxa_names(GLBRC018_OTU_bact.mock)]
head(rep_set.GLBRC018_bact_OTU_mock)
length(rep_set.GLBRC018_bact_OTU_mock)
#202


write.fasta(sequences =rep_set.GLBRC018_bact_OTU_mock, names = names(rep_set.GLBRC018_bact_OTU_mock), 
            file.out =here::here("Bacterial_mock","rep_set.GLBRC018_bact_OTU_mock.fna"))

#####Mock USEARCH_global####

#Code run in Usearch

#cd HardDrive/GLBRC-TeamMicrobiome/SG_Microbiome/Bacterial_mock/

#The program needs for the codons to be uppercase
#awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' rep_set.GLBRC018_bact_OTU_mock.fna > CAP_rep_set.GLBRC018_bact_OTU_mock.fna


#~/HardDrive/Sciencey_Program/Old_usearch_v/usearch10.0.240_win32.exe -usearch_global CAP_rep_set.GLBRC018_bact_OTU_mock.fna -db zymo_mock.fasta -id 0.97 -strand both -maxaccepts 20 -maxrejects 50 -matched match_OTU_GLBRC018_bact_mock.fa -notmatched NOT_match_OTU_GLBRC018_bact_mock.fa -userout match_OTU_GLBRC018_bact_mock.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir


#00:01 10Mb    100.0% Searching, 0.7% matched


#~/HardDrive/Sciencey_Program/Old_usearch_v/usearch10.0.240_win32.exe -usearch_global CAP_rep_set.GLBRC018_bact_OTU_mock.fna -db zymo_mock.fasta -id 0.8 -strand both -maxaccepts 20 -maxrejects 50 -matched match_OTU_GLBRC018_bact_mock80.fa -notmatched NOT_match_OTU_GLBRC018_bact_mock80.fa -userout match_OTU_GLBRC018_bact_mock80.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir

#00:02 16Mb    100.0% Searching, 50.7% matched

#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")

matchs_mocks_OTU_bact=read.table(here::here("Bacterial_mock","match_OTU_GLBRC018_bact_mock.txt"),header = F)
head(matchs_mocks_OTU_bact)
colnames(matchs_mocks_OTU_bact)=user_out_col
nrow(matchs_mocks_OTU_bact)
#60
length(unique(matchs_mocks_OTU_bact$target))
#46

GLBRC018_OTU_bact.mock_match<-prune_taxa(as.character(unique(matchs_mocks_OTU_bact$query)),
                                     GLBRC018_OTU_bact.mock)
sum(taxa_sums(GLBRC018_OTU_bact.mock_match))/sum(taxa_sums(GLBRC018_OTU_bact.mock))
#0.9981522 match to mock community sequences

#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")

matchs_mocks_OTU_bact_80=read.table(here::here("Bacterial_mock","match_OTU_GLBRC018_bact_mock80.txt"),header = F)
head(matchs_mocks_OTU_bact_80)
colnames(matchs_mocks_OTU_bact_80)=user_out_col
nrow(matchs_mocks_OTU_bact_80)
#12237
length(unique(matchs_mocks_OTU_bact_80$target))
#46


GLBRC018_OTU_bact.mock_match80<-prune_taxa(as.character(unique(matchs_mocks_OTU_bact_80$query)),
                                             GLBRC018_OTU_bact.mock)
ntaxa(GLBRC018_OTU_bact.mock_match80)
sum(taxa_sums(GLBRC018_OTU_bact.mock_match80))/sum(taxa_sums(GLBRC018_OTU_bact.mock))
#0.999759

mock_comp_OTU=data.frame("match_sum80"=sample_sums(GLBRC018_OTU_bact.mock_match80),
                         "match_sum97"=sample_sums(GLBRC018_OTU_bact.mock_match),
                         "tot_sum"=sample_sums(GLBRC018_OTU_bact.mock), 
                         estimate_richness(GLBRC018_OTU_bact.mock_match,measures = "Observed"),
                         estimate_richness(GLBRC018_OTU_bact.mock_match80,measures = "Observed"),
                         estimate_richness(GLBRC018_OTU_bact.mock,measures = "Observed"))
colnames(mock_comp_OTU)[4:6]=c("match_rich97","match_rich80","tot_rich")

mock_comp_OTU$prop_read_match97=mock_comp_OTU$match_sum97/mock_comp_OTU$tot_sum
mock_comp_OTU$prop_read_match80=mock_comp_OTU$match_sum80/mock_comp_OTU$tot_sum

mock_comp_OTU$prop_taxa_match97=mock_comp_OTU$match_rich97/mock_comp_OTU$tot_rich
mock_comp_OTU$prop_taxa_match80=mock_comp_OTU$match_rich80/mock_comp_OTU$tot_rich
mock_comp_OTU$sampleID=row.names(mock_comp_OTU)

mock_comp_OTU_m=melt(mock_comp_OTU[,c("sampleID","prop_read_match97","prop_read_match80","prop_taxa_match97","prop_taxa_match80")])

var_order= c("prop_read_match97","prop_taxa_match97","prop_read_match80","prop_taxa_match80")
otu_var_names=c("Reads matching\nat 97%","OTU matching\nat 97%","Reads matching\nat 80%","OTU matching\nat 80%")

ggplot(mock_comp_OTU_m, aes(x=factor(variable, levels = var_order),y=value, group=sampleID, color=sampleID))+geom_point(size=4)+
  scale_x_discrete(name=NULL,label=otu_var_names)+ylab("Precision\n(proportion of sample)")+theme_bw()+theme(axis.text = element_text(size=16),
                                                                                                             axis.title = element_text(size=20))
mock_comp_OTU_m$pipeline=rep("OTU")







#How frequent are these taxa in the other samples?
#97% match to mock community sequences

#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")
matchs_mocks_OTU=read.table(here::here("Bacterial_mock","match_OTU_GLBRC018_bact_mock.txt"),header = F)
head(matchs_mocks_OTU)
colnames(matchs_mocks_OTU)=user_out_col
nrow(matchs_mocks_OTU)
#60
length(unique(matchs_mocks_OTU$target))
#46

GLBRC018_OTU_bact_phyl_match<-prune_taxa(as.character(unique(matchs_mocks_OTU$query)),
                                     GLBRC018_OTU_bact_phyl)



#Remove the mock samples and samples with no reads
GLBRC018_OTU_bact.NO_mock_match=subset_samples(GLBRC018_OTU_bact_phyl_match,seq_samp_type!="mock")
nsamples(GLBRC018_OTU_bact.NO_mock_match)
#1686
GLBRC018_OTU_bact.NO_mock_match=prune_samples(sample_sums(GLBRC018_OTU_bact.NO_mock_match) > 0, GLBRC018_OTU_bact.NO_mock_match)
nsamples(GLBRC018_OTU_bact.NO_mock_match)
#1335
ntaxa(GLBRC018_OTU_bact.NO_mock_match)
#11

sort((sample_sums(GLBRC018_OTU_bact.NO_mock_match)),decreasing = T)

GLBRC018_OTU_bact.NO_mock_match_SUM_DF=data.frame(sample_sums(GLBRC018_OTU_bact.NO_mock_match))
colnames(GLBRC018_OTU_bact.NO_mock_match_SUM_DF)="mock_reads"

hist(sample_sums(GLBRC018_OTU_bact.NO_mock_match))

ggplot(GLBRC018_OTU_bact.NO_mock_match_SUM_DF,aes(x=mock_reads))+geom_histogram()+
  scale_x_continuous(name = "Numbers of reads matching mock per sample")+ylab("Frequency")+
  theme_classic()+theme(axis.title = element_text(size = 20),axis.text = element_text(size = 16))



GLBRC018_OTU_bact_phyl_match_SUM_DF=data.frame(sample_sums(GLBRC018_OTU_bact_phyl))
colnames(GLBRC018_OTU_bact_phyl_match_SUM_DF)="tot_reads"


GLBRC018_OTU_bact_phyl_mock_match_SUM_DF=merge(GLBRC018_OTU_bact_phyl_match_SUM_DF,GLBRC018_OTU_bact.NO_mock_match_SUM_DF,by="row.names")
summary(GLBRC018_OTU_bact_phyl_mock_match_SUM_DF)

GLBRC018_OTU_bact_phyl_mock_match_SUM_DF$prop_mock=GLBRC018_OTU_bact_phyl_mock_match_SUM_DF$mock_reads/GLBRC018_OTU_bact_phyl_mock_match_SUM_DF$tot_reads


row.names(GLBRC018_OTU_bact_phyl_mock_match_SUM_DF)=GLBRC018_OTU_bact_phyl_mock_match_SUM_DF$Row.names

GLBRC018_OTU_bact_phyl_mock_match_SUM_DF[order(GLBRC018_OTU_bact_phyl_mock_match_SUM_DF$prop_mock,decreasing = T),][1:100,]


#Here are the three weird samples
#Amp279 Amp531 Amp528
sample_data(GLBRC018_OTU_bact_phyl)
sample_sums(subset_samples(GLBRC018_OTU_bact_phyl, SampleID=="Amp279"|SampleID=="Amp531"|SampleID=="Amp528"))








#####Mock based taxa control#####

#Let's try to figure out how to eliminate the Contamination and Crap OTUs

ntaxa(GLBRC018_OTU_bact.mock)
#202
sum(otu_table(GLBRC018_OTU_bact.mock))
#1145163

max(sample_sums(GLBRC018_OTU_bact.mock))
#153643
min(sample_sums(GLBRC018_OTU_bact.mock))
#19353

sort(sample_sums(GLBRC018_OTU_bact.mock))
#ZymoMockC ZymoMockA ZymoMockB    Amp661    Amp589     Amp31    Amp429     Amp23    Amp280    Amp395    Amp526    Amp491    Amp689    Amp127 
#19353     19516     20124     34188     41661     41757     43351     46198     60380     67092     75867     76351     86545    102877 
#Amp334    Amp201    Amp120 
#123584    132676    153643 





#Let's remove taxa that only occur in 2 sample 

GLBRC018_OTU_bact.mock_pa=transform_sample_counts(GLBRC018_OTU_bact.mock,pa)
max(otu_table(GLBRC018_OTU_bact.mock_pa))

GLBRC018_OTU_bact.mock_v2=prune_taxa(taxa_sums(GLBRC018_OTU_bact.mock_pa) > 2, GLBRC018_OTU_bact.mock)
ntaxa(GLBRC018_OTU_bact.mock_v2)
#39
sum(otu_table(GLBRC018_OTU_bact.mock_v2))
#1144928

max(sample_sums(GLBRC018_OTU_bact.mock_v2))
#153633
min(sample_sums(GLBRC018_OTU_bact.mock_v2))
#19343
otu_table(GLBRC018_OTU_bact.mock_v2)
tax_table(GLBRC018_OTU_bact.mock_v2)


#Let's ordinate them
GLBRC018_OTU_bact.mock_v2_ord=ordinate(GLBRC018_OTU_bact.mock_v2,method = "NMDS")
#*** Solution reached

head(sample_data(GLBRC018_OTU_bact.mock_v2))

plot_ordination(GLBRC018_OTU_bact.mock_v2,GLBRC018_OTU_bact.mock_v2_ord)+geom_point(aes(color = DNA_plate), size=3)+geom_text(aes(label = SampleID))



adonis(distance(GLBRC018_OTU_bact.mock_v2,method = "bray")~sample_data(GLBRC018_OTU_bact.mock_v2)$DNA_plate,permutations = 9999)
#Error in G * t(hat) : non-conformable arrays
#Pretty bad

#Effects of library

plot_ordination(GLBRC018_OTU_bact.mock_v2,GLBRC018_OTU_bact.mock_v2_ord)+geom_point(aes(color = Run ), size=3)+geom_text(aes(label = SampleID))

adonis(distance(GLBRC018_OTU_bact.mock_v2,method = "bray")~sample_data(GLBRC018_OTU_bact.mock_v2)$Run,permutations = 9999)
#sample_data(GLBRC018_OTU_bact.mock_v2)$Run  4   0.77790 0.194475  4.2243 0.58473 0.0038 **
#not super strong


#Rarefying the community

GLBRC018_OTU_bact.mock_v2_rar=rarefy_even_depth(GLBRC018_OTU_bact.mock_v2,replace=F,sample.size=19000,rngseed=99)
#3OTUs were removed because they are no longer 
#present in any sample after random subsampling
ntaxa(GLBRC018_OTU_bact.mock_v2_rar)
#36
sum(otu_table(GLBRC018_OTU_bact.mock_v2_rar))
#323000

max(sample_sums(GLBRC018_OTU_bact.mock_v2_rar))
#19000
min(sample_sums(GLBRC018_OTU_bact.mock_v2_rar))
#19000
otu_table(GLBRC018_OTU_bact.mock_v2_rar)
tax_table(GLBRC018_OTU_bact.mock_v2_rar)


#Lets ordinate them
GLBRC018_OTU_bact.mock_v2_rar_ord=ordinate(GLBRC018_OTU_bact.mock_v2_rar,method = "NMDS")
#*** Solution reached
#0.0561177

#
#Effects of DNA plate 

plot_ordination(GLBRC018_OTU_bact.mock_v2_rar,GLBRC018_OTU_bact.mock_v2_rar_ord)+geom_point(aes(color = DNA_plate), size=3)+geom_text(aes(label = SampleID))

adonis(distance(GLBRC018_OTU_bact.mock_v2_rar,method = "bray")~sample_data(GLBRC018_OTU_bact.mock_v2_rar)$DNA_plate,permutations = 9999)
#Error in G * t(hat) : non-conformable arrays
#Not great still

#Effects of library

plot_ordination(GLBRC018_OTU_bact.mock_v2_rar,GLBRC018_OTU_bact.mock_v2_rar_ord)+geom_point(aes(color = Run), size=3)+geom_text(aes(label = SampleID))

adonis(distance(GLBRC018_OTU_bact.mock_v2_rar,method = "bray")~sample_data(GLBRC018_OTU_bact.mock_v2_rar)$Run,permutations = 9999)
#Error in `colnames<-`(`*tmp*`, value = colnames(lhs)) : 
#attempt to set 'colnames' on an object with less than two dimensions
#not super strong

#####MMPRNT Only ######
head(sample_data(GLBRC018_OTU_bact_phyl))
unique(sample_data(GLBRC018_OTU_bact_phyl)$Project)
nsamples(GLBRC018_OTU_bact_phyl)
#1703
GLBRC018_OTU_bact_MMPRNT_mock=subset_samples(GLBRC018_OTU_bact_phyl,Project=="MMPRNT"|
                                          seq_samp_type=="mock")

GLBRC018_OTU_bact_MMPRNT_mock=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock) > 0, GLBRC018_OTU_bact_MMPRNT_mock)
nsamples(GLBRC018_OTU_bact_MMPRNT_mock)
#1212

ntaxa(GLBRC018_OTU_bact_MMPRNT_mock)
# 47034
sum(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#33652520
mean(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#27766.11
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#1
max(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#153643

unique(sample_data(GLBRC018_OTU_bact_MMPRNT_mock)$Sample_or_Control)


#Let's remove the controls and composites from the dataset

GLBRC018_OTU_bact_MMPRNT_mock=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock,Root_soil=="Root"|Root_soil=="Soil"|seq_samp_type=="mock")

GLBRC018_OTU_bact_MMPRNT_mock=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock) > 0, GLBRC018_OTU_bact_MMPRNT_mock)
nsamples(GLBRC018_OTU_bact_MMPRNT_mock)
#1152
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock)
# 46948
sum(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#33346041
mean(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#28946.22
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#21
max(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock))
#153643


#Let's choose the best 
unique(sample_data(GLBRC018_OTU_bact_MMPRNT_mock)$rep_pcr) 


GLBRC018_OTU_bact_MMPRNT_mock_re_PCR=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock,rep_pcr!="none")
sample_data(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR)
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR)
#56

unique(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR)$sampleID_long) 


GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_diversity=merge(estimate_richness(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR),
                                                data.frame(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR)),by="row.names")
head(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_diversity)

colnames(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_diversity)[colnames(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_diversity)==
                                                      "sample_sums.GLBRC018_OTU_bact_MMPRNT_mock_re_PCR."]="sample_sum"
GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta=merge(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_diversity,
                                               data.frame(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR)),
                                               by.x = "Row.names",by.y = "row.names")
head(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta)

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta,aes(x=sampleID_long,y=sample_sum,fill=rep_pcr))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")


ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta,aes(x=sampleID_long,y=Observed,fill=rep_pcr))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")


#let's see if the most reads means it is the PCR

GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk=GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta%>%group_by(sampleID_long) %>%
  mutate(PCR_rank = order(order(sample_sum, decreasing=TRUE)))

unique(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk$PCR_rank)

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=sample_sum,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")


#how does the diversity look

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=Observed,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=Shannon ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=InvSimpson ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

#They look okay overall. 
#I want to look at if Mock reads somehow made it into these samples 


#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")
matchs_mocks_OTU=read.table(here::here("Bacterial_mock","match_OTU_GLBRC018_bact_mock.txt"),header = F)
head(matchs_mocks_OTU)
colnames(matchs_mocks_OTU)=user_out_col
nrow(matchs_mocks_OTU)
#60
length(unique(matchs_mocks_OTU$target))
#46

GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_match<-prune_taxa(as.character(unique(matchs_mocks_OTU$query)),
                                                  GLBRC018_OTU_bact_MMPRNT_mock_re_PCR)


hist(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_match))

sort(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_match),decreasing = T)


GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2=merge(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk,
                                                     sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_match),
                                                     by.x = "Row.names",by.y = "row.names")
head(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2)

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2,aes(x=sampleID_long,y=y ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2,aes(x=sampleID_long,y=y ,fill=factor(rep_pcr)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

#Amp279 has the most reads overall and all of them are mocks

subset(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2,Row.names=="Amp279")

#Let's take a look at rePCR of MMPRNT-1539


ggplot(subset(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2, sampleID_long=="MMPRNT-1539"),
       aes(x=sampleID_long,y= sample_sum,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")
ggplot(subset(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2, sampleID_long=="MMPRNT-1539"),
       aes(x=sampleID_long,y= Observed,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

#Amp352 looks good though

ggplot(subset(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2,Row.names!="Amp279"),aes(x=sampleID_long,y=y ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(subset(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2,Row.names!="Amp279"),aes(x=sampleID_long,y=y ,fill=factor(rep_pcr)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

#Let's make a list of samples that we want to keep

prune_PCR_bact_samp=subset(GLBRC018_OTU_bact_MMPRNT_mock_re_PCR_div_meta_rnk_2,PCR_rank!=1)$Row.names
length(prune_PCR_bact_samp)
#31

#Need to replace Amp352 with Amp279

prune_PCR_bact_samp=c(prune_PCR_bact_samp[prune_PCR_bact_samp!="Amp352"],"Amp279")
length(prune_PCR_bact_samp)
#31
length(sample_names(GLBRC018_OTU_bact_MMPRNT_mock))
#1152

final_soil_samples=sample_names(GLBRC018_OTU_bact_MMPRNT_mock)[sample_names(GLBRC018_OTU_bact_MMPRNT_mock) %w/o%prune_PCR_bact_samp]
length(final_soil_samples)
#1121


GLBRC018_OTU_bact_MMPRNT_mock_fin=prune_samples(final_soil_samples,GLBRC018_OTU_bact_MMPRNT_mock)
GLBRC018_OTU_bact_MMPRNT_mock_fin=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin) > 0, GLBRC018_OTU_bact_MMPRNT_mock_fin)
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_fin)
#1121

ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_fin)
# 46926
sum(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin))
#33089683
mean(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin))
#29518

min(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin))
#162
max(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin))
#153643

#Let's check to make sure we have removed the rePCR samples


nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_fin,Root_soil=="Root"&seq_samp_type!="mock"))
#234
length(unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_fin,Root_soil=="Root"&seq_samp_type!="mock"))$sampleID_long))
#234
unique(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_fin)$Root_soil)

nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_fin,Root_soil=="Soil"&seq_samp_type!="mock"))
#870
length(unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_fin,Root_soil=="Soil"&seq_samp_type!="mock"))$sampleID_long))
#870




#####Singletons removal ######

#Let remove the taxon that occur in less than 5 samples 
GLBRC018_OTU_bact_MMPRNT_mock_fin_pa=transform_sample_counts(GLBRC018_OTU_bact_MMPRNT_mock_fin,pa)


GLBRC018_OTU_bact_MMPRNT_mock_fin_fil=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_pa) > 4, GLBRC018_OTU_bact_MMPRNT_mock_fin)
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil)
# 32318
sum(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil))
#33020693

mean(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil))
#29456.46
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil))
#1
max(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil))
#161
sort(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil))[1:40]

#MMPRNT1365 MMPRNT1932 MMPRNT1644 MMPRNT1371 MMPRNT2007 MMPRNT1724      Amp42 MMPRNT1916 MMPRNT1456 MMPRNT1919 MMPRNT1921 MMPRNT1277 MMPRNT1920 
#161        433       1015       1630       2316       4426       7281       8106       8961       9112       9279       9903      10268 
#MMPRNT1913 MMPRNT1528 MMPRNT1918 MMPRNT1741 MMPRNT1739 MMPRNT1455      Amp50 MMPRNT1512     Amp779 MMPRNT1917 MMPRNT1274 MMPRNT1738 MMPRNT1622 
#10862      11016      11230      11281      11388      11706      11837      11911      12567      12707      12823      13138      13600 
#MMPRNT1561 MMPRNT1914 MMPRNT1121 MMPRNT1145 MMPRNT1401 MMPRNT1915 MMPRNT1153 MMPRNT1497 MMPRNT1154 MMPRNT1898 MMPRNT1626      Amp44 MMPRNT1155 
#13615      13685      13710      13845      13910      13938      13989      14104      14137      14165      14205      14331      14573 
#MMPRNT1306 
#14781

#Let's remove the samples that are less than 10000 reads 

GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2=prune_samples(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil)>10000,GLBRC018_OTU_bact_MMPRNT_mock_fin_fil)

GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2) > 0, GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2)
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2)
# 32318
GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2_pa=transform_sample_counts(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2,pa)
GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2_pa) > 4, GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2)
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2)
# 32302
sum(taxa_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))
#32957920
mean(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))
#29718.59
min(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))
#10268
max(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))
#153643
sort(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))[1:10]
#MMPRNT1920 MMPRNT1913 MMPRNT1528 MMPRNT1918 MMPRNT1741 MMPRNT1739 MMPRNT1455      Amp50 MMPRNT1512     Amp779 
#10268      10862      11016      11230      11281      11388      11706      11837      11911      12567





#####Mock community with singletons removed######

GLBRC018_OTU_bact.mock_f2=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2,seq_samp_type=="mock")
nsamples(GLBRC018_OTU_bact.mock_f2)
#17
GLBRC018_OTU_bact.mock_f2=prune_taxa(taxa_sums(GLBRC018_OTU_bact.mock_f2) > 0, GLBRC018_OTU_bact.mock_f2)
ntaxa(GLBRC018_OTU_bact.mock_f2)
#195


sum(otu_table(GLBRC018_OTU_bact.mock_f2))
#1145151

max(sample_sums(GLBRC018_OTU_bact.mock_f2))
#153643
min(sample_sums(GLBRC018_OTU_bact.mock_f2))
#19353

sort(sample_sums(GLBRC018_OTU_bact.mock_f2))
#ZymoMockC ZymoMockA ZymoMockB    Amp661    Amp589     Amp31    Amp429     Amp23    Amp280    Amp395    Amp526    Amp491    Amp689    Amp127 
#19353     19516     20124     34188     41661     41757     43351     46198     60380     67088     75867     76349     86544    102875 
#Amp334    Amp201    Amp120 
#123581    132676    153643 


#Let's ordinate this 
nsamples(GLBRC018_OTU_bact.mock_f2)
#17

GLBRC018_OTU_bact.mock_f2_ord=ordinate(GLBRC018_OTU_bact.mock_f2,method = "NMDS")
#*** Solution reached


plot_ordination(GLBRC018_OTU_bact.mock_f2,GLBRC018_OTU_bact.mock_f2_ord,color = "Run",label = "SampleID")
#not great but the best I can do for now. 

#####Rarefy the community####

GLBRC018_OTU_bact_MMPRNT_mock_rar=rarefy_even_depth(GLBRC018_OTU_bact.mock_f2,replace=F,sample.size=10000,rngseed=99)
#122OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(GLBRC018_OTU_bact_MMPRNT_mock_rar, file = here::here("R_files","GLBRC018_BACT_OTU_MMPRNT_mock_rar.RData"))
load(here::here("R_files","GLBRC018_BACT_OTU_MMPRNT_mock_rar.RData"))



#####Mock community with Rarefy the community######

GLBRC018_OTU_bact_mock_rar=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_rar,seq_samp_type=="mock")
nsamples(GLBRC018_OTU_bact_mock_rar)
#17
GLBRC018_OTU_bact_mock_rar=prune_taxa(taxa_sums(GLBRC018_OTU_bact_mock_rar) > 0, GLBRC018_OTU_bact_mock_rar)
ntaxa(GLBRC018_OTU_bact_mock_rar)
#73


sum(otu_table(GLBRC018_OTU_bact_mock_rar))
#170000

max(sample_sums(GLBRC018_OTU_bact_mock_rar))
#10000
min(sample_sums(GLBRC018_OTU_bact_mock_rar))
#10000




#Let's ordinate this after removing the wonky samples
nsamples(GLBRC018_OTU_bact_mock_rar)
#17


GLBRC018_OTU_bact_mock_rar_ord=ordinate(GLBRC018_OTU_bact_mock_rar,method = "NMDS")
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
 # stress is (nearly) zero: you may have insufficient data


plot_ordination(GLBRC018_OTU_bact_mock_rar,GLBRC018_OTU_bact_mock_rar_ord,color = "Run",label = "SampleID")
#not great but the best I can do for now. 

#####Remove non bactal reads####
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2)
#32302

sum(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))
#32957920
head(tax_table(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))
unique(data.frame(tax_table(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2))$Domain)
GLBRC018_OTU_bact_MMPRNT_mock_bact=subset_taxa(GLBRC018_OTU_bact_MMPRNT_mock_fin_fil2,Domain!="")
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_bact)
# 1109
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_bact)
#32274
sum(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_bact))
#32949839


GLBRC018_OTU_bact_MMPRNT_mock_bact<-subset_taxa(GLBRC018_OTU_bact_MMPRNT_mock_bact,Class!="c:Chloroplast")
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_bact)
#32064
sum(otu_table(GLBRC018_OTU_bact_MMPRNT_mock_bact))
#32657828

GLBRC018_OTU_bact_MMPRNT_mock_bact<-subset_taxa(GLBRC018_OTU_bact_MMPRNT_mock_bact,Family!="f:Mitochondria")
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_bact)
#1109
ntaxa(GLBRC018_OTU_bact_MMPRNT_mock_bact)
#31950
sum(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_bact))
#32275688
sort(sample_sums(GLBRC018_OTU_bact_MMPRNT_mock_bact))[1:20]
#MMPRNT1920 MMPRNT1913      Amp50 MMPRNT1528 MMPRNT1918 MMPRNT1741 MMPRNT1739 MMPRNT1455 MMPRNT1512     Amp779 MMPRNT1917 MMPRNT1274 MMPRNT1738 
#10215      10796      10958      10986      11176      11235      11263      11660      11878      12470      12582      12695      13015 
#Amp90 MMPRNT1622 MMPRNT1561 MMPRNT1914 MMPRNT1121 MMPRNT1401 MMPRNT1145 
#13328      13552      13586      13605      13649      13729      13780 

#####Rarefy the community####

GLBRC018_OTU_bact_MMPRNT_mock_bact_rar=rarefy_even_depth(GLBRC018_OTU_bact_MMPRNT_mock_bact,replace=F,sample.size=10000,rngseed=99)
#111OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar, file = here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData")
load(here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
#####END: Pre-processing Bacterial community####


#####START: Pre-processing Fungal community####

#Read in the fasta file

rep_set.GLBRC018_fung_OTU<- read.fasta(here::here("Fung_HPCC_out","rep_set_combined_merged_ITS1_2_GLBRC_otus.fa"), as.string = TRUE, set.attributes = FALSE)
head(rep_set.GLBRC018_fung_OTU)

#Load seq mapping
GLBRC018_fung_map=read.csv(here::here("Fung_HPCC_out","Mapping_Fungi_GLBRC_2018.csv"), header = T)
summary(GLBRC018_fung_map)
nrow(GLBRC018_fung_map)
#1563
head(GLBRC018_fung_map)
#Let's combine the map with the metadata
#Metadata not included in this GIT directory. Let me know if you need the raw files
MMPRNT_16s.metadata_mapp_raw=read.table(file = "D:/MMPRNT_16S_016-018/ABIOTIC_data/MMPRNT.metadata_mapp_raw_comp_OUTLIER_filter.txt")

MMPRNT.metadata_mapp_trunc=MMPRNT_16s.metadata_mapp_raw[,c("sampleID_long","collectionDate","siteID","plotType","FertStatus","plotRep","year","dis_Mid","UTM_Lat_Cord","UTM_Lon_Cord")]

MMPRNT.metadata_mapp_trunc[1:10,1:10]
nrow(MMPRNT.metadata_mapp_trunc)
#2700
colnames(GLBRC018_fung_map)[colnames(GLBRC018_fung_map)=="field_sampleID"]="sampleID_long"

GLBRC018_fung_map_metadata=merge(GLBRC018_fung_map,MMPRNT.metadata_mapp_trunc, by="sampleID_long", all.x = T)
head(GLBRC018_fung_map_metadata)
summary(GLBRC018_fung_map_metadata)
colnames(GLBRC018_fung_map_metadata)
row.names(GLBRC018_fung_map_metadata)=GLBRC018_fung_map_metadata$sampleID_seq 
GLBRC018_fung_map_metadata[1:10,1:10]

#OTU table
otu_GLBRC018_fung=otu_table(read.table(here::here("Fung_HPCC_out","combined_merged_ITS1_2_GLBRC_OTU_table.txt"),sep = "\t", header = T, row.names = 1), taxa_are_rows = T)
otu_GLBRC018_fung[1:10,1:10]


#Load Taxa table
#####CONSTAX IS BUGGY####
#I am going to use CONSTAX for taxonomy since it goes deeper and catches some crap sequences
GLBRC018_fung_raw_CONSTAX_UNITE8.2_all = read.delim(here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","combined_taxonomy.txt"),
                                                    sep = "\t",header = T,fill=T, row.names = 1)
head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all)
#8039

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all[GLBRC018_fung_raw_CONSTAX_UNITE8.2_all==""|
                                         is.na(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all)]="Unknown"
head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all)
unique(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Kingdom_SINTAX)

#let's look into the strength of classifications

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phyla_unkn=ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_RDP=="Unknown"|
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_BLAST=="Unknown"|
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_SINTAX=="Unknown", "TRUE","FALSE")


nrow(subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all,Phyla_unkn==TRUE))
#4365

#How many TAXA have a classification at Phyla across the three measures but differ in their phyla classification

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all,Phyla_unkn==FALSE)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW)
#3674


GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phyla_diff=ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phylum_RDP!=
                                                               GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phylum_BLAST|
                                                               GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phylum_RDP!=
                                                               GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phylum_SINTAX|
                                                               GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phylum_BLAST!=
                                                               GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW$Phylum_SINTAX, "TRUE","FALSE")




nrow(subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW,Phyla_diff==TRUE))
#17

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW,Phyla_diff==TRUE)



sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_RDP==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_BLAST,0,1))
#4
sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_SINTAX==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_BLAST,0,1))
#17
sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_RDP==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_SINTAX,0,1))
#13

sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_SINTAX!=
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_BLAST&
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_SINTAX==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif$Phylum_RDP,1,0))
#4

unique(with(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif,interaction(Phylum_SINTAX,Phylum_RDP,Phylum_BLAST)))


data.frame(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif) %>% group_by(Phylum_SINTAX,Phylum_BLAST,Phylum_RDP)%>%summarise_at(vars(Kingdom_BLAST),~n())

#   Phylum_SINTAX   Phylum_BLAST      Phylum_RDP        Kingdom_BLAST
#<chr>           <chr>             <chr>                     <int>
#1 Ascomycota      Basidiomycota     Ascomycota                    2
#2 Ascomycota      Basidiomycota     Basidiomycota                 4
#3 Ascomycota      Glomeromycota     Ascomycota                    1
#4 Ascomycota      Mortierellomycota Mortierellomycota             8
#5 Basidiomycota   Ascomycota        Ascomycota                    1
#6 Chytridiomycota Mortierellomycota Chytridiomycota               1


#Let's output the sintax buggy results 
GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif, Phylum_SINTAX!=Phylum_RDP&Phylum_SINTAX!=Phylum_BLAST)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX)
#13
write.csv(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX, here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_combined_taxonomy.csv"))


#Let's load in the consensus taxonomy to out output the strange OTUs 

GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus = read.delim(here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa/constax_taxonomy.txt"),
                                                          sep = "\t",header = T,fill=T, row.names = 1)
head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus)
#8039


GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus_SINTAX=GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus[row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus)%in%row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX),]
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus_SINTAX)
#13


write.csv(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus_SINTAX, here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_constax_taxonomy.csv"))

#now let's load in the actual Syntax database

GLBRC018_fung_raw_SINTAX_UNITE8.2 = read.delim( here::here("Fung_HPCC_out","taxonomy_combined_merged_ITS1_2_GLBRC_otus.SINTAX"),
                                               sep = "\t",header = F,fill=T, row.names = 1)

head(GLBRC018_fung_raw_SINTAX_UNITE8.2)
nrow(GLBRC018_fung_raw_SINTAX_UNITE8.2)
#8039


GLBRC018_fung_raw_SINTAX_UNITE8.2_SINTAX=GLBRC018_fung_raw_SINTAX_UNITE8.2[row.names(GLBRC018_fung_raw_SINTAX_UNITE8.2)%in%row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX),]
nrow(GLBRC018_fung_raw_SINTAX_UNITE8.2_SINTAX)
#13

write.csv(GLBRC018_fung_raw_SINTAX_UNITE8.2_SINTAX, here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_taxonomy_combined_merged_ITS1_2_GLBRC_otus.csv"))



#total number of phyla classification that differ with unknowns included

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phyla_diff=ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_RDP!=
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_BLAST|
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_RDP!=
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_SINTAX|
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_BLAST!=
                                                           GLBRC018_fung_raw_CONSTAX_UNITE8.2_all$Phylum_SINTAX, "TRUE","FALSE")
nrow(subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all,Phyla_diff==TRUE))
#2276

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all,Phyla_diff==TRUE)



sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_RDP==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_BLAST,0,1))
#1922
sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_BLAST,0,1))
#2001
sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_RDP==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX,0,1))
#648

sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX!=
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_BLAST&
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_SINTAX==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif$Phylum_RDP,1,0))
#1628


#Let's look at the times when RDP is UNKNOWN

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif, 
                                                        Phylum_RDP=="Unknown"&Phylum_BLAST!="Unknown"&
                                                          Phylum_SINTAX!="Unknown")
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B)
#108

sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B$Phylum_SINTAX==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B$Phylum_BLAST,0,1))
#4

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B_SINTAX=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B,
                                                               Phylum_SINTAX!=Phylum_BLAST)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B_SINTAX)
#4


#Let's look at the times when BLAST is UNKNOWN

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif, 
                                                        Phylum_RDP!="Unknown"&Phylum_BLAST=="Unknown"&
                                                          Phylum_SINTAX!="Unknown")
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R)
#329

sum(ifelse(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R$Phylum_SINTAX==
             GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R$Phylum_RDP,0,1))
#10

GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R_SINTAX=subset(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R,
                                                               Phylum_SINTAX!=Phylum_RDP)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R_SINTAX)
#10
GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX=rbind(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_B_SINTAX,GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_S_R_SINTAX)
write.csv(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX, here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_bug_UNKNOWN_combined_taxonomy.csv"))
head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX)

#Let's load in the consensus taxonomy to out output the strange OTUs 


head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus)
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus)
#8039


GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus_SINTAX_UNK=GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus[row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus)%in%row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX),]
nrow(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus_SINTAX_UNK)
#14


write.csv(GLBRC018_fung_raw_CONSTAX_UNITE8.2_concensus_SINTAX_UNK, here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_UNK_bug_constax_taxonomy.csv"))

#now let's load in the actual Syntax database



head(GLBRC018_fung_raw_SINTAX_UNITE8.2)
nrow(GLBRC018_fung_raw_SINTAX_UNITE8.2)
#8039


GLBRC018_fung_raw_SINTAX_UNITE8.2_SINTAX_UNK=GLBRC018_fung_raw_SINTAX_UNITE8.2[row.names(GLBRC018_fung_raw_SINTAX_UNITE8.2)%in%row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX),]
nrow(GLBRC018_fung_raw_SINTAX_UNITE8.2_SINTAX_UNK)
#14

write.csv(GLBRC018_fung_raw_SINTAX_UNITE8.2_SINTAX_UNK, here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","SINTAX_UNK_bug_taxonomy_combined_merged_ITS1_2_GLBRC_otus.csv"))


#I am going to output the rep seqs for these weird OTUs
length(rep_set.GLBRC018_fung_OTU)
#8039

rep_set.GLBRC018_fung_OTU_cons_SINTAX_bug=rep_set.GLBRC018_fung_OTU[names(rep_set.GLBRC018_fung_OTU) %in% 
                                                                      c(row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_P_dif_UNK_SINTAX),
                                                                        row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_all_KNW_P_dif_SINTAX))]
head(rep_set.GLBRC018_fung_OTU_cons_SINTAX_bug)
length(rep_set.GLBRC018_fung_OTU_cons_SINTAX_bug)
#27


write.fasta(sequences =rep_set.GLBRC018_fung_OTU_cons_SINTAX_bug, names = names(rep_set.GLBRC018_fung_OTU_cons_SINTAX_bug), 
            file.out =here::here("Fung_HPCC_out","TEST_constax_V2.1_classification_all_v4.2.2020_taxa","rep_set.GLBRC018_fung_OTU_cons_SINTAX_bug.fna"))

#####CONSTAX IS TROUBLING ME####

#I am going to use the sintax classification for now 


GLBRC018_fung_raw_SINTAX_UNITE8.2 = read.delim( here::here("Fung_HPCC_out","taxonomy_combined_merged_ITS1_2_GLBRC_otus.SINTAX"),
                                                sep = "\t",header = F,fill=T, row.names = 1)

head(GLBRC018_fung_raw_SINTAX_UNITE8.2)
nrow(GLBRC018_fung_raw_SINTAX_UNITE8.2)
#8039
head(GLBRC018_fung_raw_SINTAX_UNITE8.2)
GLBRC018_fung_SINTAX_UNITE8.2_80c=GLBRC018_fung_raw_SINTAX_UNITE8.2[,c(1,4)]
head(GLBRC018_fung_SINTAX_UNITE8.2_80c)

GLBRC018_fung_SINTAX_UNITE8.2_80c_sep=GLBRC018_fung_SINTAX_UNITE8.2_80c %>% separate(V4, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ",")
row.names(GLBRC018_fung_SINTAX_UNITE8.2_80c_sep)=GLBRC018_fung_SINTAX_UNITE8.2_80c_sep$V1
GLBRC018_fung_SINTAX_UNITE8.2_80c_sep$V1=NULL
GLBRC018_fung_SINTAX_UNITE8.2_80c_sep[is.na(GLBRC018_fung_SINTAX_UNITE8.2_80c_sep)] <- "Unknown"
GLBRC018_fung_SINTAX_UNITE8.2_80c_sep_mat=as.matrix(GLBRC018_fung_SINTAX_UNITE8.2_80c_sep)
head(GLBRC018_fung_SINTAX_UNITE8.2_80c_sep_mat)
unique(data.frame(GLBRC018_fung_SINTAX_UNITE8.2_80c_sep_mat)$Domain)
TAXA_80C_UNITE8.2_GLBRC018=tax_table(GLBRC018_fung_SINTAX_UNITE8.2_80c_sep_mat)
head(TAXA_80C_UNITE8.2_GLBRC018)




GLBRC018_OTU_fung_phyl=phyloseq(otu_GLBRC018_fung,sample_data(GLBRC018_fung_map_metadata),TAXA_80C_UNITE8.2_GLBRC018)
ntaxa(GLBRC018_OTU_fung_phyl)
#8039
sum(taxa_sums(GLBRC018_OTU_fung_phyl))
#128045995
mean(sample_sums(GLBRC018_OTU_fung_phyl))
#83254.87
min(sample_sums(GLBRC018_OTU_fung_phyl))
#1
max(sample_sums(GLBRC018_OTU_fung_phyl))
#342074
sort(sample_sums(GLBRC018_OTU_fung_phyl))[1:81]

#MMPRNT2018          Amp1951          Amp1962          Amp1965          Amp1971          Amp1975          Amp1976          Amp1977          Amp1978 
#1                1                1                1                1                1                1                1                1 
#Amp1981          Amp2001          Amp2024          Amp2034          Amp2036          Amp2181       MMPRNT1665          Amp1986       MMPRNT1455 
#1                1                1                1                1                1                2                3                4 
#Amp1966     CONTROLP5B05          Amp1995    CONTROLP09B05          Amp2035          Amp1987          Amp1283          Amp1655        PCRBLANK3 
#4                5                8               27               43              113              151              165              197 
#Amp1355          Amp2119          Amp1307       MMPRNT1857          Amp1348       MMPRNT1146          Amp1368       MMPRNT1450          Amp2072 
#198              199              203              229              254              255              280              326              334 
#Amp1625          Amp1536          Amp2111          Amp1544          Amp1293          Amp1510          Amp1640          Amp1606          Amp2054 
#376              422              432              442              475              495              496              594              596 
#MMPRNT1508    CONTROLP07D07       MMPRNT1860 rPCRMMPRNT1312d1          Amp2051       MMPRNT1149 rPCRMMPRNT1309d1       MMPRNT1311 rPCRMMPRNT1456d2 
#598              604              632              636              659              662              667              707              750 
#Amp2168          Amp2197          Amp1538          Amp2080        PCRBLANK1       MMPRNT1453          Amp1560        PCRBLANK2       MMPRNT1302 
#771              771              789              792              812              828              906              936              948 
#Amp2191 rPCRMMPRNT1304d2          Amp2074          Amp2154          Amp1997 rPCRMMPRNT1309d2          Amp1948       MMPRNT1682       MMPRNT1456 
#1000             1017             1097             1209             1267             1331             1356             1371             1451 
#Amp2152    CONTROLP08G06     CONTROLP6G11       MMPRNT1304    CONTROLP08H03       MMPRNT2035          Amp2110       MMPRNT1963          Amp2109 
#1561             1599             1608             1610             1616             1640             1646             1702             1740 


#####Mock community before removing fungi######

GLBRC018_OTU_fung.mock=subset_samples(GLBRC018_OTU_fung_phyl,seq_samp_type=="mock")
nsamples(GLBRC018_OTU_fung.mock)
#13
GLBRC018_OTU_fung.mock=prune_taxa(taxa_sums(GLBRC018_OTU_fung.mock) > 0, GLBRC018_OTU_fung.mock)
ntaxa(GLBRC018_OTU_fung.mock)
#746


sum(otu_table(GLBRC018_OTU_fung.mock))
#1567679

max(sample_sums(GLBRC018_OTU_fung.mock))
#218299
min(sample_sums(GLBRC018_OTU_fung.mock))
#4

sort(sample_sums(GLBRC018_OTU_fung.mock))
#Amp1966 Amp1281 Amp1591 Amp1360 Amp1495 Amp1936 Amp1584 Amp2102 Amp1487 Amp2070   Mock2 Amp2187   Mock1 
#4   84054   93328   97774  101500  106412  115315  122289  135653  138784  161300  192967  218299 


GLBRC018_OTU_fung.mock_ord=ordinate(GLBRC018_OTU_fung.mock,method = "NMDS")
#Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data


sample_data(GLBRC018_OTU_fung.mock)$spiked_sample=ifelse(sample_data(GLBRC018_OTU_fung.mock)$mock_sampleID=="Mock1","mock","spiked_sample")

plot_ordination(GLBRC018_OTU_fung.mock,GLBRC018_OTU_fung.mock_ord,color = "spiked_sample")


sort(taxa_sums(GLBRC018_OTU_fung.mock))

rep_set.GLBRC018_fung_OTU_mock=rep_set.GLBRC018_fung_OTU[names(rep_set.GLBRC018_fung_OTU) %in% taxa_names(GLBRC018_OTU_fung.mock)]
head(rep_set.GLBRC018_fung_OTU_mock)
length(rep_set.GLBRC018_fung_OTU_mock)
#746


write.fasta(sequences =rep_set.GLBRC018_fung_OTU_mock, names = names(rep_set.GLBRC018_fung_OTU_mock), 
            file.out =here::here("Fung_HPCC_out","rep_set.GLBRC018018_fung_OTU_mock.fna"))

#####Mock USEARCH_global####

#Code run in Usearch

#cd HardDrive/GLBRC-TeamMicrobiome/SG_Microbiome/Fungal_mock/

#The program needs for the codons to be uppercase
#awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' rep_set.GLBRC018018_fung_OTU_mock.fna > CAP_rep_set.GLBRC018018_fung_OTU_mock.fna


#~/HardDrive/Sciencey_Program/Old_usearch_v/usearch10.0.240_win32.exe -usearch_global CAP_rep_set.GLBRC018018_fung_OTU_mock.fna -db Bonito_mock_amptk_synmock.fa -id 0.97 -strand both -maxaccepts 20 -maxrejects 50 -matched match_OTU_GLBRC018_fung_mock.fa -notmatched NOT_match_OTU_GLBRC018_fung_mock.fa -userout match_OTU_GLBRC018_fung_mock.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir


#1.9% matched to mock community sequences


#~/HardDrive/Sciencey_Program/Old_usearch_v/usearch10.0.240_win32.exe -usearch_global CAP_rep_set.GLBRC018018_fung_OTU_mock.fna -db Bonito_mock_amptk_synmock.fa -id 0.8 -strand both -maxaccepts 20 -maxrejects 50 -matched match_OTU_GLBRC018_fung_mock80.fa -notmatched NOT_match_OTU_GLBRC018_fung_mock80.fa -userout match_OTU_GLBRC018_fung_mock80.txt -userfields query+target+id+mid+bits+evalue+ql+ts+qlor+qhir+tlor+thir

#00:00 14Mb    100.0% Searching, 1.9% matched

#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")

matchs_mocks_fung_OTU=read.table(here::here("Fungal_mock","match_OTU_GLBRC018_fung_mock.txt"),header = F)
head(matchs_mocks_fung_OTU)
colnames(matchs_mocks_fung_OTU)=user_out_col
nrow(matchs_mocks_fung_OTU)
#14
length(unique(matchs_mocks_fung_OTU$target))
#12

GLBRC018_OTU_fung.mock_match<-prune_taxa(as.character(unique(matchs_mocks_fung_OTU$query)),
                                         GLBRC018_OTU_fung.mock)
sum(taxa_sums(GLBRC018_OTU_fung.mock_match))/sum(taxa_sums(GLBRC018_OTU_fung.mock))
#0.9315491 match to mock community sequences


#What does the the taxonomy look like 

tax_table(GLBRC018_OTU_fung.mock_match)


#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")

matchs_mocks_fung_OTU_80=read.table(here::here("Fungal_mock","match_OTU_GLBRC018_fung_mock80.txt"),header = F)
head(matchs_mocks_fung_OTU_80)
colnames(matchs_mocks_fung_OTU_80)=user_out_col
nrow(matchs_mocks_fung_OTU_80)
#20
length(unique(matchs_mocks_fung_OTU_80$target))
#12


GLBRC018_OTU_fung.mock_match80<-prune_taxa(as.character(unique(matchs_mocks_fung_OTU_80$query)),
                                           GLBRC018_OTU_fung.mock)
ntaxa(GLBRC018_OTU_fung.mock_match80)
sum(taxa_sums(GLBRC018_OTU_fung.mock_match80))/sum(taxa_sums(GLBRC018_OTU_fung.mock))
#0.9315491

mock_comp_fung_OTU=data.frame("match_sum80"=sample_sums(GLBRC018_OTU_fung.mock_match80),
                              "match_sum97"=sample_sums(GLBRC018_OTU_fung.mock_match),
                              "tot_sum"=sample_sums(GLBRC018_OTU_fung.mock), 
                              estimate_richness(GLBRC018_OTU_fung.mock_match,measures = "Observed"),
                              estimate_richness(GLBRC018_OTU_fung.mock_match80,measures = "Observed"),
                              estimate_richness(GLBRC018_OTU_fung.mock,measures = "Observed"))
colnames(mock_comp_fung_OTU)[4:6]=c("match_rich97","match_rich80","tot_rich")

mock_comp_fung_OTU$prop_read_match97=mock_comp_fung_OTU$match_sum97/mock_comp_fung_OTU$tot_sum
mock_comp_fung_OTU$prop_read_match80=mock_comp_fung_OTU$match_sum80/mock_comp_fung_OTU$tot_sum

mock_comp_fung_OTU$prop_taxa_match97=mock_comp_fung_OTU$match_rich97/mock_comp_fung_OTU$tot_rich
mock_comp_fung_OTU$prop_taxa_match80=mock_comp_fung_OTU$match_rich80/mock_comp_fung_OTU$tot_rich
mock_comp_fung_OTU$sampleID=row.names(mock_comp_fung_OTU)

mock_comp_fung_OTU_m=melt(mock_comp_fung_OTU[,c("sampleID","prop_read_match97","prop_read_match80","prop_taxa_match97","prop_taxa_match80")])

var_order= c("prop_read_match97","prop_taxa_match97","prop_read_match80","prop_taxa_match80")
otu_var_names=c("Reads matching\nat 97%","OTU matching\nat 97%","Reads matching\nat 80%","OTU matching\nat 80%")

ggplot(mock_comp_fung_OTU_m, aes(x=factor(variable, levels = var_order),y=value, group=sampleID, color=sampleID))+geom_point(size=4)+
  scale_x_discrete(name=NULL,label=otu_var_names)+ylab("Precision\n(proportion of sample)")+theme_bw()+theme(axis.text = element_text(size=16),
                                                                                                             axis.title = element_text(size=20))
mock_comp_fung_OTU_m$pipeline=rep("OTU")

#Amp1936 and Amp1966 have no matches to mock sequences

subset(mock_comp_fung_OTU_m, sampleID=="Amp1936"|sampleID=="Amp1966")
subset(mock_comp_fung_OTU, sampleID=="Amp1936"|sampleID=="Amp1966")


ggplot(subset(mock_comp_fung_OTU_m, sampleID!="Amp1936"&sampleID!="Amp1966"), 
       aes(x=factor(variable, levels = var_order),y=value, group=sampleID, color=sampleID))+geom_point(size=4)+
  scale_x_discrete(name=NULL,label=otu_var_names)+ylab("Precision\n(proportion of sample)")+geom_text(aes(label=sampleID),position = position_dodge(0.5))+
  theme_bw()+theme(axis.text = element_text(size=16),axis.title = element_text(size=20))


mean(subset(mock_comp_fung_OTU_m,sampleID!="Amp1936"&sampleID!="Amp1966"&variable=="prop_read_match97")$value)
summary(subset(mock_comp_fung_OTU_m,sampleID!="Amp1936"&sampleID!="Amp1966"&variable=="prop_read_match97")$value)

#Let remove the two samples

mock_comp_fung_OTU_m_filt=subset(mock_comp_fung_OTU_m,sampleID!="Amp1936"&sampleID!="Amp1966")

#Let's look at spike in versus pure mocks

GLBRC018_OTU_fung.mock_meta=sample_data(GLBRC018_OTU_fung.mock)
GLBRC018_OTU_fung.mock_meta$spiked_sample=ifelse(GLBRC018_OTU_fung.mock_meta$mock_sampleID=="Mock1","mock","spiked_sample")

mock_comp_fung_OTU_m_filt_s=merge(mock_comp_fung_OTU_m_filt,data.frame(GLBRC018_OTU_fung.mock_meta), by.x="sampleID",by.y="sampleID_seq")


ggplot(mock_comp_fung_OTU_m_filt_s, aes(x=spiked_sample,y=value, group=spiked_sample, color=spiked_sample))+geom_point(size=4)+
  facet_wrap(~factor(variable,label=otu_var_names, levels = var_order),scales = "free_y")+scale_x_discrete(name=NULL,label=c("Mock Only","Spiked with sample"))+
  ylab("Precision\n(proportion of sample)")+theme_bw()+theme(axis.text = element_text(size=16),axis.title = element_text(size=20))


#Let's ordinate this 
nsamples(GLBRC018_OTU_fung.mock)
#13
GLBRC018_OTU_fung.mock_fil=subset_samples(GLBRC018_OTU_fung.mock,sampleID_seq!="Amp1936"&sampleID_seq!="Amp1966")
nsamples(GLBRC018_OTU_fung.mock_fil)
#11
GLBRC018_OTU_fung.mock_fil_ord=ordinate(GLBRC018_OTU_fung.mock_fil,method = "NMDS")
#Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data


plot_ordination(GLBRC018_OTU_fung.mock_fil,GLBRC018_OTU_fung.mock_fil_ord,color = "spiked_sample",label = "sampleID_seq")


sort(taxa_sums(GLBRC018_OTU_fung.mock))



#####Taxanomic Likages####

#OTU Table from 97%

GLBRC018_OTU_fung.mock_match_TBL=data.frame(otu_table(GLBRC018_OTU_fung.mock_match))
GLBRC018_OTU_fung.mock_match_TBL$query=row.names(GLBRC018_OTU_fung.mock_match_TBL)
head(GLBRC018_OTU_fung.mock_match_TBL)

GLBRC018_OTU_fung.mock_match_CLASS=merge(GLBRC018_OTU_fung.mock_match_TBL,matchs_mocks_fung_OTU[,c("query","target","id")],
                                         by="query")

GLBRC018_OTU_fung.mock_match_CLASS$ASV_sum=rowSums(GLBRC018_OTU_fung.mock_match_CLASS[,sample_names(GLBRC018_OTU_fung.mock_match)])


GLBRC018_OTU_fung.mock_match_CLASS_sum=GLBRC018_OTU_fung.mock_match_CLASS %>% group_by(target)%>%summarise_at(vars(sample_names(GLBRC018_OTU_fung.mock_match)),~nnzero(.))


#OTU Table from 80%

GLBRC018_OTU_fung.mock_matchh80_TBL=data.frame(otu_table(GLBRC018_OTU_fung.mock_match80))
GLBRC018_OTU_fung.mock_matchh80_TBL$query=row.names(GLBRC018_OTU_fung.mock_matchh80_TBL)
head(GLBRC018_OTU_fung.mock_matchh80_TBL)

GLBRC018_OTU_fung.mock_matchh80_CLASS=merge(GLBRC018_OTU_fung.mock_matchh80_TBL,matchs_mocks_fung_OTU_80[,c("query","target","id")],
                                            by="query")

GLBRC018_OTU_fung.mock_matchh80_CLASS$ASV_sum=rowSums(GLBRC018_OTU_fung.mock_matchh80_CLASS[,sample_names(GLBRC018_OTU_fung.mock_match)])


GLBRC018_OTU_fung.mock_matchh80_CLASS_sum=GLBRC018_OTU_fung.mock_matchh80_CLASS %>% group_by(target)%>%summarise_at(vars(sample_names(GLBRC018_OTU_fung.mock_match)),~nnzero(.))


GLBRC018_OTU_fung.mock_match_over_spliting=data.frame(t(rbind(colSums(GLBRC018_OTU_fung.mock_match_CLASS_sum[,sample_names(GLBRC018_OTU_fung.mock_match)] != 0),
                                                              colSums(GLBRC018_OTU_fung.mock_match_CLASS_sum[,sample_names(GLBRC018_OTU_fung.mock_match)]),
                                                              colSums(GLBRC018_OTU_fung.mock_matchh80_CLASS_sum[,sample_names(GLBRC018_OTU_fung.mock_match)] != 0),
                                                              colSums(GLBRC018_OTU_fung.mock_matchh80_CLASS_sum[,sample_names(GLBRC018_OTU_fung.mock_match)]))))


colnames(GLBRC018_OTU_fung.mock_match_over_spliting)=c("mock_taxa97","ASV_taxa97","mock_taxa80","ASV_taxa80")


GLBRC018_OTU_fung.mock_match_over_spliting$over_split97=GLBRC018_OTU_fung.mock_match_over_spliting$ASV_taxa97/
  GLBRC018_OTU_fung.mock_match_over_spliting$mock_taxa97
GLBRC018_OTU_fung.mock_match_over_spliting$over_split80=GLBRC018_OTU_fung.mock_match_over_spliting$ASV_taxa80/
  GLBRC018_OTU_fung.mock_match_over_spliting$mock_taxa80

GLBRC018_OTU_fung.mock_match_over_spliting$sampleID=row.names(GLBRC018_OTU_fung.mock_match_over_spliting)

GLBRC018_OTU_fung.mock_match_over_spliting_m<-melt(GLBRC018_OTU_fung.mock_match_over_spliting)
GLBRC018_OTU_fung.mock_match_over_spliting_m$pipeline=rep("OTU")
#Using the same metrics that Pauvert et al 2019 used

#pull out the ASV with the max frequnence for each mock taxa


GLBRC018_OTU_fung.mock_match_CLASS_sub=GLBRC018_OTU_fung.mock_match_CLASS %>%
  group_by(target) %>%
  filter(row_number() == which.max(ASV_sum))

GLBRC018_OTU_fung.mock_matchh80_CLASS_sub=GLBRC018_OTU_fung.mock_matchh80_CLASS %>%
  group_by(target) %>%
  filter(row_number() == which.max(ASV_sum))


#Sensitivity 
#true positive rate TP/(TP + FN)
#TP is the number of true-positive OTUs (or ASVs)
#True-positive OTUs corresponded to fungal
#strains present in the mock community and identified by the bioinformatic
#approach considered
#FN is the number of falsenegative OTUs (or ASVs)
#False-negative OTUs corresponded
#to fungal strains present in the mock community but not
#detected by the bioinformatic approach considered. 

#Precision
#positive predictive value TP/(TP ? FP)
#TP is the number of true-positive OTUs (or ASVs)
#True-positive OTUs corresponded to fungal
#strains present in the mock community and identified by the bioinformatic
#approach considered
#FP is the number of false-positive OTUs (or ASVs).
#False-positive
#OTUs corresponded to all other OTUs. If several OTUs were
#assigned to the same fungal strain of the mock community (i.e.
#'split' OTUs), only the most abundant was considered to be a truepositive
#OTU, the others being considered false-positive OTUs.

GLBRC018_OTU_fung_sens_per=rbind(data.frame(GLBRC018_OTU_fung.mock_match_CLASS_sub %>% group_by(target)%>%
                                              summarise_at(vars(sample_names(GLBRC018_OTU_fung.mock_match)),~nnzero(.))%>%
                                              summarise_at(vars(sample_names(GLBRC018_OTU_fung.mock_match)),~sum(.))),
                                 colSums(GLBRC018_OTU_fung.mock_match_CLASS_sub[,c(sample_names(GLBRC018_OTU_fung.mock_match))]),
                                 data.frame(GLBRC018_OTU_fung.mock_matchh80_CLASS_sub %>% group_by(target)%>%
                                              summarise_at(vars(sample_names(GLBRC018_OTU_fung.mock_match)),~nnzero(.))%>%
                                              summarise_at(vars(sample_names(GLBRC018_OTU_fung.mock_match)),~sum(.))),
                                 colSums(GLBRC018_OTU_fung.mock_matchh80_CLASS_sub[,c(sample_names(GLBRC018_OTU_fung.mock_match))]),
                                 sample_sums(GLBRC018_OTU_fung.mock),
                                 t(estimate_richness(GLBRC018_OTU_fung.mock,measures = "Observed")))


row.names(GLBRC018_OTU_fung_sens_per)=c("tp_rich97","tp_reads97","tp_rich80","tp_reads80","tot_sum","tot_rich")



GLBRC018_OTU_fung_sens_per_t=data.frame(t(GLBRC018_OTU_fung_sens_per))


GLBRC018_OTU_fung_sens_per_t$sensitivity97=GLBRC018_OTU_fung_sens_per_t$tp_rich97/12

GLBRC018_OTU_fung_sens_per_t$precision97=GLBRC018_OTU_fung_sens_per_t$tp_rich97/GLBRC018_OTU_fung_sens_per_t$tot_rich

GLBRC018_OTU_fung_sens_per_t$TP_prop_read_match97=GLBRC018_OTU_fung_sens_per_t$tp_reads97/GLBRC018_OTU_fung_sens_per_t$tot_sum

GLBRC018_OTU_fung_sens_per_t$sensitivity80=GLBRC018_OTU_fung_sens_per_t$tp_rich80/12

GLBRC018_OTU_fung_sens_per_t$precision80=GLBRC018_OTU_fung_sens_per_t$tp_rich80/GLBRC018_OTU_fung_sens_per_t$tot_rich

GLBRC018_OTU_fung_sens_per_t$TP_prop_read_match80=GLBRC018_OTU_fung_sens_per_t$tp_reads80/GLBRC018_OTU_fung_sens_per_t$tot_sum



GLBRC018_OTU_fung_sens_per_t$sampleID=row.names(GLBRC018_OTU_fung_sens_per_t)
GLBRC018_OTU_fung_sens_per_M<-melt(GLBRC018_OTU_fung_sens_per_t[,c("sampleID","sensitivity97","precision97","TP_prop_read_match97",
                                                                   "sensitivity80","precision80","TP_prop_read_match80","tot_rich")])

GLBRC018_OTU_fung_sens_per_M$pipeline=rep("OTU")


ggplot(GLBRC018_OTU_fung_sens_per_t, aes(x=precision97,y=tot_rich, color=sampleID))+geom_point(size=4)+
  ylab("Similarity to Mock")+theme_bw()+theme(axis.text = element_text(size=16),
                                              axis.title = element_text(size=20))

#####Similarity based on Bray dist####

OTU.mock_samp_sum=sample_sums(GLBRC018_OTU_fung.mock)
length(OTU.mock_samp_sum)

sim_mock_com_OTU=data.frame(rep(OTU.mock_samp_sum[1]/12,12),
                            rep(OTU.mock_samp_sum[2]/12,12),
                            rep(OTU.mock_samp_sum[3]/12,12),
                            rep(OTU.mock_samp_sum[4]/12,12),
                            rep(OTU.mock_samp_sum[5]/12,12),
                            rep(OTU.mock_samp_sum[6]/12,12),
                            rep(OTU.mock_samp_sum[7]/12,12),
                            rep(OTU.mock_samp_sum[8]/12,12),
                            rep(OTU.mock_samp_sum[9]/12,12),
                            rep(OTU.mock_samp_sum[10]/12,12),
                            rep(OTU.mock_samp_sum[11]/12,12),
                            rep(OTU.mock_samp_sum[12]/12,12),
                            rep(OTU.mock_samp_sum[13]/12,12),
                            "query"=c(GLBRC018_OTU_fung.mock_match_CLASS_sub$query))
colnames(sim_mock_com_OTU)=c(paste("sim_",names(OTU.mock_samp_sum), sep = ""),"query")

GLBRC018_OTU_fung.mock_TBL=data.frame(otu_table(GLBRC018_OTU_fung.mock))
GLBRC018_OTU_fung.mock_TBL$query=row.names(GLBRC018_OTU_fung.mock_TBL)

sim_mock_com_v_actual_OTU=merge(sim_mock_com_OTU,GLBRC018_OTU_fung.mock_TBL, by= "query",all=T)
sim_mock_com_v_actual_OTU[is.na(sim_mock_com_v_actual_OTU)] <- 0

sim_mock_com_v_actual_OTU_dist=vegdist(t(sim_mock_com_v_actual_OTU[,c(paste("sim_",names(OTU.mock_samp_sum), sep = ""),
                                                                      names(OTU.mock_samp_sum))]),
                                       method="bray")

sim_mock_com_v_actual_OTU_dist_M <- matrixConvert(sim_mock_com_v_actual_OTU_dist, 
                                                  colname = c("sample1", "sample2", "bray"))


sim_mock_com_v_actual_OTU_dist_M$sample1_2=with(sim_mock_com_v_actual_OTU_dist_M,
                                                interaction(sample1,sample2))

sim_mock_com_v_actual_OTU_dist_M_sub=subset(sim_mock_com_v_actual_OTU_dist_M,
                                            sample1_2=="sim_Mock1.Mock1"|sample1_2=="sim_Mock2.Mock2"|
                                              sample1_2=="sim_Amp1281.Amp1281"|sample1_2=="sim_Amp1360.Amp1360"|
                                              sample1_2=="sim_Amp1487.Amp1487"|sample1_2=="sim_Amp1495.Amp1495"|
                                              sample1_2=="sim_Amp1584.Amp1584"|sample1_2=="sim_Amp1591.Amp1591"|
                                              sample1_2=="sim_Amp1936.Amp1936"|sample1_2=="sim_Amp1966.Amp1966"|
                                              sample1_2=="sim_Amp2070.Amp2070"|sample1_2=="sim_Amp2102.Amp2102"|sample1_2=="sim_Amp2187.Amp2187")
nrow(sim_mock_com_v_actual_OTU_dist_M_sub)
#13
sim_mock_com_v_actual_OTU_dist_M_sub$similarity97=1-sim_mock_com_v_actual_OTU_dist_M_sub$bray



ggplot(sim_mock_com_v_actual_OTU_dist_M_sub, aes(x=sample2,y=similarity97, color=sample2))+geom_point(size=4)+
  ylab("Similarity to Mock")+theme_bw()+theme(axis.text = element_text(size=16),
                                              axis.title = element_text(size=20))

#####How frequent are mock taxa in the other samples?#####
#97% match to mock community sequences

#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")
matchs_mocks_fung_OTU=read.table(here::here("Fungal_mock","match_OTU_GLBRC018_fung_mock.txt"),header = F)
head(matchs_mocks_fung_OTU)
colnames(matchs_mocks_fung_OTU)=user_out_col
nrow(matchs_mocks_fung_OTU)
#14
length(unique(matchs_mocks_fung_OTU$target))
#12

GLBRC018_OTU_fung_phyl_match<-prune_taxa(as.character(unique(matchs_mocks_fung_OTU$query)),
                                         GLBRC018_OTU_fung_phyl)



#Remove the mock samples and samples with no reads
GLBRC018_OTU_fung.NO_mock_match=subset_samples(GLBRC018_OTU_fung_phyl_match,seq_samp_type!="mock")
nsamples(GLBRC018_OTU_fung.NO_mock_match)
#1525
GLBRC018_OTU_fung.NO_mock_match=prune_samples(sample_sums(GLBRC018_OTU_fung.NO_mock_match) > 0, GLBRC018_OTU_fung.NO_mock_match)
nsamples(GLBRC018_OTU_fung.NO_mock_match)
#533
ntaxa(GLBRC018_OTU_fung.NO_mock_match)
#14

hist(sample_sums(GLBRC018_OTU_fung.NO_mock_match))

sort(sample_sums(GLBRC018_OTU_fung.NO_mock_match),decreasing = T)[1:50]
#Amp2009       MMPRNT1921       MMPRNT1115       MMPRNT1818          Amp1578       MMPRNT1922          Amp2123        PCRBLANK2       MMPRNT2027 rPCRMMPRNT1456d2     CONTROLP6G11       MMPRNT1738 
#121218               22               12               11               11               10                9                7                6                6                5                5 
#MMPRNT1811 rPCRMMPRNT1453d2 rPCRMMPRNT1456d1 rPCRMMPRNT1545d2          Amp1596          Amp1885          Amp1888          Amp2030          Amp2112     CONTROLP5A10       MMPRNT1824       MMPRNT2004 
#5                5                5                5                5                5                5                5                5                4                4                4 
#MMPRNT2060 rPCRMMPRNT1364d1 rPCRMMPRNT1481d2          Amp2147          Amp2152          Amp2213     CONTROLP6F09       MMPRNT1709       MMPRNT1804       MMPRNT1810       MMPRNT1817       MMPRNT1825 
#4                4                4                4                4                4                3                3                3                3                3                3 
#MMPRNT1306       MMPRNT1456    CONTROLP08B09    CONTROLP09A07       MMPRNT1538       MMPRNT1539       MMPRNT2023       MMPRNT2025       MMPRNT2031       MMPRNT2037       MMPRNT2063       MMPRNT2066 
#3                3                3                3                3                3                3                3                3                3                3                3 
#MMPRNT2072       MMPRNT2073 
#3                3 

GLBRC018_OTU_fung.NO_mock_match_SUM_DF=data.frame(sample_sums(GLBRC018_OTU_fung.NO_mock_match))
colnames(GLBRC018_OTU_fung.NO_mock_match_SUM_DF)="mock_reads"


GLBRC018_OTU_fung._SUM_DF=data.frame(sample_sums(GLBRC018_OTU_fung_phyl))
colnames(GLBRC018_OTU_fung._SUM_DF)="tot_reads"

GLBRC018_OTU_fung.prop_mock_SUM_DF=merge(GLBRC018_OTU_fung.NO_mock_match_SUM_DF,GLBRC018_OTU_fung._SUM_DF,by="row.names")
head(GLBRC018_OTU_fung.prop_mock_SUM_DF)

GLBRC018_OTU_fung.prop_mock_SUM_DF$prop_mock=GLBRC018_OTU_fung.prop_mock_SUM_DF$mock_reads/GLBRC018_OTU_fung.prop_mock_SUM_DF$tot_reads

GLBRC018_OTU_fung.prop_mock_SUM_DF[order(GLBRC018_OTU_fung.prop_mock_SUM_DF$prop_mock,decreasing = T),][1:40,]




ggplot(GLBRC018_OTU_fung.prop_mock_SUM_DF,aes(x=prop_mock))+geom_histogram()+
  scale_x_continuous(name = "Proportion of reads matching mock per sample")+ylab("Frequency")+
  theme_classic()+theme(axis.title = element_text(size = 20),axis.text = element_text(size = 16))


GLBRC018_OTU_fung.prop_mock_SUM_DF_meta=merge(GLBRC018_OTU_fung.prop_mock_SUM_DF,sample_data(GLBRC018_OTU_fung.NO_mock_match),by.x="Row.names",by.y="row.names")
GLBRC018_OTU_fung.prop_mock_SUM_DF_meta[order(GLBRC018_OTU_fung.prop_mock_SUM_DF_meta$prop_mock,decreasing = T),][1:40,]




#Here are the three weird samples
#
sample_sums(subset_samples(GLBRC018_OTU_fung_phyl, sampleID_seq=="Amp1936"|sampleID_seq=="Amp1966"|sampleID_seq=="Amp2009"))
sample_sums(subset_samples(GLBRC018_OTU_fung_phyl_match, sampleID_seq=="Amp1936"|sampleID_seq=="Amp1966"|sampleID_seq=="Amp2009"))

#Let's remove the Amp2009
GLBRC018_OTU_fung.prop_mock_SUM_DF_meta_fil=subset(GLBRC018_OTU_fung.prop_mock_SUM_DF_meta,Row.names!="Amp2009")


ggplot(GLBRC018_OTU_fung.prop_mock_SUM_DF_meta_fil,aes(x=prop_mock))+geom_histogram()+
  scale_x_continuous(name = "Proportion of reads matching mock per sample")+ylab("Frequency")+
  theme_classic()+theme(axis.title = element_text(size = 20),axis.text = element_text(size = 16))


#####Mock based taxa control#####

#Let's try to figure out how to eliminate the Contamination and Crap OTUs

ntaxa(GLBRC018_OTU_fung.mock)
#746
sum(otu_table(GLBRC018_OTU_fung.mock))
#1567679

max(sample_sums(GLBRC018_OTU_fung.mock))
#218299
min(sample_sums(GLBRC018_OTU_fung.mock))
#4

sort(sample_sums(GLBRC018_OTU_fung.mock))
#Amp1966 Amp1281 Amp1591 Amp1360 Amp1495 Amp1936 Amp1584 Amp2102 Amp1487 Amp2070   Mock2 Amp2187   Mock1 
#4   84054   93328   97774  101500  106412  115315  122289  135653  138784  161300  192967  218299 


GLBRC018_OTU_fung.mock_fil=subset_samples(GLBRC018_OTU_fung.mock,sampleID_seq!="Amp1936"&sampleID_seq!="Amp1966")
nsamples(GLBRC018_OTU_fung.mock_fil)
#11
GLBRC018_OTU_fung.mock_fil=prune_taxa(taxa_sums(GLBRC018_OTU_fung.mock_fil) > 0, GLBRC018_OTU_fung.mock_fil)

sum(otu_table(GLBRC018_OTU_fung.mock_fil))
#1461263

max(sample_sums(GLBRC018_OTU_fung.mock_fil))
#218299
min(sample_sums(GLBRC018_OTU_fung.mock_fil))
#84054

sort(sample_sums(GLBRC018_OTU_fung.mock_fil))
#Amp1281 Amp1591 Amp1360 Amp1495 Amp1584 Amp2102 Amp1487 Amp2070   Mock2 Amp2187   Mock1 
#84054   93328   97774  101500  115315  122289  135653  138784  161300  192967  218299 

otu_table(GLBRC018_OTU_fung.mock_fil)


#Let's remove taxa that only occur in 2 sample 

GLBRC018_OTU_fung.mock_fil_pa=transform_sample_counts(GLBRC018_OTU_fung.mock_fil,pa)
max(otu_table(GLBRC018_OTU_fung.mock_fil_pa))

GLBRC018_OTU_fung.mock_fil_v2=prune_taxa(taxa_sums(GLBRC018_OTU_fung.mock_fil_pa) > 2, GLBRC018_OTU_fung.mock_fil)
ntaxa(GLBRC018_OTU_fung.mock_fil_v2)
#44
sum(otu_table(GLBRC018_OTU_fung.mock_fil_v2))
#1460674

max(sample_sums(GLBRC018_OTU_fung.mock_fil_v2))
#218071
min(sample_sums(GLBRC018_OTU_fung.mock_fil_v2))
#84046
otu_table(GLBRC018_OTU_fung.mock_fil_v2)
tax_table(GLBRC018_OTU_fung.mock_fil_v2)


#Let's ordinate them
GLBRC018_OTU_fung.mock_fil_v2_ord=ordinate(GLBRC018_OTU_fung.mock_fil_v2,method = "NMDS")
#*** Solution reached

#Effects of spike in
sample_data(GLBRC018_OTU_fung.mock_fil_v2)$spiked_sample=ifelse(sample_data(GLBRC018_OTU_fung.mock_fil_v2)$mock_sampleID=="Mock1","mock","spiked_sample")

plot_ordination(GLBRC018_OTU_fung.mock_fil_v2,GLBRC018_OTU_fung.mock_fil_v2_ord)+geom_point(aes(color = spiked_sample), size=3)+geom_text(aes(label = sampleID_seq))

adonis(distance(GLBRC018_OTU_fung.mock_fil_v2,method = "bray")~sample_data(GLBRC018_OTU_fung.mock_fil_v2)$spiked_sample,permutations = 9999)
#not significant
#sample_data(GLBRC018_OTU_fung.mock_fil_v2)$spiked_sample  1  0.003145 0.0031449 0.13121 0.01437 0.7899

#Effects of DNA plate 

plot_ordination(GLBRC018_OTU_fung.mock_fil_v2,GLBRC018_OTU_fung.mock_fil_v2_ord)+geom_point(aes(color = DNA_plate), size=3)+geom_text(aes(label = sampleID_seq))

adonis(distance(GLBRC018_OTU_fung.mock_fil_v2,method = "bray")~sample_data(GLBRC018_OTU_fung.mock_fil_v2)$DNA_plate,permutations = 9999)
#sample_data(GLBRC018_OTU_fung.mock_fil_v2)$DNA_plate  5   0.18684 0.037367  5.8349 0.85369 0.0331 *
#not super strong

#Effects of library

plot_ordination(GLBRC018_OTU_fung.mock_fil_v2,GLBRC018_OTU_fung.mock_fil_v2_ord)+geom_point(aes(color = Library_name), size=3)+geom_text(aes(label = sampleID_seq))

adonis(distance(GLBRC018_OTU_fung.mock_fil_v2,method = "bray")~sample_data(GLBRC018_OTU_fung.mock_fil_v2)$Library_name,permutations = 9999)
#sample_data(GLBRC018_OTU_fung.mock_fil_v2)$Library_name  2  0.144851 0.072425  7.8293 0.66186 0.0134 *
#not super strong


#Rarefying the community

GLBRC018_OTU_fung.mock_fil_v2_rar=rarefy_even_depth(GLBRC018_OTU_fung.mock_fil_v2,replace=F,sample.size=84000,rngseed=99)
ntaxa(GLBRC018_OTU_fung.mock_fil_v2_rar)
#43
sum(otu_table(GLBRC018_OTU_fung.mock_fil_v2_rar))
#924000

max(sample_sums(GLBRC018_OTU_fung.mock_fil_v2_rar))
#84000
min(sample_sums(GLBRC018_OTU_fung.mock_fil_v2_rar))
#84000
otu_table(GLBRC018_OTU_fung.mock_fil_v2_rar)
tax_table(GLBRC018_OTU_fung.mock_fil_v2_rar)


#Less ordinate them
GLBRC018_OTU_fung.mock_fil_v2_rar_ord=ordinate(GLBRC018_OTU_fung.mock_fil_v2_rar,method = "NMDS")
#*** Solution reached

#Effects of spike in
sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$spiked_sample=ifelse(sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$mock_sampleID=="Mock1","mock","spiked_sample")

plot_ordination(GLBRC018_OTU_fung.mock_fil_v2_rar,GLBRC018_OTU_fung.mock_fil_v2_rar_ord)+geom_point(aes(color = spiked_sample), size=3)+geom_text(aes(label = sampleID_seq))

adonis(distance(GLBRC018_OTU_fung.mock_fil_v2_rar,method = "bray")~sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$spiked_sample,permutations = 9999)
#not significant
#sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$spiked_sample  1 0.0006139 0.00061392  1.1306 0.1116 0.3368

#Effects of DNA plate 

plot_ordination(GLBRC018_OTU_fung.mock_fil_v2_rar,GLBRC018_OTU_fung.mock_fil_v2_rar_ord)+geom_point(aes(color = DNA_plate), size=3)+geom_text(aes(label = sampleID_seq))

adonis(distance(GLBRC018_OTU_fung.mock_fil_v2_rar,method = "bray")~sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$DNA_plate,permutations = 9999)
#sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$DNA_plate  5 0.0037897 0.00075794  2.2146 0.68892 0.0184 *
#not super strong

#Effects of library

plot_ordination(GLBRC018_OTU_fung.mock_fil_v2_rar,GLBRC018_OTU_fung.mock_fil_v2_rar_ord)+geom_point(aes(color = Library_name), size=3)+geom_text(aes(label = sampleID_seq))

adonis(distance(GLBRC018_OTU_fung.mock_fil_v2_rar,method = "bray")~sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$Library_name,permutations = 9999)
#sample_data(GLBRC018_OTU_fung.mock_fil_v2_rar)$Library_name  2 0.0024329 0.0012165   3.172 0.44228 0.0048 **
#not super strong

#####MMPRNT Only ######
head(sample_data(GLBRC018_OTU_fung_phyl))
unique(sample_data(GLBRC018_OTU_fung_phyl)$project)

GLBRC018_OTU_fung_MMPRNT_mock=subset_samples(GLBRC018_OTU_fung_phyl,project=="Field2018"|
                                               seq_samp_type=="mock")
GLBRC018_OTU_fung_MMPRNT_mock=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock,tar_region=="ITS1"&Sample_or_Control!="control")
GLBRC018_OTU_fung_MMPRNT_mock=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock) > 0, GLBRC018_OTU_fung_MMPRNT_mock)
nsamples(GLBRC018_OTU_fung_MMPRNT_mock)
#1161

ntaxa(GLBRC018_OTU_fung_MMPRNT_mock)
# 6997
sum(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock))
#108439086
mean(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock))
#93401.45
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock))
#1
max(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock))
#342074


#Let's choose the most the best 
unique(sample_data(GLBRC018_OTU_fung_MMPRNT_mock)$rep_pcr) 


GLBRC018_OTU_fung_MMPRNT_mock_re_PCR=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock,rep_pcr!="none")
sample_data(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)
#76

unique(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)$sampleID_long) 


GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_diversity=merge(estimate_richness(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR),
                                                     data.frame(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)),by="row.names")
head(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_diversity)

colnames(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_diversity)[colnames(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_diversity)==
                                                           "sample_sums.GLBRC018_OTU_fung_MMPRNT_mock_re_PCR."]="sample_sum"
GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta=merge(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_diversity,
                                                    data.frame(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)),
                                                    by.x = "Row.names",by.y = "row.names")
head(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta)

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta,aes(x=sampleID_long,y=sample_sum,fill=rep_pcr))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")


ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta,aes(x=sampleID_long,y=Observed,fill=rep_pcr))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")


#let's see if the most reads means it is the PCR

GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk=GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta%>%group_by(sampleID_long) %>%
  mutate(PCR_rank = order(order(sample_sum, decreasing=TRUE)))

unique(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk$PCR_rank)

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=sample_sum,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")


#how does the diversity look

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=Observed,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=Shannon ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk,aes(x=sampleID_long,y=InvSimpson ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

#They look okay overall. 
#I want to look at if Mock reads somehow made it into these samples 


#Hits
#https://drive5.com/usearch/manual/userfields.html
user_out_col=c("query","target","id","mid","bits","evalue","ql","ts","qlor","qhir","tlor","thir")
matchs_mocks_fung_OTU=read.table(here::here("Fungal_mock","match_OTU_GLBRC018_fung_mock.txt"),header = F)
head(matchs_mocks_fung_OTU)
colnames(matchs_mocks_fung_OTU)=user_out_col
nrow(matchs_mocks_fung_OTU)
#14
length(unique(matchs_mocks_fung_OTU$target))
#12

GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_match<-prune_taxa(as.character(unique(matchs_mocks_fung_OTU$query)),
                                                       GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)


hist(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_match))

sort(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_match),decreasing = T)


GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk_2=merge(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk,
                                                          sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_match),
                                                          by.x = "Row.names",by.y = "row.names")
head(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk_2)

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk_2,aes(x=sampleID_long,y=y ,fill=factor(PCR_rank)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

ggplot(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk_2,aes(x=sampleID_long,y=y ,fill=factor(rep_pcr)))+
  geom_bar(stat = "identity",position = position_dodge(0.5))+facet_wrap(~Root_soil, scales = "free")

#seems okay


#Let's make a list of samples that we want to keep

prune_PCR_fung_samp=subset(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_div_meta_rnk_2,PCR_rank!=1)$Row.names
length(prune_PCR_fung_samp)
#49
length(sample_names(GLBRC018_OTU_fung_MMPRNT_mock))
#1161

final_soil_fung_samples=sample_names(GLBRC018_OTU_fung_MMPRNT_mock)[sample_names(GLBRC018_OTU_fung_MMPRNT_mock) %w/o%prune_PCR_fung_samp]
length(final_soil_fung_samples)
#1112


GLBRC018_OTU_fung_MMPRNT_mock_fin=prune_samples(final_soil_fung_samples,GLBRC018_OTU_fung_MMPRNT_mock)
GLBRC018_OTU_fung_MMPRNT_mock_fin=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin) > 0, GLBRC018_OTU_fung_MMPRNT_mock_fin)
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fin)
#1112

ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fin)
# 6884
sum(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin))
#106773889
mean(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin))
#96019.68
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin))
#1
max(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin))
#342074

#Let's check to make sure we have removed the rePCR samples


nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fin,Root_soil=="Root"&seq_samp_type!="mock"))
#234
length(unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fin,Root_soil=="Root"&seq_samp_type!="mock"))$sampleID_long))
#234

nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fin,Root_soil=="Soil"&seq_samp_type!="mock"))
#865
length(unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fin,Root_soil=="Soil"&seq_samp_type!="mock"))$sampleID_long))
#865




#####Singletons removal ######

#Let remove the taxon that occur in less than 5 samples 
GLBRC018_OTU_fung_MMPRNT_mock_fin_pa=transform_sample_counts(GLBRC018_OTU_fung_MMPRNT_mock_fin,pa)


GLBRC018_OTU_fung_MMPRNT_mock_fin_fil=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_pa) > 4, GLBRC018_OTU_fung_MMPRNT_mock_fin)
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil)
# 6562
sum(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil))
#106718810
mean(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil))
#95970.15
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil))
#1
max(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil))
#342073
sort(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil))[1:40]

#MMPRNT2018       MMPRNT1665       MMPRNT1455          Amp1966       MMPRNT1857       MMPRNT1146       MMPRNT1450       MMPRNT1860       MMPRNT1149 
#1                2                4                4              229              255              326              632              662 
#Amp1560       MMPRNT1682       MMPRNT2035       MMPRNT1963       MMPRNT1309 rPCRMMPRNT1456d1       MMPRNT1306          Amp1345 rPCRMMPRNT1424d1 
#906             1371             1640             1702             2439             2725             4413             5287             5349 
#rPCRMMPRNT1304d1          Amp1654       MMPRNT1312          Amp1333 rPCRMMPRNT1508d1       MMPRNT1114          Amp1343          Amp1346          Amp1540 
#6042             6341             7586             8017             9128             9890            16186            20172            22175 

#Let's remove the samples that are less than 10000 reads 

GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2=prune_samples(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil)>10000,GLBRC018_OTU_fung_MMPRNT_mock_fin_fil)

GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2) > 0, GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2)
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2)
# 6562
GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2_pa=transform_sample_counts(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2,pa)
GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2_pa) > 4, GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2)
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2)
# 6559
sum(taxa_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))
#106642973
mean(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))
#98017.44
min(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))
#16186
max(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))
#342073
sort(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))[1:10]
#Amp1343    Amp1346    Amp1540 MMPRNT1965    Amp1514    Amp1656 MMPRNT1273    Amp1334    Amp1552    Amp1347 
#16186      20172      22175      25128      28959      29312      29332      30298      30900      32043 





#####Mock community with singletons removed######

GLBRC018_OTU_fung.mock_f2=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2,seq_samp_type=="mock")
nsamples(GLBRC018_OTU_fung.mock_f2)
#12
GLBRC018_OTU_fung.mock_f2=prune_taxa(taxa_sums(GLBRC018_OTU_fung.mock_f2) > 0, GLBRC018_OTU_fung.mock_f2)
ntaxa(GLBRC018_OTU_fung.mock_f2)
#717


sum(otu_table(GLBRC018_OTU_fung.mock_f2))
#1567641

max(sample_sums(GLBRC018_OTU_fung.mock_f2))
#218266
min(sample_sums(GLBRC018_OTU_fung.mock_f2))
#84054

sort(sample_sums(GLBRC018_OTU_fung.mock_f2))
#Amp1281 Amp1591 Amp1360 Amp1495 Amp1936 Amp1584 Amp2102 Amp1487 Amp2070   Mock2 Amp2187   Mock1 
#84054   93328   97774  101500  106412  115315  122289  135653  138784  161299  192967  218266 


#Let's ordinate this 
nsamples(GLBRC018_OTU_fung.mock_f2)
#12
GLBRC018_OTU_fung.mock_f2_fil=subset_samples(GLBRC018_OTU_fung.mock_f2,sampleID_seq!="Amp1936"&sampleID_seq!="Amp1966")
nsamples(GLBRC018_OTU_fung.mock_f2_fil)
#11
GLBRC018_OTU_fung.mock_f2_fil_ord=ordinate(GLBRC018_OTU_fung.mock_f2_fil,method = "NMDS")
#Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data


plot_ordination(GLBRC018_OTU_fung.mock_f2_fil,GLBRC018_OTU_fung.mock_f2_fil_ord,color = "spiked_sample",label = "sampleID_seq")
#not great but the best I can do for now. 

#####Rarefy the community####

GLBRC018_OTU_fung_MMPRNT_mock_rar=rarefy_even_depth(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2,replace=F,sample.size=10000,rngseed=99)
#46OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(GLBRC018_OTU_fung_MMPRNT_mock_rar, file = here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_rar.RData"))



#####Mock community with Rarefy the community######

GLBRC018_OTU_fung_mock_rar=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_rar,seq_samp_type=="mock")
nsamples(GLBRC018_OTU_fung_mock_rar)
#12
GLBRC018_OTU_fung_mock_rar=prune_taxa(taxa_sums(GLBRC018_OTU_fung_mock_rar) > 0, GLBRC018_OTU_fung_mock_rar)
ntaxa(GLBRC018_OTU_fung_mock_rar)
#109


sum(otu_table(GLBRC018_OTU_fung_mock_rar))
#120000

max(sample_sums(GLBRC018_OTU_fung_mock_rar))
#10000
min(sample_sums(GLBRC018_OTU_fung_mock_rar))
#10000




#Let's ordinate this after removing the wonky samples
nsamples(GLBRC018_OTU_fung_mock_rar)
#12
GLBRC018_OTU_fung_mock_rar_fil=subset_samples(GLBRC018_OTU_fung_mock_rar,sampleID_seq!="Amp1936"&sampleID_seq!="Amp1966")
nsamples(GLBRC018_OTU_fung_mock_rar_fil)
#11
GLBRC018_OTU_fung_mock_rar_fil_ord=ordinate(GLBRC018_OTU_fung_mock_rar_fil,method = "NMDS")
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
# stress is (nearly) zero: you may have insufficient data


plot_ordination(GLBRC018_OTU_fung_mock_rar_fil,GLBRC018_OTU_fung_mock_rar_fil_ord,color = "spiked_sample",label = "sampleID_seq")
#not great but the best I can do for now. 

#####Remove non fungal reads####
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2)
#6559
sum(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))
#106642973
head(tax_table(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2))
GLBRC018_OTU_fung_MMPRNT_mock_fung=subset_taxa(GLBRC018_OTU_fung_MMPRNT_mock_fin_fil2,Domain=="d:Fungi")
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung)
#1088
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fung)
#3099
sum(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fung))
#74298553
GLBRC018_OTU_fung_MMPRNT_mock_fung=subset_taxa(GLBRC018_OTU_fung_MMPRNT_mock_fung,Class!="c:Malasseziomycetes")
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung)
#1088
ntaxa(GLBRC018_OTU_fung_MMPRNT_mock_fung)
#3099
sum(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fung))
#74298553
sort(sample_sums(GLBRC018_OTU_fung_MMPRNT_mock_fung))[1:20]


#####Rarefy the community####

GLBRC018_OTU_fung_MMPRNT_mock_fung_rar=rarefy_even_depth(GLBRC018_OTU_fung_MMPRNT_mock_fung,replace=F,sample.size=10000,rngseed=99)
#2 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
#Amp1343Amp1599

#9OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar, file = here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
#####END: Pre-processing Fungal community####





#####All Sites Analysis#####


#save(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar, file = here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)
#1109
head(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar))

#I am going to define All Sites as the JUly overlap in root and soil sampling

unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"&siteID=="LUX"))$collectionDate)
#"5/29/2018" "9/17/2018" "8/20/2018" "7/30/2018" "6/25/2018" "10/3/2018"
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"&siteID=="LC"))$collectionDate)
#"7/10/2018"


GLBRC018_OTU_bact_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,siteID!="LUX"|collectionDate=="7/30/2018")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites)
#595
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LUX"))$collectionDate)

GLBRC018_OTU_bact_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,siteID!="LC"|collectionDate=="7/10/2018")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites)
#405
unique(sample_data(subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LC"))$collectionDate)

#Let's use only G5 since there is overlap in sampling between roots and soil

GLBRC018_OTU_bact_MMPRNT_All_sites_G5=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites,plotType=="G5")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
#227

#save(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar, file = "here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#1086
head(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar))

#I am going to define All Sites as the JUly overlap in root and soil sampling

unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&siteID=="LUX"))$collectionDate)
#"5/29/2018" "9/17/2018" "8/20/2018" "7/30/2018" "6/25/2018" "10/3/2018"
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&siteID=="LC"))$collectionDate)
#"7/10/2018"


GLBRC018_OTU_fung_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,siteID!="LUX"|collectionDate=="7/30/2018")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites)
#584
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LUX"))$collectionDate)

GLBRC018_OTU_fung_MMPRNT_All_sites=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,siteID!="LC"|collectionDate=="7/10/2018")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites)
#402
unique(sample_data(subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,Root_soil=="Root"&siteID=="LC"))$collectionDate)

#Let's use only G5 since there is overlap in sampling between roots and soil

GLBRC018_OTU_fung_MMPRNT_All_sites_G5=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites,plotType=="G5")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)
#223





##Combined the roots and soils

#Let's look at a rough NMDS

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,method = "NMDS")
#*** Solution reached
#0.1171765  


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,method = "NMDS")
#*** Solution reached
#0.2460191 

(all_site_fert_bact_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_ord)+
    geom_point(aes(color=FertStatus, shape=Root_soil),size=4)+theme_bw()+ggtitle("Bacterial Communities")+
    scale_shape_manual(values=c(19,17),name=NULL)+scale_color_discrete(name=NULL)+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5))))

(all_site_compart_bact_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_ord)+
    geom_point(aes(color=siteID, shape=Root_soil),size=4)+theme_bw()+
    scale_shape_manual(values=c(19,17),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none"))

(all_site_fert_fung_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_ord)+
    geom_point(aes(color=FertStatus, shape=Root_soil),size=4)+theme_bw()+ggtitle("Fungal Communities")+
    scale_shape_manual(values=c(19,17),name=NULL)+scale_color_discrete(name=NULL)+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

(all_site_compart_fung_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_ord)+
    geom_point(aes(color=siteID, shape=Root_soil),size=4)+theme_bw()+
    scale_shape_manual(values=c(19,17),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.text = element_text(size = 16)))

plot_grid(all_site_fert_bact_p,all_site_fert_fung_p,all_site_compart_bact_p,all_site_compart_fung_p,ncol = 2,
          rel_widths = c(1,1),axis = "r",align = "v")

#Roots only 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,Root_soil=="Root")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root)
#113
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,Root_soil=="Root")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root)
#111

#Let's look at a rough NMDS

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root,method = "NMDS")
#*** Solution reached
#0.1226669   


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root,method = "NMDS")
#*** Solution reached
#0.2188881  



(all_site_compart_bact_root_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Root Bacterial Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5))))


(all_site_compart_fung_root_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Root Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

plot_grid(all_site_compart_bact_root_p,all_site_compart_fung_root_p,ncol = 2,
          rel_widths = c(1,1.3))
#1500*700


#Soil only 
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil=subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,Root_soil=="Soil")
nsamples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil)
#114
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,Root_soil=="Soil")
nsamples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil)
#112

#Let's look at a rough NMDS

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil,method = "NMDS")
#*** Solution reached
#0.1199976    


GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil,method = "NMDS")
#*** Solution reached
#0.1712595   



(all_site_compart_bact_soil_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Bacterial Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5))))


(all_site_compart_fung_soil_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil_ord)+
    geom_point(aes(color=siteID, shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    scale_color_brewer(palette = "Set1",name=NULL,labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhinelander"))+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

plot_grid(all_site_compart_bact_soil_p,all_site_compart_fung_soil_p,ncol = 2,
          rel_widths = c(1,1.3))

#1500*700




#
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_dis=distance(GLBRC018_OTU_bact_MMPRNT_All_sites_G5,method = "bray")
#write.csv(as.matrix(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_dis), here::here("R_files","distance_files","bact_GLBRC018_OTU_bact_MMPRNT_All_sites_G5_dis.csv"))


GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
#write.csv(as.matrix(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map), here::here("R_files","distance_files","bact_GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map.csv"))

#adonis(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_dis~GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map$Root_soil*
#         GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map$siteID*GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map$FertStatus,permutations = 9999)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_dis=distance(GLBRC018_OTU_fung_MMPRNT_All_sites_G5,method = "bray")
#write.csv(as.matrix(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_dis), here::here("R_files","distance_files","fungi_GLBRC018_OTU_fung_MMPRNT_All_sites_G5_dis.csv"))

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)
#write.csv(as.matrix(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map), here::here("R_files","distance_files","fungi_GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map.csv"))


#####Beta-disp####



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

site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("Lux Arbor","Lake City","Escanaba", "Hancock","Rhinlander")



#Betadisp



(bact_betadisp_p=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis)+
  geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),fill=FertStatus))+
    geom_text(data = bact_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Betadispersion",limits = c(0.22,0.57))+theme(axis.text.x = element_blank()))

(fung_betadisp_p=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis)+
    geom_boxplot(aes(y=betadisp, x=factor(siteID,levels = site_order),fill=FertStatus))+
    geom_text(data = fung_betadisp_max_disp_letters, aes(x=siteID, y = 0.01 + max_value, label = sig_let), vjust=0, size=10)+
  facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
  ylab("Betadispersion")+theme(strip.text = element_blank(),strip.background = element_blank()))


#1700x1200

plot_grid(bact_betadisp_p,fung_betadisp_p,ncol = 1,labels = c('a)', 'b)'), label_size = 30)


#Stats 

#Roots Bacteria
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_root=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis,Root_soil=="Root")

root_betadisp_all_site_bact_mod=lmer(log(betadisp)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_root)
plot(root_betadisp_all_site_bact_mod)
hist(resid(root_betadisp_all_site_bact_mod))
qqPlot(resid(root_betadisp_all_site_bact_mod))
shapiro.test(resid(root_betadisp_all_site_bact_mod))
#W = 0.98237, p-value = 0.143


anova(root_betadisp_all_site_bact_mod)
#siteID            0.37223 0.093057     4 101.18  6.3060 0.0001428 ***
#FertStatus        0.01573 0.015731     1 100.31  1.0660 0.3043382    
#siteID:FertStatus 0.06075 0.015187     4 100.31  1.0292 0.3960468   



emmeans(root_betadisp_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate     SE    df t.ratio p.value
#ESC - HAN   0.0412 0.0386 103  1.068  0.8222 
#ESC - LC   -0.0910 0.0351 100 -2.594  0.0791 
#ESC - LUX   0.0166 0.0355 100  0.467  0.9901 
#ESC - RHN   0.0770 0.0351 100  2.197  0.1894 
#HAN - LC   -0.1322 0.0386 103 -3.427  0.0077 
#HAN - LUX  -0.0247 0.0388 103 -0.636  0.9688 
#HAN - RHN   0.0358 0.0386 103  0.928  0.8853 
#LC - LUX    0.1075 0.0355 100  3.031  0.0252 
#LC - RHN    0.1680 0.0351 100  4.791  0.0001 
#LUX - RHN   0.0605 0.0355 100  1.705  0.4359 


root_bact_disp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                         sig_let=c("ab","a","b","a","a"),
                         Root_soil=rep("Root"))


#Soil
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_soil=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis,Root_soil=="Soil")

soil_betadisp_all_site_bact_mod=lmer((betadisp)^-3~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis_soil)
plot(soil_betadisp_all_site_bact_mod)
hist(resid(soil_betadisp_all_site_bact_mod))
qqPlot(resid(soil_betadisp_all_site_bact_mod))
shapiro.test(resid(soil_betadisp_all_site_bact_mod))
#W = 0.97723, p-value = 0.04884


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


soil_bact_disp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("A","B","B","BC","C"),
                                  Root_soil=rep("Soil"))

root_soil_bact_disp_letters=rbind(root_bact_disp_letters,soil_bact_disp_letters)
bact_betadisp_max=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_betaDis%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(betadisp))

bact_betadisp_max_disp_letters=merge(root_soil_bact_disp_letters,bact_betadisp_max,by=c("siteID","Root_soil"))



#Roots
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_root=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis,Root_soil=="Root")

root_betadisp_all_site_fung_mod=lmer(sqrt(betadisp)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_root)
plot(root_betadisp_all_site_fung_mod)
hist(resid(root_betadisp_all_site_fung_mod))
qqPlot(resid(root_betadisp_all_site_fung_mod))
shapiro.test(resid(root_betadisp_all_site_fung_mod))
#W = 0.98262, p-value = 0.1585


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

root_fung_disp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("ab","ac","b","c","a"),
                                  Root_soil=rep("Root"))

#Soil
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_soil=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis,Root_soil=="Soil")

soil_betadisp_all_site_fung_mod=lmer(betadisp~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis_soil)
plot(soil_betadisp_all_site_fung_mod)
hist(resid(soil_betadisp_all_site_fung_mod))
qqPlot(resid(soil_betadisp_all_site_fung_mod))
shapiro.test(resid(soil_betadisp_all_site_fung_mod))
#W = 0.97726, p-value = 0.05284


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


soil_fung_disp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("A","B","C","BC","C"),
                                  Root_soil=rep("Soil"))
root_soil_fung_disp_letters=rbind(root_fung_disp_letters,soil_fung_disp_letters)

fung_betadisp_max=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_betaDis%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(betadisp))

fung_betadisp_max_disp_letters=merge(root_soil_fung_disp_letters,fung_betadisp_max,by=c("siteID","Root_soil"))


#####Diversity####

#Diversity 

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div=estimate_richness(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div=merge(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map,by="row.names")
head(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div)
unique(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div$siteID)

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div=estimate_richness(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div=merge(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map,by="row.names")
head(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div)
unique(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div$siteID)



site_order=c("LUX","LC","ESC", "HAN","RHN")
site_labels=c("Lux Arbor","Lake City","Escanaba", "Hancock","Rhinlander")



####Richness####


(bact_rich_p=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div)+
    geom_boxplot(aes(y=Observed, x=factor(siteID,levels = site_order),fill=FertStatus))+
    geom_text(data = bact_rich_max_disp_letters, aes(x=siteID, y = 10 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Richness",limits = c(700,2800))+theme(axis.text.x = element_blank()))


pairwise_fert_fung_rich=data.frame(y_bot=c(0,0,0,225,0,0,0,0,0,0),
                                   x_min=c(0,0,0,3.85,0,0,0,0,0,0),
                                   x_max=c(0,0,0,4.15,0,0,0,0,0,0),
                                   annot_text=c("","","","*","","","","","",""),
                                   Root_soil=c(Root_soil=c(rep("Root",5),rep("Soil",5))))

(fung_rich_p=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,aes(y=Observed, x=factor(siteID,levels = site_order)))+
    geom_boxplot(aes(y=Observed, x=factor(siteID,levels = site_order),fill=FertStatus))+
    geom_text(data = fung_rich_max_rich_letters, aes(x=siteID, y = 10 + max_value, label = sig_let), vjust=0, size=10)+
    geom_signif(data=pairwise_fert_fung_rich,aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Richness",limits = c(98,520))+theme(strip.text = element_blank(),strip.background = element_blank()))






#1700x1200

plot_grid(bact_rich_p,fung_rich_p,ncol = 1,labels = c('a)', 'b)'), label_size = 30)



#Stats 

#Roots
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_rich_all_site_bact_mod=lmer((Observed)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root)
plot(root_rich_all_site_bact_mod)
hist(resid(root_rich_all_site_bact_mod))
qqPlot(resid(root_rich_all_site_bact_mod))
shapiro.test(resid(root_rich_all_site_bact_mod))
#W = 0.98842, p-value = 0.4491


anova(root_rich_all_site_bact_mod)
#siteID            2843185  710796     4 100.425 40.5931 < 2.2e-16 ***
#FertStatus         287458  287458     1  99.979 16.4165 0.0001006 ***
#siteID:FertStatus  214424   53606     4  99.980  3.0614 0.0200188 *   

emmeans(root_rich_all_site_bact_mod,pairwise~FertStatus|siteID)

#$contrasts
#siteID = LC:
#contrast      estimate   SE  df t.ratio p.value
#Fert - Unfert  -192.33 54.0 100 -3.560  0.0006 
#siteID = RHN:
#contrast      estimate   SE  df t.ratio p.value
#Fert - Unfert  -202.67 54.0 100 -3.752  0.0003 

emmeans(root_rich_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN    216.4 42.1 101  5.145  <.0001 
#ESC - LC     100.0 38.2 100  2.617  0.0749 
#ESC - LUX    -44.5 38.6 100 -1.151  0.7788 
#ESC - RHN    387.0 38.2 100 10.132  <.0001 
#HAN - LC    -116.4 42.1 101 -2.768  0.0512 
#HAN - LUX   -260.9 42.3 101 -6.174  <.0001 
#HAN - RHN    170.7 42.1 101  4.057  0.0009 
#LC - LUX    -144.4 38.6 100 -3.738  0.0028 
#LC - RHN     287.1 38.2 100  7.515  <.0001 
#LUX - RHN    431.5 38.6 100 11.166  <.0001 

root_bact_rich_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("ac","b","ab","c","d"),
                                  Root_soil=rep("Root"))

#Soil
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_rich_all_site_bact_mod=lmer(Observed~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil)
plot(soil_rich_all_site_bact_mod)
hist(resid(soil_rich_all_site_bact_mod))
qqPlot(resid(soil_rich_all_site_bact_mod))
shapiro.test(resid(soil_rich_all_site_bact_mod))
#W = 0.99289, p-value = 0.8277

anova(soil_rich_all_site_bact_mod)
#siteID            1876754  469188     4 101.32 33.2386 <2e-16 ***



emmeans(soil_rich_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN     53.7 37.7 102  1.422  0.6148 
#ESC - LC     -31.8 34.3 101 -0.926  0.8863 
#ESC - LUX    -64.2 34.3 101 -1.872  0.3390 
#ESC - RHN    286.6 34.3 101  8.356  <.0001 
#HAN - LC     -85.4 37.7 102 -2.264  0.1653 
#HAN - LUX   -117.9 37.7 102 -3.124  0.0193 
#HAN - RHN    232.9 37.7 102  6.174  <.0001 
#LC - LUX     -32.5 34.3 101 -0.946  0.8780 
#LC - RHN     318.3 34.3 101  9.282  <.0001 
#LUX - RHN    350.8 34.3 101 10.228  <.0001 


soil_bact_rich_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("AB","A","AB","B","C"),
                                  Root_soil=rep("Soil"))

root_soil_bact_rich_letters=rbind(root_bact_rich_letters,soil_bact_rich_letters)
bact_rich_max=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(Observed))

bact_rich_max_disp_letters=merge(root_soil_bact_rich_letters,bact_rich_max,by=c("siteID","Root_soil"))


#Roots
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_rich_all_site_fung_mod=lmer(sqrt(Observed)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root)
plot(root_rich_all_site_fung_mod)
hist(resid(root_rich_all_site_fung_mod))
qqPlot(resid(root_rich_all_site_fung_mod))
shapiro.test(resid(root_rich_all_site_fung_mod))
#W = 0.98082, p-value = 0.1111


anova(root_rich_all_site_fung_mod)
#siteID            61.971 15.4928     4 98.607  9.4323 1.692e-06 ***
#FertStatus         1.295  1.2949     1 98.121  0.7883    0.3768    
#siteID:FertStatus 11.931  2.9828     4 98.122  1.8160    0.1318   

emmeans(root_rich_all_site_fung_mod,pairwise~FertStatus|siteID)

#$contrasts
#siteID = HAN:
#contrast      estimate    SE   df t.ratio p.value
#Fert - Unfert  -1.3057 0.604 98.0 -2.161  0.0331 



emmeans(root_rich_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN    1.411 0.409 100.8  3.449  0.0072 
#ESC - LC     0.291 0.383  98.4  0.759  0.9415 
#ESC - LUX   -0.976 0.374  98.2 -2.609  0.0765 
#ESC - RHN    0.449 0.374  98.2  1.199  0.7519 
#HAN - LC    -1.120 0.416 100.9 -2.691  0.0623 
#HAN - LUX   -2.387 0.407 101.0 -5.865  <.0001 
#HAN - RHN   -0.962 0.407 101.0 -2.363  0.1341 
#LC - LUX    -1.267 0.379  98.2 -3.348  0.0099 
#LC - RHN     0.158 0.379  98.2  0.417  0.9935 
#LUX - RHN    1.425 0.370  98.0  3.852  0.0019


root_fung_rich_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("ac","b","ab","c","ab"),
                                  Root_soil=rep("Root"))


#Soil
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_rich_all_site_fung_mod=lmer(Observed~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil)
plot(soil_rich_all_site_fung_mod)
hist(resid(soil_rich_all_site_fung_mod))
qqPlot(resid(soil_rich_all_site_fung_mod))
shapiro.test(resid(soil_rich_all_site_fung_mod))
#W = 0.99242, p-value = 0.3085


anova(soil_rich_all_site_fung_mod)
#siteID            109725 27431.3     4   102 29.8821 <2e-16 ***



emmeans(soil_rich_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN   16.153 9.61 102.0  1.680  0.4507 
#ESC - LC     0.542 8.75  99.0  0.062  1.0000 
#ESC - LUX  -71.462 8.85  99.2 -8.075  <.0001 
#ESC - RHN    5.746 8.85  99.2  0.649  0.9664 
#HAN - LC   -15.611 9.61 102.0 -1.624  0.4858 
#HAN - LUX  -87.615 9.72 101.9 -9.012  <.0001 
#HAN - RHN  -10.407 9.66 101.8 -1.077  0.8179 
#LC - LUX   -72.004 8.85  99.2 -8.137  <.0001 
#LC - RHN     5.205 8.85  99.2  0.588  0.9766 
#LUX - RHN   77.208 8.95  99.4  8.624  <.0001 



soil_fung_rich_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("A","A","A","B","A"),
                                  Root_soil=rep("Soil"))
root_soil_fung_rich_letters=rbind(root_fung_rich_letters,soil_fung_rich_letters)

fung_rich_max=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(Observed))

fung_rich_max_rich_letters=merge(root_soil_fung_rich_letters,fung_rich_max,by=c("siteID","Root_soil"))


###InvSimpson####


test_test_test=data.frame(y_bot=c(0,0,0,250,0,0,0,0,0,0),
                                   x_min=c(0.85,1.85,2.85,3.85,4.85,0.85,1.85,2.85,3.85,4.85),
                                   x_max=c(1.15,2.15,3.15,4.15,5.15,1.15,2.15,3.15,4.15,5.15),annot_text=c("","**","#","","**","","","","",""),
                                   annot_text=c("","","","*","","","","","",""),
                                   siteID=c("ESC","HAN","LC","LUX","RHN"),
                                   Root_soil=c(Root_soil=c(rep("Root",5),rep("Soil",5))))

pairwise_fert_bact_invSimp=data.frame(y_bot=c(0,95,100,0,55,0,0,0,0,0),
                                   x_min=c(0,1.85,2.85,0,4.85,0,0,0,0,0),
                                   x_max=c(0,2.15,3.15,0,5.15,0,0,0,0,0),
                                   annot_text=c(NA,"**","#",NA," ** ",NA,NA,NA,NA,NA),
                                   Root_soil=c(Root_soil=c(rep("Root",5),rep("Soil",5))))

(bact_invSimp_p=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div)+
   geom_boxplot(aes(y=InvSimpson, x=factor(siteID,levels = site_order),fill=FertStatus))+
   geom_text(data = bact_invSimp_max_disp_letters, aes(x=siteID, y = 40 + max_value, label = sig_let), vjust=0, size=10)+
    geom_signif(data=pairwise_fert_bact_invSimp,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
   facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
   scale_y_continuous(name = "Inverse Simpson",limits = c(0,520))+theme(axis.text.x = element_blank()))

pairwise_fert_fung_invsimp=data.frame(y_bot=c(0,0,17,0,0,0,0,0,0,0),
                          x_min=c(0,0,2.85,0,0,0,0,0,0,0),
                          x_max=c(0,0,3.15,0,0,0,0,0,0,0),
                          annot_text=c("","","**","","","","","","",""),
                          Root_soil=c(Root_soil=c(rep("Root",5),rep("Soil",5))))


(fung_invSimp_p=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,aes(y=Observed, x=factor(siteID,levels = site_order)))+
    geom_boxplot(aes(y=InvSimpson, x=factor(siteID,levels = site_order),fill=FertStatus))+
    geom_text(data = fung_invsimp_max_invsimp_letters, aes(x=siteID, y = 1 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    geom_signif(data=pairwise_fert_fung_invsimp,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    scale_y_continuous(name = "Inverse Simpson",limits = c(0,57))+theme(strip.text = element_blank(),strip.background = element_blank()))


#1800x1200

plot_grid(bact_invSimp_p,fung_invSimp_p,ncol = 1,labels = c('a)', 'b)'), label_size = 30)



#Roots
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_invsimp_all_site_bact_mod=lmer(log(InvSimpson)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root)
plot(root_invsimp_all_site_bact_mod)
hist(resid(root_invsimp_all_site_bact_mod))
qqPlot(resid(root_invsimp_all_site_bact_mod))
shapiro.test(resid(root_invsimp_all_site_bact_mod))
#W = 0.98644, p-value = 0.3154


anova(root_invsimp_all_site_bact_mod)
#siteID            9.6381 2.40952     4 100.846  15.613 5.656e-10 ***
#FertStatus        1.7970 1.79699     1  99.992  11.644  0.000931 ***
#siteID:FertStatus 0.8544 0.21359     4  99.992   1.384  0.244919 



emmeans(root_invsimp_all_site_bact_mod,pairwise~FertStatus|siteID)
#siteID = LC:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.4439 0.160 100 -2.768  0.0067 

#siteID = RHN:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.4283 0.160 100 -2.671  0.0088 

emmeans(root_invsimp_all_site_bact_mod,pairwise~siteID)
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN   0.2171 0.125 103  1.738  0.4157 
#ESC - LC   -0.0906 0.113 100 -0.799  0.9303 
#ESC - LUX  -0.1890 0.115 100 -1.647  0.4714 
#ESC - RHN   0.6152 0.113 100  5.425  <.0001 
#HAN - LC   -0.3077 0.125 103 -2.463  0.1072 
#HAN - LUX  -0.4060 0.125 102 -3.236  0.0138 
#HAN - RHN   0.3982 0.125 103  3.188  0.0160 
#LC - LUX   -0.0983 0.115 100 -0.857  0.9117 
#LC - RHN    0.7058 0.113 100  6.224  <.0001 
#LUX - RHN   0.8042 0.115 100  7.009  <.0001

root_bact_invSimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                     sig_let=c("ab","a","a","b","c"),
                                     Root_soil=rep("Root"))

#Soil
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_invsimp_all_site_bact_mod=lmer(InvSimpson~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil)
plot(soil_invsimp_all_site_bact_mod)
hist(resid(soil_invsimp_all_site_bact_mod))
qqPlot(resid(soil_invsimp_all_site_bact_mod))
shapiro.test(resid(soil_invsimp_all_site_bact_mod))
#W = 0.98546, p-value = 0.2568


anova(soil_invsimp_all_site_bact_mod)
#siteID            136322   34080     4 101.45  7.9600 1.273e-05 ***


emmeans(soil_invsimp_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN     27.2 20.8 102  1.306  0.6880 
#ESC - LC      10.5 18.9 101  0.554  0.9812 
#ESC - LUX    -36.2 18.9 101 -1.916  0.3154 
#ESC - RHN     67.2 18.9 101  3.559  0.0051 
#HAN - LC     -16.7 20.8 102 -0.803  0.9292 
#HAN - LUX    -63.4 20.8 102 -3.047  0.0240 
#HAN - RHN     40.1 20.8 102  1.928  0.3094 
#LC - LUX     -46.7 18.9 101 -2.470  0.1057 
#LC - RHN      56.8 18.9 101  3.005  0.0271 
#LUX - RHN    103.4 18.9 101  5.476  <.0001 

soil_bact_invSimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("AB","AC","AB","B","C"),
                                  Root_soil=rep("Soil"))



root_soil_bact_invSimp_letters=rbind(root_bact_invSimp_letters,soil_bact_invSimp_letters)
bact_invSimp_max=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(InvSimpson))

bact_invSimp_max_disp_letters=merge(root_soil_bact_invSimp_letters,bact_invSimp_max,by=c("siteID","Root_soil"))

#Roots
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_invsimp_all_site_fung_mod=lmer(log(InvSimpson)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root)
plot(root_invsimp_all_site_fung_mod)
hist(resid(root_invsimp_all_site_fung_mod))
qqPlot(resid(root_invsimp_all_site_fung_mod))
shapiro.test(resid(root_invsimp_all_site_fung_mod))
#W = 0.98796, p-value = 0.4278


anova(root_invsimp_all_site_fung_mod)
#siteID            2.2116 0.55291     4 98.721  2.7730 0.031222 * 
#FertStatus        0.1187 0.11875     1 98.285  0.5956 0.442128   
#siteID:FertStatus 3.2446 0.81116     4 98.286  4.0683 0.004287 **

emmeans(root_invsimp_all_site_fung_mod,pairwise~FertStatus|siteID)


#$contrasts
#siteID = ESC:
#contrast      estimate    SE   df t.ratio p.value
#Fert - Unfert   0.6009 0.187 98.3  3.221  0.0017 

emmeans(root_invsimp_all_site_fung_mod,pairwise~siteID)
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN -0.17822 0.142 100.8 -1.251  0.7215 
#ESC - LC   0.07710 0.133  98.5  0.578  0.9780 
#ESC - LUX  0.00213 0.130  98.2  0.016  1.0000 
#ESC - RHN -0.30240 0.130  98.2 -2.319  0.1477 
#HAN - LC   0.25533 0.145 100.9  1.762  0.4017 
#HAN - LUX  0.18035 0.142 101.0  1.272  0.7089 
#HAN - RHN -0.12418 0.142 101.0 -0.876  0.9051 
#LC - LUX  -0.07497 0.132  98.2 -0.569  0.9793 
#LC - RHN  -0.37951 0.132  98.2 -2.878  0.0386 
#LUX - RHN -0.30453 0.129  98.0 -2.363  0.1345

root_fung_invsimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                     sig_let=c("ab","ab","a","ab","b"),
                                     Root_soil=rep("Root"))

#Soil
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_invsimp_all_site_fung_mod=lmer(InvSimpson~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil)
plot(soil_invsimp_all_site_fung_mod)
hist(resid(soil_invsimp_all_site_fung_mod))
qqPlot(resid(soil_invsimp_all_site_fung_mod))
shapiro.test(resid(soil_invsimp_all_site_fung_mod))
#W = 0.99223, p-value = 0.7821


anova(soil_invsimp_all_site_fung_mod)
#siteID            970.49  242.62     4 99.823  3.2402 0.01522 *
#FertStatus          2.21    2.21     1 99.176  0.0295 0.86395  
#siteID:FertStatus 562.79  140.70     4 99.159  1.8790 0.12006 


emmeans(soil_invsimp_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN    2.007 2.75 101.1  0.729  0.9492 
#ESC - LC     8.444 2.50  99.0  3.380  0.0090 
#ESC - LUX    2.198 2.53  99.0  0.870  0.9073 
#ESC - RHN    4.228 2.53  99.1  1.673  0.4552 
#HAN - LC     6.437 2.75 101.1  2.339  0.1412 
#HAN - LUX    0.191 2.78 101.2  0.069  1.0000 
#HAN - RHN    2.222 2.76 100.6  0.804  0.9289 
#LC - LUX    -6.246 2.53  99.0 -2.472  0.1055 
#LC - RHN    -4.215 2.53  99.1 -1.668  0.4584 
#LUX - RHN    2.030 2.56  99.1  0.794  0.9318 


soil_fung_invsimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                  sig_let=c("A","AB","B","AB","AB"),
                                  Root_soil=rep("Soil"))
root_soil_fung_invsimp_letters=rbind(root_fung_invsimp_letters,soil_fung_invsimp_letters)

fung_invsimp_max=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(InvSimpson))

fung_invsimp_max_invsimp_letters=merge(root_soil_fung_invsimp_letters,fung_invsimp_max,by=c("siteID","Root_soil"))




#####Lux Arbor analyses#####

#save(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar, file = here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)
#1109
head(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar))


#Let's focus on G5 


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,plotType=="G5"&siteID=="LUX")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
#496
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5) > 0, GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)

#save(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar, file = here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#1086
head(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar))


#Let's focus on G5 


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,plotType=="G5"&siteID=="LUX")
nsamples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
#490


#Let's ordinate

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#20: scale factor of the gradient < sfgrmin
#0.06852062  


(LUX_fert_bact_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_ord)+
    geom_point(aes(color=FertStatus, shape=Root_soil),size=4)+theme_bw())

(LUX_compart_bact_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=Root_soil),size=4)+geom_text(aes(label=SampleID))+
    theme_bw())

plot_grid(LUX_fert_p,LUX_compart_p,ncol = 2)

#Amp129 is an outlier so I am going to remove it for now

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,SampleID!="Amp129")

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,method = "NMDS")
#*** Solution reached
#0.06875135   

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#3: no. of iterations >= maxit
#6: stress ratio > sratmax
#11: scale factor of the gradient < sfgrmin
#0.1462961 


(LUX_fert_bact_p2=plot_ordination(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil_ord)+
    geom_point(aes(color=FertStatus, shape=Root_soil),size=4)+theme_bw()+ggtitle("Bacterial Communities")+
    scale_shape_manual(values=c(19,17),name=NULL)+scale_color_discrete(name=NULL)+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5))))

(LUX_compart_bact_p2=plot_ordination(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=Root_soil),size=4)+theme_bw()+
    scale_shape_manual(values=c(19,17),name=NULL)+
    theme(legend.position = "none"))

(LUX_fert_fung_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_ord)+
    geom_point(aes(color=FertStatus, shape=Root_soil),size=4)+theme_bw()+ggtitle("Fungal Communities")+
    scale_shape_manual(values=c(19,17),name=NULL)+scale_color_discrete(name=NULL)+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

(LUX_compart_fung_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=Root_soil),size=4)+theme_bw()+
    scale_shape_manual(values=c(19,17),name=NULL)+
    theme(legend.text = element_text(size = 16)))

plot_grid(LUX_fert_bact_p2,LUX_fert_fung_p,LUX_compart_bact_p2,LUX_compart_fung_p,ncol = 2,
          rel_widths = c(1,1),axis = "r",align = "v")
#1400*1200

#Roots only 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_root=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,Root_soil=="Root")
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root)
#142
GLBRC018_OTU_fung_MMPRNT_LUX_G5_root=subset_samples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,Root_soil=="Root")
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root)
#141

#Let's look at a rough NMDS

GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root,method = "NMDS")
#*** Solution reached
#0.2373164   


GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root,method = "NMDS")
#*** No convergence -- monoMDS stopping criteria:
#1: no. of iterations >= maxit
#19: stress ratio > sratmax
#0.2194756   



(LUX_compart_bact_root_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_LUX_G5_root,GLBRC018_OTU_bact_MMPRNT_LUX_G5_root_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Root Bacterial Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5))))


(LUX_compart_fung_root_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_LUX_G5_root,GLBRC018_OTU_fung_MMPRNT_LUX_G5_root_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Root Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

plot_grid(LUX_compart_bact_root_p,LUX_compart_fung_root_p,ncol = 2,
          rel_widths = c(1,1.2))
#1500*700


#Soil only 
GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,Root_soil=="Soil")
nsamples(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil)
#353
GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,Root_soil=="Soil")
nsamples(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil)
#349

#Let's look at a rough NMDS

GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_ord=ordinate(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil,method = "NMDS")
#*** Solution reached
#0.2375187     


GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_ord=ordinate(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil,method = "NMDS")
#*** Solution reached
#0.1959129    



(LUX_compart_bact_soil_p=plot_ordination(GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil,GLBRC018_OTU_bact_MMPRNT_LUX_G5_soil_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Bacterial Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    theme(legend.position = "none",plot.title = (element_text(size = 30,hjust = 0.5))))


(LUX_compart_fung_soil_p=plot_ordination(GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil,GLBRC018_OTU_fung_MMPRNT_LUX_G5_soil_ord)+
    geom_point(aes(color=as.Date(collectionDate, format="%m/%d/%Y"), shape=FertStatus,fill=FertStatus),size=4,stroke=3)+theme_bw()+ggtitle("Soil Fungal Communities")+
    scale_shape_manual(values=c(19,24),name=NULL)+scale_fill_manual(values = c("black","white"),name=NULL)+
    theme(legend.text = element_text(size = 16),plot.title = (element_text(size = 30,hjust = 0.5))))

plot_grid(LUX_compart_bact_soil_p,LUX_compart_fung_soil_p,ncol = 2,
          rel_widths = c(1,1.2))

#1500*700





GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis=distance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil,method = "bray")
#write.csv(as.matrix(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis), here::here("R_files","distance_files","bact_GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis.csv"))
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
#write.csv(as.matrix(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map), here::here("R_files","distance_files","bact_GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map.csv"))


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis=distance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,method = "bray")
#write.csv(as.matrix(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis), here::here("R_files","distance_files","fungi_GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis.csv"))
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
#write.csv(as.matrix(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map), here::here("R_files","distance_files","fungi_GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map.csv"))


####Beta dispersion####

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betamod=betadisper(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis, with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,interaction(Root_soil,collectionDate,FertStatus)))
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis=as.data.frame(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betamod$distances)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,by="row.names")
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis)[colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis)==
                                                   "GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betamod$distances"]="betadisp"


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betamod=betadisper(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis, with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map,interaction(Root_soil,collectionDate,FertStatus)))
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis=as.data.frame(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betamod$distances)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map,by="row.names")
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis)[colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis)==
                                                        "GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betamod$distances"]="betadisp"


#betadisp
(lux_bact_betadisp_p=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,color=FertStatus))+
  geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
  facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+ylab("Betadispersion")+
  theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(axis.text.x = element_blank()))


(lux_fung_betadisp_p=ggplot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_betaDis, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=betadisp,color=FertStatus))+
  geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
  facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+ylab("Betadispersion")+
  theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(strip.text = element_blank(),strip.background = element_blank()))

#1700x1200

plot_grid(lux_bact_betadisp_p,lux_fung_betadisp_p,ncol = 1,labels = c('a)', 'b)'), label_size = 30)

#Stats

#Roots
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis_root=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Root")

root_betadisp_LUX_bact_mod=lmer(log(betadisp)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis_root)
plot(root_betadisp_LUX_bact_mod)
hist(resid(root_betadisp_LUX_bact_mod))
qqPlot(resid(root_betadisp_LUX_bact_mod))
shapiro.test(resid(root_betadisp_LUX_bact_mod))
#W = 0.98018, p-value = 0.03734


anova(root_betadisp_LUX_bact_mod)
#collectionDate            0.219132 0.043826     5 127.04  2.9450 0.01502 *
#FertStatus                0.000021 0.000021     1 127.03  0.0014 0.97036  
#collectionDate:FertStatus 0.088110 0.017622     5 127.04  1.1841 0.32062  

emmeans(root_betadisp_LUX_bact_mod,pairwise~FertStatus|collectionDate)

emmeans(root_betadisp_LUX_bact_mod,pairwise~collectionDate)

#Soils
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_soil=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_betaDis,Root_soil=="Soil")

soil_betadisp_LUX_bact_mod=lmer((betadisp)^-1~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_soil)
plot(soil_betadisp_LUX_bact_mod)
hist(resid(soil_betadisp_LUX_bact_mod))
qqPlot(resid(soil_betadisp_LUX_bact_mod))
shapiro.test(resid(soil_betadisp_LUX_bact_mod))
#W = 0.9939, p-value = 0.1675



anova(soil_betadisp_LUX_bact_mod, type = 3)
#collectionDate            3.0920 0.22086    14   320  2.8414 0.0004863 ***
#FertStatus                0.2301 0.23010     1   320  2.9603 0.0862989 .  
#collectionDate:FertStatus 1.1817 0.08441    14   320  1.0860 0.3693572    

emmeans(soil_betadisp_LUX_bact_mod,pairwise~FertStatus|collectionDate)






####Diversity ####

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div=estimate_richness(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map,by="row.names")
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div)
unique(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div$siteID)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div=estimate_richness(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map,by="row.names")
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div)
unique(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div$siteID)

#####Richness####

pairwise_fert_bact_rich_lux=data.frame(y_bot=c(0,0,2310,2500),
                                      x_min=as.Date(c(NA,NA,"4/25/2018","8/15/2018"), format="%m/%d/%Y"),
                                      x_max=as.Date(c(NA,NA,"5/5/2018","8/25/2018"), format="%m/%d/%Y"),
                                      annot_text=c("","","#","**"),
                                      Root_soil=c(Root_soil=c(rep("Root",2),rep("Soil",2))))

(lux_bact_rich_p=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=Observed))+
  geom_boxplot(aes(group=interaction(collectionDate,FertStatus),color=FertStatus),position = position_dodge(8))+
    geom_signif(data=pairwise_fert_bact_rich_lux,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
  facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+ylab("Richness")+
  theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(axis.text.x = element_blank()))
#

(lux_fung_rich_p=ggplot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=Observed,color=FertStatus))+
  geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
  facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+ylab("Richness")+
  theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(strip.text = element_blank(),strip.background = element_blank()))

#1700x1200

plot_grid(lux_bact_rich_p,lux_fung_rich_p,ncol = 1,labels = c('a)', 'b)'), label_size = 30)

#Stats

#Roots
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_rich_LUX_bact_mod=lmer(Observed~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root)
plot(root_rich_LUX_bact_mod)
hist(resid(root_rich_LUX_bact_mod))
qqPlot(resid(root_rich_LUX_bact_mod))
shapiro.test(resid(root_rich_LUX_bact_mod))
#W = 0.99436, p-value = 0.8554


anova(root_rich_LUX_bact_mod)
#collectionDate            524869  104974     5 127.05  5.4371 0.0001444 ***
#FertStatus                 46751   46751     1 127.04  2.4215 0.1221689    
#collectionDate:FertStatus 120034   24007     5 127.05  1.2434 0.2927484    

emmeans(root_rich_LUX_bact_mod,pairwise~FertStatus|collectionDate)
#


#Soils
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_rich_LUX_bact_mod=lmer((Observed)^2~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_rich_LUX_bact_mod)
hist(resid(soil_rich_LUX_bact_mod))
qqPlot(resid(soil_rich_LUX_bact_mod))
shapiro.test(resid(soil_rich_LUX_bact_mod))
#W = 0.997, p-value = 0.7637



anova(soil_rich_LUX_bact_mod)
#collectionDate            9.2361e+13 6.5972e+12    14 3.3998e+10 17.6458 < 2e-16 ***
#FertStatus                1.5374e+12 1.5374e+12     1 8.6013e+11  4.1121 0.04258 *  
#collectionDate:FertStatus 7.6441e+12 5.4601e+11    14 9.5095e+10  1.4604 0.11670   

emmeans(soil_rich_LUX_bact_mod,pairwise~FertStatus|collectionDate)
#collectionDate = 8/20/2018:
#contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert   782623 260722 320  3.002  0.0029 
#collectionDate = 4/30/2018:
 # contrast      estimate     SE  df t.ratio p.value
#Fert - Unfert  -464653 255285 320 -1.820  0.0697 




#Roots
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_rich_LUX_mod=lmer(Observed~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root)
plot(root_rich_LUX_mod)
hist(resid(root_rich_LUX_mod))
qqPlot(resid(root_rich_LUX_mod))
shapiro.test(resid(root_rich_LUX_mod))
#W = 0.99478, p-value = 0.8941


anova(root_rich_LUX_mod)
#collectionDate            134575 26915.1     5 126.02 23.5489 < 2e-16 ***
#FertStatus                  4088  4087.9     1 126.03  3.5767 0.06089 .  
#collectionDate:FertStatus   5012  1002.4     5 126.02  0.8771 0.49861   

emmeans(root_rich_LUX_mod,pairwise~FertStatus|collectionDate)
#collectionDate = 10/3/2018:
#contrast      estimate   SE  df t.ratio p.value
#Fert - Unfert   -31.83 13.8 126 -2.306  0.0227  


#Soils
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_rich_LUX_mod=lmer((Observed)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_rich_LUX_mod)
hist(resid(soil_rich_LUX_mod))
qqPlot(resid(soil_rich_LUX_mod))
shapiro.test(resid(soil_rich_LUX_mod))
#W = 0.99539, p-value = 0.3919



anova(soil_rich_LUX_mod, type = 3)
#collectionDate            175880 12562.9    14 316.01  8.5704 9.729e-16 ***
#FertStatus                   731   730.9     1 316.01  0.4986    0.4806    
#collectionDate:FertStatus  21023  1501.7    14 316.01  1.0244    0.4283  

emmeans(soil_rich_LUX_mod,pairwise~FertStatus|collectionDate)




####Inverse Simpson####
(lux_bact_invsimp_p=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=InvSimpson,color=FertStatus))+
  geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
  facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+xlab("Inverse Simpson")+
  theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(axis.text.x = element_blank()))

pairwise_fert_fung_invsimp_lux=data.frame(y_bot=c(22.5,18.5,0,0),
                                       x_min=as.Date(c("6/20/2018","9/29/2018",NA,NA), format="%m/%d/%Y"),
                                       x_max=as.Date(c("6/30/2018","10/8/2018",NA,NA), format="%m/%d/%Y"),
                                       annot_text=c("*","#","",""),
                                       Root_soil=c(rep("Root",2),rep("Soil",2)))

(lux_fung_invsimp_p=ggplot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=InvSimpson))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus),color=FertStatus),position = position_dodge(8))+
    geom_signif(data=pairwise_fert_fung_invsimp_lux,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    facet_wrap(~Root_soil,nrow = 1,scales = "free_x")+xlab("Inverse Simpson")+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(strip.text = element_blank(),strip.background = element_blank()))

#1700x1200

plot_grid(lux_bact_invsimp_p,lux_fung_invsimp_p,ncol = 1,labels = c('a)', 'b)'), label_size = 30)

#Stats

#Roots
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_invsimp_LUX_bact_mod=lmer(sqrt(InvSimpson)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_root)
plot(root_invsimp_LUX_bact_mod)
hist(resid(root_invsimp_LUX_bact_mod))
qqPlot(resid(root_invsimp_LUX_bact_mod))
shapiro.test(resid(root_invsimp_LUX_bact_mod))
#W = 0.99026, p-value = 0.43



anova(root_invsimp_LUX_bact_mod)
#collectionDate            84.356 16.8711     5 127.01  9.1702 1.802e-07 ***
#FertStatus                 4.706  4.7060     1 127.01  2.5579    0.1122    
#collectionDate:FertStatus  4.996  0.9993     5 127.01  0.5432    0.7433   

emmeans(root_invsimp_LUX_bact_mod,pairwise~FertStatus|collectionDate)





#Soils
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_invsimp_LUX_bact_mod=lmer((InvSimpson)^2~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_invsimp_LUX_bact_mod)
hist(resid(soil_invsimp_LUX_bact_mod))
qqPlot(resid(soil_invsimp_LUX_bact_mod))
shapiro.test(resid(soil_invsimp_LUX_bact_mod))
#W = 0.99024, p-value = 0.0191


anova(soil_invsimp_LUX_bact_mod)
#collectionDate            8.0301e+10 5735804522    14 320.01  5.1549 9.056e-09 ***
#FertStatus                4.3726e+09 4372598865     1 320.00  3.9298   0.04829 *   
#collectionDate:FertStatus 1.3600e+10  971425975    14 320.00  0.8730   0.58878

emmeans(soil_invsimp_LUX_bact_mod,pairwise~FertStatus|collectionDate)

#Roots
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Root")

root_invsimp_LUX_mod=lmer(log(InvSimpson)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_root)
plot(root_invsimp_LUX_mod)
hist(resid(root_invsimp_LUX_mod))
qqPlot(resid(root_invsimp_LUX_mod))
shapiro.test(resid(root_invsimp_LUX_mod))
#W = 0.98585, p-value = 0.1569



anova(root_invsimp_LUX_mod)
#collectionDate            6.8094 1.36188     5 126.04  7.9468 1.548e-06 ***
#FertStatus                1.8709 1.87090     1 126.06 10.9170  0.001241 ** 
#collectionDate:FertStatus 0.5566 0.11131     5 126.06  0.6495  0.662365   

emmeans(root_invsimp_LUX_mod,pairwise~FertStatus|collectionDate)
#$contrasts
#collectionDate = 6/25/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.4082 0.169 126 -2.416  0.0171 

#collectionDate = 10/3/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.3184 0.169 126 -1.884  0.0619 


#Soils
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div,Root_soil=="Soil")

soil_invsimp_LUX_mod=lmer((InvSimpson)~collectionDate*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_div_soil)
plot(soil_invsimp_LUX_mod)
hist(resid(soil_invsimp_LUX_mod))
qqPlot(resid(soil_invsimp_LUX_mod))
shapiro.test(resid(soil_invsimp_LUX_mod))
#W = 0.99593, p-value = 0.5087


anova(soil_invsimp_LUX_mod)
#collectionDate            2914.02 208.144    14 316.02  1.6894 0.05651 .
#FertStatus                  82.39  82.393     1 316.02  0.6688 0.41410  
#collectionDate:FertStatus 1674.90 119.636    14 316.02  0.9710 0.48289  






#Let's graph the pairwise distance between communities

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_mat=read.csv(here::here("R_files","distance_files","bact_GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis.csv"),row.names = 1)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis=as.dist(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M <- matrixConvert(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis, 
                                                      colname = c("sample1", "sample2", "bray"))#total distance 
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M,
                                             GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate",
                                                                                   "UTM_Lat_Cord", "UTM_Lon_Cord")],
                                             by.x = "sample1",by.y = "row.names")
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2)
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2)[4:7]=
  c("s1_Root_soil","s1_collectionDate","s1_UTM_Lat_Cord", "s1_UTM_Lon_Cord")

nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2)
#122265
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta=merge(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_2,
                                                GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate",
                                                                                      "UTM_Lat_Cord", "UTM_Lon_Cord")],
                                                by.x = "sample2",by.y = "row.names")

colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)
colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)[8:11]=
  c("s2_Root_soil","s2_collectionDate","s2_UTM_Lat_Cord", "s2_UTM_Lon_Cord")
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)
#122265
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$day_dis=
  as.numeric(abs(as.Date(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s1_collectionDate, format="%m/%d/%Y")-
                   as.Date(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s2_collectionDate, format="%m/%d/%Y")))


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$spatial_distance=
  sqrt((GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s1_UTM_Lat_Cord-GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s2_UTM_Lat_Cord)^2+
         (GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s1_UTM_Lon_Cord-GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta$s2_UTM_Lon_Cord)^2)

summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta)

#Let's limit it to only within compartment comparisons

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta,s1_Root_soil==s2_Root_soil)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub)
#72139
unique(with(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,interaction(s1_Root_soil,s2_Root_soil)))


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis=as.dist(read.csv(here::here("R_files","distance_files","fungi_GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis.csv"), row.names = 1,header = T))
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M <- matrixConvert(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis, 
                                                           colname = c("sample1", "sample2", "bray"))#total distance 
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M)
#119805
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M,
                                                  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate",
                                                                                             "UTM_Lat_Cord", "UTM_Lon_Cord")],
                                                  by.x = "sample1",by.y = "row.names")
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2)
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2)[4:7]=
  c("s1_Root_soil","s1_collectionDate","s1_UTM_Lat_Cord", "s1_UTM_Lon_Cord")

nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2)
#119805
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta=merge(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_2,
                                                     GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_map[,c("Root_soil","collectionDate",
                                                                                                "UTM_Lat_Cord", "UTM_Lon_Cord")],
                                                     by.x = "sample2",by.y = "row.names")

colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)
colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)[8:11]=
  c("s2_Root_soil","s2_collectionDate","s2_UTM_Lat_Cord", "s2_UTM_Lon_Cord")
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)
#119805
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$day_dis=
  as.numeric(abs(as.Date(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s1_collectionDate, format="%m/%d/%Y")-
                   as.Date(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s2_collectionDate, format="%m/%d/%Y")))


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$spatial_distance=
  sqrt((GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s1_UTM_Lat_Cord-GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s2_UTM_Lat_Cord)^2+
         (GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s1_UTM_Lon_Cord-GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta$s2_UTM_Lon_Cord)^2)

summary(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta)

#Let's limit it to only within compartment comparisons

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta,s1_Root_soil==s2_Root_soil)
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub)
#70596
unique(with(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,interaction(s1_Root_soil,s2_Root_soil)))




#####Change in time####

(chng_time_bact_p=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,aes(x=day_dis,y=bray, color=s1_Root_soil))+
  geom_point(alpha=0.05,position=position_dodge(1))+geom_smooth(method = "lm",size=2)+
  scale_x_continuous(name = "Days")+scale_y_continuous(name = "Pairwis community distance\n(Bray-Curtis)",limits = c(0.12,1))+
  theme_bw(base_size=24)+scale_color_discrete(name=NULL)+theme(legend.position = "none"))

(chng_time_fung_p=ggplot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,aes(x=day_dis,y=bray, color=s1_Root_soil))+
  geom_point(alpha=0.05,position=position_dodge(1))+geom_smooth(method = "lm",size=2)+
  scale_x_continuous(name = "Days")+scale_y_continuous(name = "",limits = c(0.12,1))+
  theme_bw(base_size=24)+scale_color_discrete(name=NULL)+
    theme(axis.text.y = element_blank()))

plot_grid(chng_time_bact_p,chng_time_fung_p,ncol = 2,labels = c('a)', 'b)'), label_size = 30, align = "v",axis = "b",rel_widths = c(1,1.1))
#1500*700



#bacterial

chng_time_B_mod=lm(bray~day_dis*s1_Root_soil,data=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub)
hist(resid(chng_time_B_mod))
qqPlot(resid(chng_time_B_mod))
shapiro.test(resid(chng_time_B_mod))
#Error in shapiro.test(resid(chng_time_B_mod)) : 
#sample size must be between 3 and 5000

summary(chng_time_B_mod)

anova(chng_time_B_mod)

#Split by organ

R_chng_time_B_mod=lm(bray~day_dis,data=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Root"))
hist(resid(R_chng_time_B_mod))
qqPlot(resid(R_chng_time_B_mod))
shapiro.test(resid(R_chng_time_B_mod))
#Error in shapiro.test(resid(chng_time_B_mod)) : 
#sample size must be between 3 and 5000

summary(R_chng_time_B_mod)

anova(R_chng_time_B_mod)


S_chng_time_B_mod=lm(bray~day_dis,data=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Soil"))
hist(resid(S_chng_time_B_mod))
qqPlot(resid(S_chng_time_B_mod))
shapiro.test(resid(S_chng_time_B_mod))
#Error in shapiro.test(resid(chng_time_B_mod)) : 
#sample size must be between 3 and 5000

summary(S_chng_time_B_mod)

anova(S_chng_time_B_mod)

chng_time_F_mod=lm(bray~day_dis*s1_Root_soil,data=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub)
hist(resid(chng_time_F_mod))
qqPlot(resid(chng_time_F_mod))
shapiro.test(resid(chng_time_F_mod))
#Error in shapiro.test(resid(chng_time_F_mod)) : 
#sample size must be between 3 and 5000

summary(chng_time_F_mod)

anova(chng_time_F_mod)

#Split by organ

R_chng_time_F_mod=lm(bray~day_dis,data=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Root"))
hist(resid(R_chng_time_F_mod))
qqPlot(resid(R_chng_time_F_mod))
shapiro.test(resid(R_chng_time_F_mod))
#Error in shapiro.test(resid(chng_time_F_mod)) : 
#sample size must be between 3 and 5000

summary(R_chng_time_F_mod)

anova(R_chng_time_F_mod)


S_chng_time_F_mod=lm(bray~day_dis,data=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Soil"))
hist(resid(S_chng_time_F_mod))
qqPlot(resid(S_chng_time_F_mod))
shapiro.test(resid(S_chng_time_F_mod))
#Error in shapiro.test(resid(chng_time_F_mod)) : 
#sample size must be between 3 and 5000

summary(S_chng_time_F_mod)

anova(S_chng_time_F_mod)



#####Change in distance####

(chng_dist_bact_p=ggplot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,aes(x=spatial_distance,y=bray, color=s1_Root_soil))+
  geom_point(alpha=0.05)+geom_smooth(method = "lm",size=2)+
  scale_x_continuous(name = "Distance (m)")+scale_y_continuous(name = "Pairwis community distance\n(Bray-Curtis)",limits = c(0.12,1))+
  theme_bw(base_size=24)+scale_color_discrete(name=NULL)+theme(legend.position = "none"))

(chng_dist_fung_p=ggplot(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,aes(x=spatial_distance,y=bray, color=s1_Root_soil))+
  geom_point(alpha=0.05)+geom_smooth(method = "lm",size=2)+
  scale_x_continuous(name = "Distance (m)")+scale_y_continuous(name = "",limits = c(0.12,1))+
  theme_bw(base_size=24)+scale_color_discrete(name=NULL)+
    theme(axis.text.y = element_blank()))

plot_grid(chng_dist_bact_p,chng_dist_fung_p,ncol = 2,labels = c('a)', 'b)'), label_size = 30, align = "v",axis = "b",rel_widths = c(1,1.1))

#1500*700

#Bacteria 

chng_dist_F_mod=lm(bray~spatial_distance*s1_Root_soil,data=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub)
hist(resid(chng_dist_F_mod))
qqPlot(resid(chng_dist_F_mod))
shapiro.test(resid(chng_dist_F_mod))
#Error in shapiro.test(resid(chng_dist_F_mod)) : 
#sample size must be between 3 and 5000

summary(chng_dist_F_mod)

anova(chng_dist_F_mod)

#Split by organ

R_chng_dist_B_mod=lm(bray~spatial_distance,data=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Root"))
hist(resid(R_chng_dist_B_mod))
qqPlot(resid(R_chng_dist_B_mod))
shapiro.test(resid(R_chng_dist_B_mod))
#Error in shapiro.test(resid(chng_dist_F_mod)) : 
#sample size must be between 3 and 5000

summary(R_chng_dist_B_mod)

anova(R_chng_dist_B_mod)


S_chng_dist_B_mod=lm(bray~spatial_distance,data=subset(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Soil"))
hist(resid(S_chng_dist_B_mod))
qqPlot(resid(S_chng_dist_B_mod))
shapiro.test(resid(S_chng_dist_B_mod))
#Error in shapiro.test(resid(chng_dist_F_mod)) : 
#sample size must be between 3 and 5000

summary(S_chng_dist_B_mod)

anova(S_chng_dist_B_mod)


#Fungi 


chng_dist_F_mod=lm(bray~spatial_distance*s1_Root_soil,data=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub)
hist(resid(chng_dist_F_mod))
qqPlot(resid(chng_dist_F_mod))
shapiro.test(resid(chng_dist_F_mod))
#Error in shapiro.test(resid(chng_dist_F_mod)) : 
#sample size must be between 3 and 5000

summary(chng_dist_F_mod)

anova(chng_dist_F_mod)

#Split by organ

R_chng_dist_F_mod=lm(bray~spatial_distance,data=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Root"))
hist(resid(R_chng_dist_F_mod))
qqPlot(resid(R_chng_dist_F_mod))
shapiro.test(resid(R_chng_dist_F_mod))
#Error in shapiro.test(resid(chng_dist_F_mod)) : 
#sample size must be between 3 and 5000

summary(R_chng_dist_F_mod)

anova(R_chng_dist_F_mod)


S_chng_dist_F_mod=lm(bray~spatial_distance,data=subset(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_dis_M_meta_sub,s1_Root_soil=="Soil"))
hist(resid(S_chng_dist_F_mod))
qqPlot(resid(S_chng_dist_F_mod))
shapiro.test(resid(S_chng_dist_F_mod))
#Error in shapiro.test(resid(chng_dist_F_mod)) : 
#sample size must be between 3 and 5000

summary(S_chng_dist_F_mod)

anova(S_chng_dist_F_mod)



#####GDM of Space and Time####

#save(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar, file = here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
nsamples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)
#1109
head(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar))


#Let's focus on G5 


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,plotType=="G5"&siteID=="LUX")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
#496
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5) > 0, GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
ntaxa(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5)
#23991

#Bacteria
#Roots 




GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,Root_soil=="Root")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root)
#143
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root) > 0, GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root)
ntaxa(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root)
#13289


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root)[,c("SampleID","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)
#143
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)[colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_metadata=data.frame(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root,method = "bray")
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF=data.frame(as.matrix(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray))
row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site=data.frame("site"=row.names(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF),GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site[1:10,1:10]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site)
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm<- formatsitepair(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_bray_DF_site, 3, siteColumn="site", XColumn="UTM_Lon_Cord", YColumn="UTM_Lat_Cord",
                                                       predData=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_Space_date_metadata)


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod<- gdm(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_TabForm, geo=TRUE)
summary(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod)
plot(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod, plot.layout=c(2,2))
#Percent Deviance Explained:  12.9080179972579

#Space explained 
(sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients[1:3])/sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients))*GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod$explained

#Time explained 
(sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients[4:6])/sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients))*GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod$explained



GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline=isplineExtract(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod)

head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline$x)
root_bact_splin_x=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline$x
summary(root_bact_splin_x)
colnames(root_bact_splin_x)
root_bact_splin_x_M=melt(root_bact_splin_x)
root_bact_splin_x_M[,1]=NULL
colnames(root_bact_splin_x_M)=c("Factor","x")
root_bact_splin_x_M[,"Root_soil"]="Root"

head(root_bact_splin_x_M)
summary(root_bact_splin_x_M)
root_bact_splin_y=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline$y
summary(root_bact_splin_y)
colnames(root_bact_splin_y)

root_bact_splin_y_M=melt(root_bact_splin_y)
root_bact_splin_y_M[,1]=NULL
colnames(root_bact_splin_y_M)=c("Factor","y")
head(root_bact_splin_y_M)
summary(root_bact_splin_y_M)
root_bact_comb_splin=merge(root_bact_splin_x_M,root_bact_splin_y_M, by="row.names")
colnames(root_bact_comb_splin)[2]="Factor"
root_bact_comb_splin$Factor.y=NULL
root_bact_comb_splin$Row.names=NULL
head(root_bact_comb_splin)
colnames(root_bact_comb_splin)


#Soils


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil=subset_samples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5,Root_soil=="Soil")
nsamples(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil)
#353
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil) > 0, GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil)
ntaxa(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil)
#22687


GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data=sample_data(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil)[,c("SampleID","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)
#353
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)[colnames(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)=="SampleID"]="site"
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data)

GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_metadata=data.frame(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil,method = "bray")
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
#Percent Deviance Explained:  15.5395701889918

#Space explained 
(sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients[1:3])/sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients))*GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod$explained

#Time explained 
(sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients[4:6])/sum(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients))*GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod$explained



GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline=isplineExtract(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod)

head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline$x)
soil_bact_splin_x=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline$x
summary(soil_bact_splin_x)
colnames(soil_bact_splin_x)
soil_bact_splin_x_M=melt(soil_bact_splin_x)
soil_bact_splin_x_M[,1]=NULL
colnames(soil_bact_splin_x_M)=c("Factor","x")
soil_bact_splin_x_M[,"Root_soil"]="Soil"

head(soil_bact_splin_x_M)
summary(soil_bact_splin_x_M)
soil_bact_splin_y=GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline$y
summary(soil_bact_splin_y)
colnames(soil_bact_splin_y)

soil_bact_splin_y_M=melt(soil_bact_splin_y)
soil_bact_splin_y_M[,1]=NULL
colnames(soil_bact_splin_y_M)=c("Factor","y")
head(soil_bact_splin_y_M)
summary(soil_bact_splin_y_M)
soil_bact_comb_splin=merge(soil_bact_splin_x_M,soil_bact_splin_y_M, by="row.names")
colnames(soil_bact_comb_splin)[2]="Factor"
soil_bact_comb_splin$Factor.y=NULL
soil_bact_comb_splin$Row.names=NULL
head(soil_bact_comb_splin)
colnames(soil_bact_comb_splin)

root_soil_bact_comb_splin=rbind(root_bact_comb_splin,soil_bact_comb_splin)



(root_soil_bact_GDM_p <- ggplot(root_soil_bact_comb_splin,aes(x=x,y=y,color=Root_soil)) +
    geom_smooth(se=FALSE,method = 'loess',size=2, aes(linetype=Root_soil))+scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)")+
    scale_linetype_manual(values=c("solid","longdash"))+
    scale_color_manual(values=c("red","black"))+facet_wrap(~Factor,scales = "free")+
    theme_bw(base_size=24)+theme())

#Factor
#Fungi


#save(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar, file = here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData"))
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)
#1086
head(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar))


#Let's focus on G5 


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,plotType=="G5"&siteID=="LUX")
nsamples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5)
#490



#Roots 




GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root=subset_samples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,Root_soil=="Root")
nsamples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root)
#141
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root) > 0, GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root)
ntaxa(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root)
#2138


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data=sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root)[,c("sampleID_seq","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)
#141
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)[colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_metadata=data.frame(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root,method = "bray")
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
#Percent Deviance Explained:  7.20588721452468

#Space explained 
(sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients[1:3])/sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients))*GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod$explained

#Time explained 
(sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients[4:6])/sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod$coefficients))*GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod$explained



GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline=isplineExtract(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod)

head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline$x)
root_fung_splin_x=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline$x
summary(root_fung_splin_x)
colnames(root_fung_splin_x)
root_fung_splin_x_M=melt(root_fung_splin_x)
root_fung_splin_x_M[,1]=NULL
colnames(root_fung_splin_x_M)=c("Factor","x")
root_fung_splin_x_M[,"Root_soil"]="Root"

head(root_fung_splin_x_M)
summary(root_fung_splin_x_M)
root_fung_splin_y=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_root_gdmTabMod_spline$y
summary(root_fung_splin_y)
colnames(root_fung_splin_y)

root_fung_splin_y_M=melt(root_fung_splin_y)
root_fung_splin_y_M[,1]=NULL
colnames(root_fung_splin_y_M)=c("Factor","y")
head(root_fung_splin_y_M)
summary(root_fung_splin_y_M)
root_fung_comb_splin=merge(root_fung_splin_x_M,root_fung_splin_y_M, by="row.names")
colnames(root_fung_comb_splin)[2]="Factor"
root_fung_comb_splin$Factor.y=NULL
root_fung_comb_splin$Row.names=NULL
head(root_fung_comb_splin)
colnames(root_fung_comb_splin)



#Soils




GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil=subset_samples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5,Root_soil=="Soil")
nsamples(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil)
#349
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil) > 0, GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil)
ntaxa(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil)
#2539


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data=sample_data(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil)[,c("sampleID_seq","collectionDate","UTM_Lat_Cord","UTM_Lon_Cord")]
nrow(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)
#349
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate_N=as.numeric(as.Date(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data$collectionDate,format = "%m/%d/%Y"))

colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)[colnames(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)=="sampleID_seq"]="site"
head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data)

GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_metadata=data.frame(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_Space_date_data[,c("site","UTM_Lat_Cord","UTM_Lon_Cord","collectionDate_N")])






GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_bray=phyloseq::distance(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil,method = "bray")
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
#Percent Deviance Explained:  15.0598568158892

#Space explained 
(sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients[1:3])/sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients))*GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod$explained

#Time explained 
(sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients[4:6])/sum(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod$coefficients))*GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod$explained


GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline=isplineExtract(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod)

head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline$x)
soil_fung_splin_x=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline$x
summary(soil_fung_splin_x)
colnames(soil_fung_splin_x)
soil_fung_splin_x_M=melt(soil_fung_splin_x)
soil_fung_splin_x_M[,1]=NULL
colnames(soil_fung_splin_x_M)=c("Factor","x")
soil_fung_splin_x_M[,"Root_soil"]="Soil"

head(soil_fung_splin_x_M)
summary(soil_fung_splin_x_M)
soil_fung_splin_y=GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5_soil_gdmTabMod_spline$y
summary(soil_fung_splin_y)
colnames(soil_fung_splin_y)

soil_fung_splin_y_M=melt(soil_fung_splin_y)
soil_fung_splin_y_M[,1]=NULL
colnames(soil_fung_splin_y_M)=c("Factor","y")
head(soil_fung_splin_y_M)
summary(soil_fung_splin_y_M)
soil_fung_comb_splin=merge(soil_fung_splin_x_M,soil_fung_splin_y_M, by="row.names")
colnames(soil_fung_comb_splin)[2]="Factor"
soil_fung_comb_splin$Factor.y=NULL
soil_fung_comb_splin$Row.names=NULL
head(soil_fung_comb_splin)
colnames(soil_fung_comb_splin)

root_soil_fung_comb_splin=rbind(root_fung_comb_splin,soil_fung_comb_splin)



(root_soil_fung_GDM_p <- ggplot(root_soil_fung_comb_splin,aes(x=x,y=y,color=Root_soil)) +
    geom_smooth(se=FALSE,method = 'loess',size=2, aes(linetype=Root_soil))+scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)")+
    scale_linetype_manual(values=c("solid","longdash"))+
    scale_color_manual(values=c("red","black"))+facet_wrap(~Factor,scales = "free")+
    theme_bw(base_size=24)+theme())


#Change over time 

(Time_root_soil_bact_GDM_p <- ggplot(subset(root_soil_bact_comb_splin, Factor=="collectionDate_N"),aes(x=as.Date(x,origin="1970-01-01"),y=y,color=Root_soil)) +
    geom_smooth(se=FALSE,method = 'loess',size=2, aes(linetype=Root_soil))+scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)",limits = c(-0.001,0.25))+
    scale_linetype_manual(values=c("solid","longdash"),name=NULL)+
    scale_color_discrete(name=NULL)+scale_x_date(name= NULL)+
    theme_bw(base_size=24)+theme(legend.position = "none"))

(Time_root_soil_fung_GDM_p <- ggplot(subset(root_soil_fung_comb_splin, Factor=="collectionDate_N"),aes(x=as.Date(x,origin="1970-01-01"),y=y,color=Root_soil)) +
    geom_smooth(se=FALSE,method = 'loess',size=2, aes(linetype=Root_soil))+scale_y_continuous(name = "",limits = c(-0.001,0.25))+
    scale_linetype_manual(values=c("solid","longdash"),name=NULL)+
    scale_color_discrete(name=NULL)+scale_x_date(name= NULL)+
    theme_bw(base_size=24)+theme(axis.text.y = element_blank()))


plot_grid(Time_root_soil_bact_GDM_p,Time_root_soil_fung_GDM_p,ncol = 2,labels = c('a)', 'b)'), label_size = 30, align = "v",axis = "b",rel_widths = c(1,1.1))


#Change over time 

(Dist_root_soil_bact_GDM_p <- ggplot(subset(root_soil_bact_comb_splin, Factor=="Geographic"),aes(x=x,y=y,color=Root_soil)) +
    geom_smooth(se=FALSE,method = 'loess',size=2, aes(linetype=Root_soil))+scale_y_continuous(name = "Partial ecological distance \n (Bray-Curtis)",limits = c(-0.001,0.34))+
    scale_linetype_manual(values=c("solid","longdash"),name=NULL)+
    scale_color_discrete(name=NULL)+scale_x_continuous(name= "Distance (m)")+
    theme_bw(base_size=24)+theme(legend.position = "none"))

(Dist_root_soil_fung_GDM_p <- ggplot(subset(root_soil_fung_comb_splin, Factor=="Geographic"),aes(x=x,y=y,color=Root_soil)) +
    geom_smooth(se=FALSE,method = 'loess',size=2, aes(linetype=Root_soil))+scale_y_continuous(name = "",limits = c(-0.001,0.34))+
    scale_linetype_manual(values=c("solid","longdash"),name=NULL)+
    scale_color_discrete(name=NULL)+scale_x_continuous(name= "Distance (m)")+
    theme_bw(base_size=24)+theme())

#

plot_grid(Dist_root_soil_bact_GDM_p,Dist_root_soil_fung_GDM_p,ncol = 2,labels = c('a)', 'b)'), label_size = 30, align = "v",axis = "b",rel_widths = c(1,1.1))
