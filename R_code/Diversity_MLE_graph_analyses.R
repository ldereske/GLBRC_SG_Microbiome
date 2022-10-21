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

#####Calculate Diversity metrics####
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_map=sample_data(GLBRC018_OTU_bact_MMPRNT_All_sites_G5)
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_map=sample_data(GLBRC018_OTU_fung_MMPRNT_All_sites_G5)

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
site_labels=c("LUX"="Lux Arbor","LC"="Lake City","ESC"="Escanaba","HAN"= "Hancock","RHN"="Rhinelander")


####Richness####

#Stats 

#Roots
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_rich_all_site_bact_mod=lmer((Observed)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root)
plot(root_rich_all_site_bact_mod)
hist(resid(root_rich_all_site_bact_mod))
qqPlot(resid(root_rich_all_site_bact_mod))
shapiro.test(resid(root_rich_all_site_bact_mod))
#W = 0.98842, p-value = 0.4491
simulateResiduals(fittedModel = root_rich_all_site_bact_mod, plot = T)

anova(root_rich_all_site_bact_mod)
#siteID            2837112  709278     4 92.997 39.6794 < 2.2e-16 ***
#FertStatus         291606  291606     1 91.919 16.3134 0.0001112 ***
#siteID:FertStatus  209668   52417     4 91.921  2.9324 0.0248494 * 

emmeans(root_rich_all_site_bact_mod,pairwise~FertStatus|siteID)

#$contrasts
#siteID = LC:
#contrast      estimate   SE   df t.ratio p.value
#Fert - Unfert  -192.33 54.6 92.0  -3.524  0.0007
#siteID = RHN:
#contrast      estimate   SE   df t.ratio p.value
#Fert - Unfert  -202.67 54.6 92.0  -3.713  0.0004

emmeans(root_rich_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN    209.4 42.2 95.2   4.963  <.0001
#ESC - LC     100.0 38.6 92.0   2.590  0.0806
#ESC - LUX    -46.3 39.1 92.3  -1.186  0.7594
#ESC - RHN    387.0 38.6 92.0  10.028  <.0001
#HAN - LC    -109.5 42.2 95.2  -2.594  0.0795
#HAN - LUX   -255.8 42.5 94.5  -6.019  <.0001
#HAN - RHN    177.6 42.2 95.2   4.209  0.0005
#LC - LUX    -146.3 39.1 92.3  -3.744  0.0028
#LC - RHN     287.1 38.6 92.0   7.438  <.0001
#LUX - RHN    433.4 39.1 92.3  11.090  <.0001


#Soil
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_rich_all_site_bact_mod=lmer(Observed~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil)
plot(soil_rich_all_site_bact_mod)
hist(resid(soil_rich_all_site_bact_mod))
qqPlot(resid(soil_rich_all_site_bact_mod))
shapiro.test(resid(soil_rich_all_site_bact_mod))
#W = 0.99289, p-value = 0.8277
simulateResiduals(fittedModel = soil_rich_all_site_bact_mod, plot = T)

anova(soil_rich_all_site_bact_mod)
#siteID            1877527  469382     4 93.992 32.6514 <2e-16 ***



emmeans(soil_rich_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN     57.1 37.9 95.2   1.505  0.5618
#ESC - LC     -31.8 34.6 93.0  -0.917  0.8896
#ESC - LUX    -64.2 34.6 93.0  -1.855  0.3488
#ESC - RHN    286.6 34.6 93.0   8.280  <.0001
#HAN - LC     -88.8 37.9 95.2  -2.343  0.1407
#HAN - LUX   -121.3 37.9 95.2  -3.199  0.0157
#HAN - RHN    229.5 37.9 95.2   6.055  <.0001
#LC - LUX     -32.5 34.6 93.0  -0.938  0.8814
#LC - RHN     318.3 34.6 93.0   9.197  <.0001
#LUX - RHN    350.8 34.6 93.0  10.135  <.0001




#Roots
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_rich_all_site_fung_mod=lmer(sqrt(Observed)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root)
plot(root_rich_all_site_fung_mod)
hist(resid(root_rich_all_site_fung_mod))
qqPlot(resid(root_rich_all_site_fung_mod))
shapiro.test(resid(root_rich_all_site_fung_mod))
#W = 0.98082, p-value = 0.1111
simulateResiduals(fittedModel = root_rich_all_site_fung_mod, plot = T)

anova(root_rich_all_site_fung_mod)
#siteID            61.781 15.4453     4   101  9.3761 1.741e-06 *** 

emmeans(root_rich_all_site_fung_mod,pairwise~FertStatus|siteID)

#$contrasts
#siteID = HAN:
#contrast      estimate    SE   df t.ratio p.value
#Fert - Unfert  -1.3057 0.605 90.1  -2.158  0.0336



emmeans(root_rich_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN    1.404 0.406 94.1   3.457  0.0072
#ESC - LC     0.291 0.384 91.4   0.760  0.9414
#ESC - LUX   -0.975 0.375 90.6  -2.601  0.0786
#ESC - RHN    0.450 0.375 90.6   1.199  0.7517
#HAN - LC    -1.113 0.411 96.0  -2.706  0.0604
#HAN - LUX   -2.380 0.403 95.0  -5.910  <.0001
#HAN - RHN   -0.955 0.403 95.0  -2.370  0.1326
#LC - LUX    -1.267 0.379 90.9  -3.340  0.0104
#LC - RHN     0.158 0.379 90.9   0.417  0.9935
#LUX - RHN    1.425 0.371 90.1   3.847  0.0020





#Soil
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_rich_all_site_fung_mod=lmer(sqrt(Observed)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil)
plot(soil_rich_all_site_fung_mod)
hist(resid(soil_rich_all_site_fung_mod))
qqPlot(resid(soil_rich_all_site_fung_mod))
shapiro.test(resid(soil_rich_all_site_fung_mod))
#W = 0.99242, p-value = 0.3085
simulateResiduals(fittedModel = soil_rich_all_site_fung_mod, plot = T)

anova(soil_rich_all_site_fung_mod)
#siteID            82.237 20.5594     4   102 28.2327 8.618e-16 ***



emmeans(soil_rich_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate    SE   df t.ratio p.value
#ESC - HAN   0.4673 0.268 96.0   1.746  0.4114
#ESC - LC    0.0161 0.246 91.1   0.065  1.0000
#ESC - LUX  -1.9435 0.249 91.6  -7.794  <.0001
#ESC - RHN   0.1712 0.249 91.6   0.687  0.9589
#HAN - LC   -0.4513 0.268 96.0  -1.686  0.4476
#HAN - LUX  -2.4109 0.271 96.5  -8.908  <.0001
#HAN - RHN  -0.2962 0.270 95.1  -1.097  0.8079
#LC - LUX   -1.9596 0.249 91.6  -7.859  <.0001
#LC - RHN    0.1551 0.249 91.6   0.622  0.9712
#LUX - RHN   2.1147 0.252 92.0   8.378  <.0001



###InvSimpson####


#Roots
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_invsimp_all_site_bact_mod=lmer(log(InvSimpson)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_root)
plot(root_invsimp_all_site_bact_mod)
hist(resid(root_invsimp_all_site_bact_mod))
qqPlot(resid(root_invsimp_all_site_bact_mod))
shapiro.test(resid(root_invsimp_all_site_bact_mod))
#W = 0.98672, p-value = 0.3322



anova(root_invsimp_all_site_bact_mod)
#siteID            9.6757 2.41893     4   103  15.247 7.99e-10 ***
#FertStatus        1.7852 1.78523     1   103  11.253 0.001114 ** 
#siteID:FertStatus 0.8643 0.21608     4   103   1.362 0.252315    



emmeans(root_invsimp_all_site_bact_mod,pairwise~FertStatus|siteID)
#siteID = LC:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.4439 0.163 92.1  -2.730  0.0076

#siteID = RHN:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.4283 0.163 92.1  -2.634  0.0099

#$contrasts
#siteID = ESC:
# contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert  -0.3085 0.163 92.1  -1.897  0.0609

emmeans(root_invsimp_all_site_bact_mod,pairwise~siteID)
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN   0.2326 0.125 97.0   1.862  0.3448
#ESC - LC   -0.0906 0.115 92.1  -0.788  0.9335
#ESC - LUX  -0.1868 0.116 92.6  -1.605  0.4978
#ESC - RHN   0.6152 0.115 92.1   5.351  <.0001
#HAN - LC   -0.3233 0.125 97.0  -2.587  0.0807
#HAN - LUX  -0.4195 0.126 96.1  -3.328  0.0106
#HAN - RHN   0.3826 0.125 97.0   3.062  0.0233
#LC - LUX   -0.0962 0.116 92.6  -0.827  0.9217
#LC - RHN    0.7058 0.115 92.1   6.139  <.0001
#LUX - RHN   0.8021 0.116 92.6   6.892  <.0001



#Soil
GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_invsimp_all_site_bact_mod=lmer(InvSimpson~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div_soil)
plot(soil_invsimp_all_site_bact_mod)
hist(resid(soil_invsimp_all_site_bact_mod))
qqPlot(resid(soil_invsimp_all_site_bact_mod))
shapiro.test(resid(soil_invsimp_all_site_bact_mod))
#W = 0.98578, p-value = 0.273




anova(soil_invsimp_all_site_bact_mod)
#siteID            134024   33506     4 93.833  7.3166 3.575e-05 ***


emmeans(soil_invsimp_all_site_bact_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN     21.3 21.3 97.3   1.003  0.8534
#ESC - LC      10.5 19.5 93.1   0.536  0.9834
#ESC - LUX    -36.2 19.5 93.1  -1.853  0.3500
#ESC - RHN     67.2 19.5 93.1   3.441  0.0076
#HAN - LC     -10.9 21.3 97.3  -0.511  0.9861
#HAN - LUX    -57.5 21.3 97.3  -2.704  0.0607
#HAN - RHN     45.9 21.3 97.3   2.156  0.2055
#LC - LUX     -46.7 19.5 93.1  -2.389  0.1276
#LC - RHN      56.8 19.5 93.1   2.906  0.0361
#LUX - RHN    103.4 19.5 93.1   5.294  <.0001


#Roots
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Root")

root_invsimp_all_site_fung_mod=lmer(log(InvSimpson)~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_root)
plot(root_invsimp_all_site_fung_mod)
hist(resid(root_invsimp_all_site_fung_mod))
qqPlot(resid(root_invsimp_all_site_fung_mod))
shapiro.test(resid(root_invsimp_all_site_fung_mod))
#W = 0.98727, p-value = 0.3799


anova(root_invsimp_all_site_fung_mod)
#siteID            2.2187 0.55467     4 90.896  2.9090 0.02581 * 
#FertStatus        0.1326 0.13265     1 89.554  0.6957 0.40646   
#siteID:FertStatus 3.2841 0.82103     4 89.570  4.3060 0.00312 **

emmeans(root_invsimp_all_site_fung_mod,pairwise~FertStatus|siteID)


#$contrasts
#siteID = ESC:
#contrast      estimate    SE   df t.ratio p.value
#Fert - Unfert   0.6081 0.183 90.8   3.330  0.0013

emmeans(root_invsimp_all_site_fung_mod,pairwise~siteID)
#contrast  estimate    SE    df t.ratio p.value
#ESC - HAN -0.17800 0.139 93.3  -1.285  0.7010
#ESC - LC   0.07510 0.131 91.1   0.575  0.9784
#ESC - LUX -0.00146 0.128 90.4  -0.011  1.0000
#ESC - RHN -0.30599 0.128 90.4  -2.398  0.1253
#HAN - LC   0.25310 0.140 95.0   1.802  0.3786
#HAN - LUX  0.17655 0.137 94.1   1.284  0.7014
#HAN - RHN -0.12798 0.137 94.1  -0.931  0.8841
#LC - LUX  -0.07655 0.129 90.6  -0.593  0.9758
#LC - RHN  -0.38109 0.129 90.6  -2.952  0.0320
#LUX - RHN -0.30453 0.126 90.1  -2.416  0.1205



#Soil
GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil=subset(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,Root_soil=="Soil")

soil_invsimp_all_site_fung_mod=lmer(InvSimpson~siteID*FertStatus+(1|plotRep),data = GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div_soil)
plot(soil_invsimp_all_site_fung_mod)
hist(resid(soil_invsimp_all_site_fung_mod))
qqPlot(resid(soil_invsimp_all_site_fung_mod))
shapiro.test(resid(soil_invsimp_all_site_fung_mod))
#W = 0.99158, p-value = 0.7256


anova(soil_invsimp_all_site_fung_mod)
#siteID            978.07 244.518     4 92.431  3.3455 0.01323 *
#FertStatus          1.90   1.898     1 91.297  0.0260 0.87235  
#siteID:FertStatus 558.02 139.506     4 91.274  1.9087 0.11566  


emmeans(soil_invsimp_all_site_fung_mod,pairwise~siteID)
#$contrasts
#contrast  estimate   SE    df t.ratio p.value
#ESC - HAN    1.927 2.70 94.4   0.715  0.9527
#ESC - LC     8.444 2.47 91.0   3.421  0.0081
#ESC - LUX    2.142 2.50 91.3   0.857  0.9115
#ESC - RHN    4.225 2.50 91.4   1.691  0.4448
#HAN - LC     6.517 2.70 94.4   2.417  0.1198
#HAN - LUX    0.216 2.73 94.8   0.079  1.0000
#HAN - RHN    2.298 2.72 93.7   0.846  0.9154
#LC - LUX    -6.301 2.50 91.3  -2.522  0.0946
#LC - RHN    -4.219 2.50 91.4  -1.689  0.4461
#LUX - RHN    2.082 2.53 91.7   0.823  0.9229


#Graphing diversity plots####
site_colors=c("ESC"="#E41A1C", "HAN"="#377EB8", "LC"="#4DAF4A", 
              "LUX"="#984EA3", "RHN"="#FF7F00")
pairwise_fert_bact_rich=data.frame(y_bot=c(1660,1230),
                                   x_min=c(1.8,4.8),
                                   x_max=c(2.2,5.2),
                                   annot_text=c("***"," *** "),
                                   Root_soil=c(rep("Root",2)))

root_bact_rich_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("a","bc","ab","c","d"),
                                  Root_soil=rep("Root"))

soil_bact_rich_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("A","AB","AB","B","C"),
                                  Root_soil=rep("Soil"))

root_soil_bact_rich_letters=rbind(root_bact_rich_letters,soil_bact_rich_letters)
bact_rich_max=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(Observed))

bact_rich_max_disp_letters=merge(root_soil_bact_rich_letters,bact_rich_max,by=c("siteID","Root_soil"))


root_fung_rich_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("a","bc","ac","b","bc"),
                                  Root_soil=rep("Root"))

soil_fung_rich_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                  sig_let=c("A","B","B","B","B"),
                                  Root_soil=rep("Soil"))
root_soil_fung_rich_letters=rbind(root_fung_rich_letters,soil_fung_rich_letters)

fung_rich_max=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(Observed))

fung_rich_max_rich_letters=merge(root_soil_fung_rich_letters,fung_rich_max,by=c("siteID","Root_soil"))



pairwise_fert_bact_invSimp=data.frame(y_bot=c(0,95,100,0,55,0,0,0,0,0),
                                      x_min=c(0,1.85,2.85,0,4.85,0,0,0,0,0),
                                      x_max=c(0,2.15,3.15,0,5.15,0,0,0,0,0),
                                      annot_text=c(NA,"**","#",NA," ** ",NA,NA,NA,NA,NA),
                                      Root_soil=c(Root_soil=c(rep("Root",5),rep("Soil",5))))


root_bact_invSimp_letters=data.frame(siteID=c("LUX","LC","ESC", "HAN","RHN"),
                                     sig_let=c("a","ab","ab","b","c"),
                                     Root_soil=rep("Root"))

soil_bact_invSimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                     sig_let=c("A","AB","A","A","B"),
                                     Root_soil=rep("Soil"))



root_soil_bact_invSimp_letters=rbind(root_bact_invSimp_letters,soil_bact_invSimp_letters)
bact_invSimp_max=GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(InvSimpson))

bact_invSimp_max_disp_letters=merge(root_soil_bact_invSimp_letters,bact_invSimp_max,by=c("siteID","Root_soil"))


pairwise_fert_fung_invsimp=data.frame(y_bot=c(0,0,17,0,0,0,0,0,0,0),
                                      x_min=c(0,0,2.85,0,0,0,0,0,0,0),
                                      x_max=c(0,0,3.15,0,0,0,0,0,0,0),
                                      annot_text=c("","","**","","","","","","",""),
                                      Root_soil=c(Root_soil=c(rep("Root",5),rep("Soil",5))))

root_fung_invsimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                     sig_let=c("ab","ab","a","ab","b"),
                                     Root_soil=rep("Root"))

soil_fung_invsimp_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                                     sig_let=c("A","AB","B","AB","AB"),
                                     Root_soil=rep("Soil"))
root_soil_fung_invsimp_letters=rbind(root_fung_invsimp_letters,soil_fung_invsimp_letters)

fung_invsimp_max=GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div%>%group_by(siteID,Root_soil)%>%summarise(max_value=max(InvSimpson))

fung_invsimp_max_invsimp_letters=merge(root_soil_fung_invsimp_letters,fung_invsimp_max,by=c("siteID","Root_soil"))



(bact_rich_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div)+
    geom_boxplot(aes(y=Observed, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = bact_rich_max_disp_letters, aes(x=siteID, y = 250 + max_value, label = sig_let), vjust=0, size=10)+
    geom_signif(data=pairwise_fert_bact_rich,aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 10,manual = T)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Bacterial\nRichness",limits = c(700,2900))+
    theme(axis.text.x = element_blank(),legend.position = "none",
          strip.background = element_rect(fill = "white",color = "black")))


(fung_rich_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,aes(y=Observed, x=factor(siteID,levels = site_order)))+
    geom_boxplot(aes(y=Observed, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = fung_rich_max_rich_letters, aes(x=siteID, y = 10 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    scale_y_continuous(name = "Fungal\nRichness",limits = c(98,520))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none",axis.text.x = element_blank()))



(bact_invSimp_p3=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div)+
    geom_boxplot(aes(y=InvSimpson, x=factor(siteID,levels = site_order),
                     fill=FertStatus,color=FertStatus))+
    geom_text(data = bact_invSimp_max_disp_letters, aes(x=siteID, y = 40 + max_value, label = sig_let), vjust=0, size=10)+
    geom_signif(data=pairwise_fert_bact_invSimp,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    scale_fill_manual(values = c("black","grey"), name=NULL)+scale_color_manual(values = c("darkgrey","black"), name=NULL)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Bacterial\nInverse Simpson",limits = c(0,520))+
    theme(axis.text.x = element_blank(),legend.key.size = unit(3, 'cm'),legend.text = element_text(size=32)))

(bact_invSimp_p2=ggplot(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_div)+
    geom_boxplot(aes(y=InvSimpson, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = bact_invSimp_max_disp_letters, aes(x=siteID, y = 40 + max_value, label = sig_let), vjust=0, size=10)+
    geom_signif(data=pairwise_fert_bact_invSimp,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    scale_y_continuous(name = "Bacterial\nInverse Simpson",limits = c(0,520))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none",axis.text.x = element_blank()))


(fung_invSimp_p2=ggplot(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_div,aes(y=Observed, x=factor(siteID,levels = site_order)))+
    geom_boxplot(aes(y=InvSimpson, x=factor(siteID,levels = site_order),
                     fill=siteID,alpha=FertStatus),color="black")+
    geom_text(data = fung_invsimp_max_invsimp_letters, aes(x=siteID, y = 1 + max_value, label = sig_let), vjust=0, size=10)+
    facet_wrap(~Root_soil,nrow = 1)+scale_x_discrete(labels=site_labels,name=NULL)+theme_cowplot(font_size = 24)+
    geom_signif(data=pairwise_fert_fung_invsimp,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    scale_fill_manual(values = site_colors, name=NULL)+scale_alpha_manual(values = c(1,0.1))+
    scale_y_continuous(name = "Fungal\nInverse Simpson",limits = c(0,57))+
    theme(strip.text = element_blank(),strip.background = element_blank(),legend.position = "none"))




(rich_invSimp_4panel=plot_grid(bact_rich_p2,fung_rich_p2,bact_invSimp_p2,fung_invSimp_p2,ncol = 1,
                               labels = c('a)', 'b)','c)','d)'), label_size = 30,align = "v"))



plot_grid(rich_invSimp_4panel,get_legend(bact_invSimp_p3),ncol = 2,rel_widths = c(4,1.5))



ggsave(plot_grid(rich_invSimp_4panel,get_legend(bact_invSimp_p3),ncol = 2,rel_widths = c(4,1.5)),
       filename = "richness_invSimpson_boxplot_All_Sites_p.png",path = here::here("Manuscript","MLE_comm_figs"),width = 25,height =20)

ggsave(plot_grid(rich_invSimp_4panel,get_legend(bact_invSimp_p3),ncol = 2,rel_widths = c(4,1.5)),
       filename = "richness_invSimpson_boxplot_All_Sites_p.svg",path = here::here("Manuscript","MLE_comm_figs"),width = 25,height =20)

