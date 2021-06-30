#####Code for the analysis of the Metadata from the GLBRC Scale up experiments####
here::i_am("R_code/Bell-Dereske&Evans_Rain_analyses_figures.R")
library(plyr); library(dplyr)
library(reshape2)
library(tidyr)
library(lme4)
library(car)
library(ggplot2)
library(vegan)
library(otuSummary)
library(stringr)
library(cowplot)
library(lmerTest)
library(emmeans)
library(factoextra)
library(ggsignif)
library("viridis") 



#Outlier filtered metadata
#Metadata not included in this GIT directory. Let me know if you need the raw files
MMPRNT_16s.metadata_mapp_raw=read.table(file = "D:/MMPRNT_16S_016-018/ABIOTIC_data/MMPRNT.metadata_mapp_raw_comp_OUTLIER_filter.txt")

colnames(MMPRNT_16s.metadata_mapp_raw)
unique(MMPRNT_16s.metadata_mapp_raw$year)

MMPRNT_2018.metadata_mapp_raw=subset(MMPRNT_16s.metadata_mapp_raw,year==2018)
nrow(MMPRNT_2018.metadata_mapp_raw)
#870
colSums(is.na(MMPRNT_2018.metadata_mapp_raw))


MMPRNT_2018.metadata_mapp_coverage=MMPRNT_2018.metadata_mapp_raw[,c("SampleID","collectionDate","siteID","plotType","FertStatus","plotRep","xCord",
                                                                    "yCord","abs_yCord","dis_Mid","UTM_Lat_Cord","UTM_Lon_Cord",
                                                                    "percent_soil_moisture_dry_weight","dry_matter_yield_mg_ha_mean",
                                                                    "stem.fresh.weight","stem.displacement","extra.leaves.fresh.weight",
                                                                    "top.leaf.fresh.weight","total.shoot.fresh.weight","stem.dry.weight",
                                                                    "extra.leaves.dry.weight","top.leaf.dry.weight","total.shoot.dry.weight",
                                                                    "stem.bottom.10cm.dry.weight","plant.height","top.leaf.area",
                                                                    "specific.leaf.area","specific.stem.density","shoot.water.percentage",
                                                                    "core_root_mass_subplot","pH_raw_MMPRNT","pH_MMPRNT_subplot",
                                                                    "vwc_avg_mean","vwc_2_avg_mean","soil_temp_1_avg_mean","soil_temp_2_avg_mean",
                                                                    "airtc_avg_mean","rh_mean","par_den_avg_mean","par_tot_tot_mean",
                                                                    "ws_ms_avg_mean","rain_mm_tot_sum","slrmj_tot_mean","vwc_avg_CV",
                                                                    "vwc_2_avg_CV","soil_temp_1_avg_CV","soil_temp_2_avg_CV","airtc_avg_CV",
                                                                    "rh_CV","par_den_avg_CV","par_tot_tot_CV","ws_ms_avg_CV",
                                                                    "slrmj_tot_CV","ugN_fixed_gdrysoil_day","TOC","TON","MBC_NA","MBN_NA",
                                                                    "MBC_MBN_ratio","TOC_TON_ratio","ph","lime_index","p_ppm",
                                                                    "k_ppm","ca_ppm","mg_ppm","event_age_1mm","ugN_NH4_g_dry_soil_na",
                                                                    "ugN_NO3_g_dry_soil_na","potential_NH4_mineralization_rate_na",
                                                                    "potential_NO3_nitrification_rate_na","past_7d_rain","event_age")]

summary(MMPRNT_2018.metadata_mapp_coverage)


#####All Site Abiotic####



MMPRNT_2018.metadata_mapp_coverage_all_sites=subset(MMPRNT_2018.metadata_mapp_coverage,siteID!="LUX"|collectionDate=="7/30/2018")
nrow(MMPRNT_2018.metadata_mapp_coverage_all_sites)
#486


MMPRNT_2018.metadata_mapp_coverage_all_sites=subset(MMPRNT_2018.metadata_mapp_coverage_all_sites,siteID!="LC"|collectionDate=="7/10/2018")
nrow(MMPRNT_2018.metadata_mapp_coverage_all_sites)
#294

colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_all_sites))

unique(MMPRNT_2018.metadata_mapp_coverage_all_sites$siteID)
#"LC"  "HAN" "RHN" "ESC" "LUX"
unique(MMPRNT_2018.metadata_mapp_coverage_all_sites$plotType)

MMPRNT_2018.metadata_mapp_coverage_all_sites_G5=subset(MMPRNT_2018.metadata_mapp_coverage_all_sites,plotType=="G5")
nrow(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5)
#114
colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_all_sites_G5))


MMPRNT_2018.metadata_mapp_all_sites_G5=MMPRNT_2018.metadata_mapp_coverage_all_sites_G5
MMPRNT_2018.metadata_mapp_all_sites_G5[,c("event_age","pH_raw_MMPRNT","pH_MMPRNT_subplot","rain_mm_tot_sum")]=NULL

colnames(MMPRNT_2018.metadata_mapp_all_sites_G5)

(MMPRNT_2018.metadata_mapp_all_sites_G5)


#PCA 
summary(MMPRNT_2018.metadata_mapp_all_sites_G5[,13:69])

#need to remove NA's
unique(subset(MMPRNT_2018.metadata_mapp_all_sites_G5,!is.na(vwc_avg_mean))$siteID)
#LC has no VWC values
row.names(MMPRNT_2018.metadata_mapp_all_sites_G5)=MMPRNT_2018.metadata_mapp_all_sites_G5$SampleID

res.pca <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c(13:69)]), scale = TRUE)

fviz_eig(res.pca)

fviz_pca_ind(res.pca, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_ind(res.pca, habillage=na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c(3,13:69)])$siteID,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)


fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c(3,13:69)])$siteID,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)


#Let's remove VWC so we can include the pionts from Lake City
colnames(MMPRNT_2018.metadata_mapp_all_sites_G5)

res.pca_v2 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c(13:30,33:40,43:69)]), scale = TRUE)

fviz_eig(res.pca_v2)

fviz_pca_ind(res.pca_v2, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_ind(res.pca_v2, habillage=na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c(3,13:30,33:40,43:69)])$siteID,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)


fviz_pca_var(res.pca_v2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca_v2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)

fviz_pca_biplot(res.pca_v2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c(3,13:30,33:40,43:69)])$siteID,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)



#Plant data only


colnames(MMPRNT_2018.metadata_mapp_all_sites_G5)

res.pca_v3 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("dry_matter_yield_mg_ha_mean","total.shoot.dry.weight",
                                                                       "specific.leaf.area","plant.height","core_root_mass_subplot"
                                                                       )]), scale = TRUE)

fviz_eig(res.pca_v3)

fviz_pca_ind(res.pca_v3, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Plant data")


fviz_pca_ind(res.pca_v3, habillage=
               na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("dry_matter_yield_mg_ha_mean","total.shoot.dry.weight","siteID",
                                                                 "specific.leaf.area","plant.height","core_root_mass_subplot")])$siteID,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Plant data")


fviz_pca_var(res.pca_v3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Plant data")


fviz_pca_biplot(res.pca_v3, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)+ggtitle("Plant data")

fviz_pca_biplot(res.pca_v3, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("dry_matter_yield_mg_ha_mean","total.shoot.dry.weight","siteID",
                                                                    "specific.leaf.area","plant.height","core_root_mass_subplot")])$siteID,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Plant data")


value_max_shoot_bio<-MMPRNT_2018.metadata_mapp_all_sites_G5 %>% group_by(siteID) %>% summarize(max_value = max(total.shoot.dry.weight))
value_max_yield<-MMPRNT_2018.metadata_mapp_all_sites_G5 %>% group_by(siteID) %>% summarize(max_value = max(dry_matter_yield_mg_ha_mean))


(shoot_bio_p<-ggplot(MMPRNT_2018.metadata_mapp_all_sites_G5,aes(x=siteID, y=total.shoot.dry.weight))+
  geom_boxplot(data = MMPRNT_2018.metadata_mapp_all_sites_G5,aes(x=siteID, y=total.shoot.dry.weight,color=FertStatus))+
  geom_text(data = value_max_shoot_bio, aes(x=siteID, y = 0.1 + max_value, label = shoot_bio_letters$sig_let), vjust=0, size=10)+
  geom_signif(y_position = c(13, 6.5), xmin = c(3.85, 4.85), xmax = c(4.15, 5.15),annotation = c("**", "**"), 
              tip_length = 0,textsize = 8)+theme_classic(base_size = 25)+theme(legend.position = "none"))

(yield_prod_p<-ggplot(MMPRNT_2018.metadata_mapp_all_sites_G5,aes(x=siteID, y=dry_matter_yield_mg_ha_mean))+
  geom_boxplot(data = MMPRNT_2018.metadata_mapp_all_sites_G5_yield,aes(x=siteID, y=dry_matter_yield_mg_ha_mean,color=FertStatus))+
  geom_text(data = value_max_yield, aes(x=siteID, y = 0.1 + max_value, label = yield_letters$sig_let), vjust=0, size=10)+
  geom_signif(y_position = c(9.3, 6.5), xmin = c(0.85,4.85), xmax = c(1.15,5.15),annotation = c("*"), 
            tip_length = 0,textsize = 8)+theme_classic(base_size = 25))

plot_grid(shoot_bio_p,yield_prod_p,rel_widths = c(1,1.45))

#Shoot biomass

all_site_shoot_biomass_mod=lmer(log(total.shoot.dry.weight)~FertStatus*siteID+(1|plotRep),data = MMPRNT_2018.metadata_mapp_all_sites_G5)
plot(all_site_shoot_biomass_mod)
hist(resid(all_site_shoot_biomass_mod))
qqPlot(resid(all_site_shoot_biomass_mod))
shapiro.test(resid(all_site_shoot_biomass_mod))
##W = 0.98515, p-value = 0.2417
anova(all_site_shoot_biomass_mod)

emmeans(all_site_shoot_biomass_mod,pairwise~siteID,adjust="fdr")

shoot_bio_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                             sig_let=c("a","b","bd","c","cd"))
emmeans(all_site_shoot_biomass_mod,pairwise~FertStatus|siteID)


#Shoot biomass
MMPRNT_2018.metadata_mapp_all_sites_G5_yield=MMPRNT_2018.metadata_mapp_all_sites_G5%>%group_by(FertStatus,siteID,plotRep)%>%
  summarise(dry_matter_yield_mg_ha_mean=mean(dry_matter_yield_mg_ha_mean))
all_site_yield_mod=lmer((dry_matter_yield_mg_ha_mean)~FertStatus*siteID+(1|plotRep),data = MMPRNT_2018.metadata_mapp_all_sites_G5_yield)
plot(all_site_yield_mod)
hist(resid(all_site_yield_mod))
qqPlot(resid(all_site_yield_mod))
shapiro.test(resid(all_site_yield_mod))
#W = 0.94341, p-value = 0.05409
anova(all_site_yield_mod)

emmeans(all_site_yield_mod,pairwise~siteID,adjust="fdr")

yield_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                             sig_let=c("a","b","c","c","c"))
emmeans(all_site_yield_mod,pairwise~FertStatus|siteID)

#Soil abiotic variables 


colnames(MMPRNT_2018.metadata_mapp_all_sites_G5)

res.pca_v4 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                       "percent_soil_moisture_dry_weight",
                                                                       "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")]), scale = TRUE)

fviz_eig(res.pca_v4)

fviz_pca_ind(res.pca_v4, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Soil Nutrients")


fviz_pca_ind(res.pca_v4, habillage=
               na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("siteID","ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                 "percent_soil_moisture_dry_weight",
                                                                 "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")])$siteID,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Soil Nutrients")


fviz_pca_var(res.pca_v4,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Soil Nutrients")


fviz_pca_biplot(res.pca_v4, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)+ggtitle("Soil Nutrients")

fviz_pca_biplot(res.pca_v4, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5[,c("siteID","ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                    "percent_soil_moisture_dry_weight",
                                                                    "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")])$siteID,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Soil Nutrients")
MMPRNT_2018.metadata_mapp_all_sites_G5_soil=MMPRNT_2018.metadata_mapp_all_sites_G5%>%group_by(FertStatus,siteID,plotRep)%>%
  summarise(p_ppm=mean(p_ppm),ph=mean(ph))

value_max_soil_P<-MMPRNT_2018.metadata_mapp_all_sites_G5 %>% group_by(siteID) %>% summarize(max_value = max(p_ppm))
value_max_soil_pH<-MMPRNT_2018.metadata_mapp_all_sites_G5 %>% group_by(siteID) %>% summarize(max_value = max(ph))


(soil_phos_p<-ggplot(MMPRNT_2018.metadata_mapp_all_sites_G5_soil,aes(x=siteID, y=p_ppm))+
    geom_boxplot(data = MMPRNT_2018.metadata_mapp_all_sites_G5_soil,aes(x=siteID, y=p_ppm,color=FertStatus))+
    geom_text(data = value_max_soil_P, aes(x=siteID, y = 0.1 + max_value, label = soil_P_letters$sig_let), vjust=0, size=10)+
  theme_classic(base_size = 25)+theme(legend.position = "none"))

(soil_PH_p<-ggplot(MMPRNT_2018.metadata_mapp_all_sites_G5_soil,aes(x=siteID, y=ph))+
    geom_boxplot(data = MMPRNT_2018.metadata_mapp_all_sites_G5_soil,aes(x=siteID, y=ph,color=FertStatus))+
    geom_text(data = value_max_soil_pH, aes(x=siteID, y = 0.1 + max_value, label = soil_pH_letters$sig_let), vjust=0, size=10)+
    theme_classic(base_size = 25))

plot_grid(soil_phos_p,soil_PH_p,rel_widths = c(1,1.45))
#1400x600
#Soil phosporus

all_site_soil_P_mod=lmer((p_ppm)~FertStatus*siteID+(1|plotRep),data = MMPRNT_2018.metadata_mapp_all_sites_G5_soil)
plot(all_site_soil_P_mod)
hist(resid(all_site_soil_P_mod))
qqPlot(resid(all_site_soil_P_mod))
shapiro.test(resid(all_site_soil_P_mod))
#W = 0.97191, p-value = 0.4455
anova(all_site_soil_P_mod)


emmeans(all_site_soil_P_mod,pairwise~siteID,adjust="fdr")

soil_P_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                             sig_let=c("a","b","c","ac","d"))


#Soil pH


all_site_soil_pH_mod=lmer((ph)~FertStatus*siteID+(1|plotRep),data = MMPRNT_2018.metadata_mapp_all_sites_G5_soil)
plot(all_site_soil_pH_mod)
hist(resid(all_site_soil_pH_mod))
qqPlot(resid(all_site_soil_pH_mod))
shapiro.test(resid(all_site_soil_pH_mod))
#W = 0.9772, p-value = 0.618
anova(all_site_soil_pH_mod)


emmeans(all_site_soil_pH_mod,pairwise~siteID,adjust="fdr")

soil_pH_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                          sig_let=c("a","b","c","d","e"))


#Microbial traits
MMPRNT_2018.metadata_mapp_all_sites_G5_Msub=subset(MMPRNT_2018.metadata_mapp_all_sites_G5,MBN_NA<300)

colnames(MMPRNT_2018.metadata_mapp_all_sites_G5_Msub)

res.pca_v5 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5_Msub[,c("ugN_fixed_gdrysoil_day","TOC","TON","MBC_NA","MBN_NA",
                                                                       "MBC_MBN_ratio","TOC_TON_ratio","potential_NH4_mineralization_rate_na",
                                                                       "potential_NO3_nitrification_rate_na"
)]), scale = TRUE)

fviz_eig(res.pca_v5)

fviz_pca_ind(res.pca_v5, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Microbial traits")


fviz_pca_ind(res.pca_v5, habillage=
               na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5_Msub[,c("siteID","ugN_fixed_gdrysoil_day","TOC","TON","MBC_NA","MBN_NA",
                                                                 "MBC_MBN_ratio","TOC_TON_ratio","potential_NH4_mineralization_rate_na",
                                                                 "potential_NO3_nitrification_rate_na")])$siteID,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Microbial traits")


fviz_pca_var(res.pca_v5,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Plant data")


fviz_pca_biplot(res.pca_v5, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)+ggtitle("Microbial traits")

fviz_pca_biplot(res.pca_v5, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  na.omit(MMPRNT_2018.metadata_mapp_all_sites_G5_Msub[,c("siteID","ugN_fixed_gdrysoil_day","TOC","TON","MBC_NA","MBN_NA",
                                                                    "MBC_MBN_ratio","TOC_TON_ratio","potential_NH4_mineralization_rate_na",
                                                                    "potential_NO3_nitrification_rate_na")])$siteID,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Microbial traits")


value_max_N_min<-MMPRNT_2018.metadata_mapp_all_sites_G5_Msub %>% group_by(siteID) %>% 
  summarize(max_value = max(potential_NH4_mineralization_rate_na,na.rm = T))
value_max_MBN<-MMPRNT_2018.metadata_mapp_all_sites_G5_Msub %>% group_by(siteID) %>% summarize(max_value = max(MBN_NA))


(N_min_p<-ggplot(MMPRNT_2018.metadata_mapp_all_sites_G5_Msub,aes(x=siteID, y=potential_NH4_mineralization_rate_na))+
    geom_boxplot(data = MMPRNT_2018.metadata_mapp_all_sites_G5_Msub,aes(x=siteID, y=potential_NH4_mineralization_rate_na,color=FertStatus))+
    geom_text(data = value_max_N_min, aes(x=siteID, y = 0.1 + max_value, label = N_min_letters$sig_let), vjust=0, size=10)+
    geom_signif(y_position = c(-0.02), xmin = c(1.85), xmax = c(2.15),annotation = c("**"), 
                tip_length = 0,textsize = 8)+theme_classic(base_size = 25)+theme(legend.position = "none"))

(MBN_p<-ggplot(MMPRNT_2018.metadata_mapp_all_sites_G5_Msub,aes(x=siteID, y=MBN_NA))+
    geom_boxplot(data = MMPRNT_2018.metadata_mapp_all_sites_G5_Msub,aes(x=siteID, y=MBN_NA,color=FertStatus))+
    geom_text(data = value_max_MBN, aes(x=siteID, y = 0.1 + max_value, label = microb_N_letters$sig_let), vjust=0, size=10)+
    theme_classic(base_size = 25))

plot_grid(N_min_p,MBN_p,rel_widths = c(1,1.45))

#Nit mineralization

all_site_N_min_mod=lmer((potential_NH4_mineralization_rate_na)~FertStatus*siteID+(1|plotRep),data = MMPRNT_2018.metadata_mapp_all_sites_G5_Msub)
plot(all_site_N_min_mod)
hist(resid(all_site_N_min_mod))
qqPlot(resid(all_site_N_min_mod))
shapiro.test(resid(all_site_N_min_mod))
#W = 0.91577, p-value = 3.677e-06
anova(all_site_N_min_mod)

emmeans(all_site_N_min_mod,pairwise~siteID,adjust="fdr")

N_min_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                             sig_let=c("a","b","c","c","d"))
emmeans(all_site_N_min_mod,pairwise~FertStatus|siteID)


#Microbial nitrogen

all_site_microb_N_mod=lmer(sqrt(MBN_NA)~FertStatus*siteID+(1|plotRep),data = MMPRNT_2018.metadata_mapp_all_sites_G5_Msub)
plot(all_site_microb_N_mod)
hist(resid(all_site_microb_N_mod))
qqPlot(resid(all_site_microb_N_mod))
shapiro.test(resid(all_site_microb_N_mod))
#W = 0.97517, p-value = 0.03483
anova(all_site_microb_N_mod)

emmeans(all_site_microb_N_mod,pairwise~siteID,adjust="fdr")

microb_N_letters=data.frame(siteID=c("ESC","HAN","LC","LUX","RHN"),
                         sig_let=c("a","b","c","c","c"))


#####LUX G5####



MMPRNT_2018.metadata_mapp_coverage_LUX=subset(MMPRNT_2018.metadata_mapp_coverage,siteID=="LUX"&plotType=="G5")
nrow(MMPRNT_2018.metadata_mapp_coverage_LUX)
#360


colSums(is.na(MMPRNT_2018.metadata_mapp_coverage_LUX))

unique(MMPRNT_2018.metadata_mapp_coverage_LUX$siteID)
#"LC"  "HAN" "RHN" "ESC" "LUX"
unique(MMPRNT_2018.metadata_mapp_coverage_LUX$plotType)



MMPRNT_2018.metadata_mapp_LUX=MMPRNT_2018.metadata_mapp_coverage_LUX
MMPRNT_2018.metadata_mapp_LUX[,c("event_age","pH_raw_MMPRNT","pH_MMPRNT_subplot")]=NULL

colnames(MMPRNT_2018.metadata_mapp_LUX)


#PCA 
summary(MMPRNT_2018.metadata_mapp_LUX[,13:70])

#need to remove NA's
unique(subset(MMPRNT_2018.metadata_mapp_LUX,!is.na(vwc_avg_mean))$collectionDate)
MMPRNT_2018.metadata_mapp_LUX%>%group_by(collectionDate)%>%summarise(sum_na=sum(is.na(total.shoot.dry.weight)))
#3 11/5/2018          24
#4 3/19/2018          24
#5 4/30/2018          24

MMPRNT_2018.metadata_mapp_LUX%>%group_by(collectionDate)%>%summarise(sum_na=sum(is.na(plant.height)))
#3 11/5/2018          24
#4 3/19/2018          24
#5 4/30/2018          24
#6 5/15/2018          24

MMPRNT_2018.metadata_mapp_LUX%>%group_by(collectionDate)%>%summarise(sum_na=sum(is.na(core_root_mass_subplot)))
#1 10/15/2018          3
#2 10/3/2018           3
#3 11/5/2018           6
#4 3/19/2018           6
#5 4/30/2018           3
#7 5/29/2018           3
#8 6/11/2018           3
#10 7/30/2018           6
#11 7/9/2018            3
#12 8/20/2018           3
#15 9/4/2018            9

MMPRNT_2018.metadata_mapp_LUX%>%group_by(collectionDate)%>%summarise(sum_na=sum(is.na(vwc_avg_mean)))
#3/19/2018          24

MMPRNT_2018.metadata_mapp_LUX%>%group_by(collectionDate)%>%summarise(sum_na=sum(is.na(MBC_MBN_ratio)))
#1 10/15/2018          5
#3 11/5/2018           1
#4 3/19/2018           5
#5 4/30/2018           3
#6 5/15/2018           3
#8 6/11/2018           7
#9 6/25/2018           7
#10 7/30/2018           4
#11 7/9/2018            1
#14 9/17/2018           1



MMPRNT_2018.metadata_mapp_LUX%>%group_by(collectionDate)%>%summarise(sum_na=sum(is.na(TOC_TON_ratio)))
#1 10/15/2018          3
#3 11/5/2018           1
#4 3/19/2018           5
#5 4/30/2018           2
#6 5/15/2018           4
#8 6/11/2018           9
#9 6/25/2018           5
#11 7/9/2018            1
#12 8/20/2018           3
#13 8/8/2018            1


row.names(MMPRNT_2018.metadata_mapp_LUX)=MMPRNT_2018.metadata_mapp_LUX$SampleID

res.pca_LUX <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(13:70)]), scale = TRUE)

fviz_eig(res.pca_LUX)



fviz_pca_ind(res.pca_LUX, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_ind(res.pca_LUX, 
             habillage=as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(2,13:70)])$collectionDate,format = "%m/%d/%Y"),
             addEllipses=TRUE, ellipse.level=0.95, palette = viridis(15)
)


fviz_pca_var(res.pca_LUX,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca_LUX, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)

fviz_pca_biplot(res.pca_LUX, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(2,13:70)])$collectionDate,format = "%m/%d/%Y"),
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)

fviz_pca_biplot(res.pca_LUX, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(6,13:70)])$plotRep,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)

colnames(MMPRNT_2018.metadata_mapp_LUX)

(lux_par_tot_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=par_den_avg_mean))+
    geom_line()+geom_point(shape=8,size=6)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(legend.position = "none"))

(lux_par_tot_CV_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=par_den_avg_CV))+
    geom_line()+geom_point(shape=8,size=6)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(legend.position = "none"))

(lux_soil_temp_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=soil_temp_1_avg_mean))+
    geom_line()+geom_point(shape=8,size=6)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(legend.position = "none"))

(lux_shoot_water_pot_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),
                                                                 y=shoot.water.percentage,color=FertStatus))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme())

plot_grid(lux_par_tot_p,lux_par_tot_CV_p,lux_soil_temp_p,lux_shoot_water_pot_p,nrow = 2,labels = c('a)', 'b)','c)','d)'), label_size = 30)

#Let's remove remove plant data and and met data to see as many dates as possible
colnames(MMPRNT_2018.metadata_mapp_LUX)
colnames(MMPRNT_2018.metadata_mapp_LUX[,c(13:14,30,52:64,66:69)])
res.pca_LUX_v2 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(13:14,30,52:64,66:69)]), scale = TRUE)

fviz_eig(res.pca_LUX_v2)

fviz_pca_ind(res.pca_LUX_v2, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_ind(res.pca_LUX_v2, habillage=as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(2,13:14,30,52:64,66:69)])$collectionDate,format = "%m/%d/%Y"),
             addEllipses=TRUE, ellipse.level=0.95, palette = viridis(15)
)


fviz_pca_var(res.pca_LUX_v2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_biplot(res.pca_LUX_v2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)

fviz_pca_biplot(res.pca_LUX_v2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(2,13:14,30,52:64,66:69)])$collectionDate,format = "%m/%d/%Y"),
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)

fviz_pca_biplot(res.pca_LUX_v2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=na.omit(MMPRNT_2018.metadata_mapp_LUX[,c(6,13:14,30,52:64,66:69)])$plotRep,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)

#Plant data only


colnames(MMPRNT_2018.metadata_mapp_LUX)

res.pca_LUX_v3 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("dry_matter_yield_mg_ha_mean","total.shoot.dry.weight",
                                                                       "specific.leaf.area","plant.height","core_root_mass_subplot"
)]), scale = TRUE)

fviz_eig(res.pca_LUX_v3)

fviz_pca_ind(res.pca_LUX_v3, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Plant data")


fviz_pca_ind(res.pca_LUX_v3, habillage=
               as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","dry_matter_yield_mg_ha_mean","total.shoot.dry.weight",
                                                                "specific.leaf.area","plant.height","core_root_mass_subplot")])$collectionDate,format = "%m/%d/%Y"),
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Plant data")


fviz_pca_var(res.pca_LUX_v3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Plant data")


fviz_pca_biplot(res.pca_LUX_v3, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)+ggtitle("Plant data")

fviz_pca_biplot(res.pca_LUX_v3, repel = TRUE,
                col.var = "black", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","dry_matter_yield_mg_ha_mean","total.shoot.dry.weight",
                                                                    "specific.leaf.area","plant.height","core_root_mass_subplot")])$collectionDate,format = "%m/%d/%Y"),
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Plant data")

fviz_pca_biplot(res.pca_LUX_v3, repel = TRUE,
                col.var = "black", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("plotRep","dry_matter_yield_mg_ha_mean","total.shoot.dry.weight",
                                                                   "specific.leaf.area","plant.height","core_root_mass_subplot")])$plotRep,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Plant data")

#value_max_shoot_bio<-MMPRNT_2018.metadata_mapp_LUX %>% group_by(siteID) %>% summarize(max_value = max(total.shoot.dry.weight))
#value_max_yield<-MMPRNT_2018.metadata_mapp_all_sites_G5 %>% group_by(siteID) %>% summarize(max_value = max(dry_matter_yield_mg_ha_mean))


pairwise_fert_shootbio_lux=data.frame(y_bot=c(11.5,1,16.5,13,12),
                                       x_min=as.Date(c("10/11/2018","5/11/2018","7/26/2018","8/4/2018","9/13/2018"), format="%m/%d/%Y"),
                                       x_max=as.Date(c("10/19/2018","5/19/2018","8/3/2018","8/12/2018","9/21/2018"), format="%m/%d/%Y"),
                                       annot_text=c("**"," ** ","*","#","  **  "))

#

(lux_yield_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=FertStatus,y=dry_matter_yield_mg_ha_mean))+
    geom_boxplot()+geom_text(aes(label=plotRep, color=factor(plotRep)),size=20)+
    theme_cowplot(font_size = 24)+scale_x_discrete(name = NULL)+theme(legend.position = "none"))
#

(lux_shoot_biomass_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=total.shoot.dry.weight))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus),color=FertStatus),position = position_dodge(8))+
    geom_signif(data=pairwise_fert_shootbio_lux,mapping = aes(y_position = y_bot, xmin = x_min, xmax = x_max,annotations = annot_text), 
                tip_length = 0,textsize = 8,manual = T)+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(strip.text = element_blank(),strip.background = element_blank()))

(lux_shoot_biomass_facet_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=total.shoot.dry.weight,color=FertStatus))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus)),position = position_dodge(8))+
    facet_wrap(~plotRep,nrow = 2,scales = "free_x")+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme())

#1700x1200
(top_row_bio=plot_grid(lux_yield_p,lux_shoot_biomass_p,nrow = 1,labels = c('a)', 'b)'), label_size = 30))
plot_grid(top_row_bio,lux_shoot_biomass_facet_p,nrow = 2,labels = c('', 'c)'), label_size = 30)





#Shoot biomass

lux_shoot_biomass_mod=lmer(log(total.shoot.dry.weight)~FertStatus*collectionDate+(1|plotRep),data = MMPRNT_2018.metadata_mapp_LUX)
plot(lux_shoot_biomass_mod)
hist(resid(lux_shoot_biomass_mod))
qqPlot(resid(lux_shoot_biomass_mod))
shapiro.test(resid(lux_shoot_biomass_mod))
#W = 0.99532, p-value = 0.5414
anova(lux_shoot_biomass_mod)

emmeans(lux_shoot_biomass_mod,pairwise~collectionDate,adjust="fdr")


emmeans(lux_shoot_biomass_mod,pairwise~FertStatus|collectionDate)
#collectionDate = 10/15/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert   0.4048 0.148 259  2.739  0.0066 

#collectionDate = 5/15/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert   0.3983 0.148 259  2.695  0.0075

#collectionDate = 7/30/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert   0.3643 0.148 259  2.465  0.0144 

#collectionDate = 8/8/2018:
#contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert   0.2616 0.148 259  1.770  0.0779 

#collectionDate = 9/17/2018:
#  contrast      estimate    SE  df t.ratio p.value
#Fert - Unfert   0.4497 0.151 259  2.976  0.0032 



#yeild
MMPRNT_2018.metadata_mapp_LUX_yield=MMPRNT_2018.metadata_mapp_LUX%>%group_by(FertStatus,plotRep)%>%
  summarise(dry_matter_yield_mg_ha_mean=mean(dry_matter_yield_mg_ha_mean))
LUX_yield_mod=lmer((dry_matter_yield_mg_ha_mean)~FertStatus+(1|plotRep),data = MMPRNT_2018.metadata_mapp_LUX_yield)
plot(LUX_yield_mod)
hist(resid(LUX_yield_mod))
qqPlot(resid(LUX_yield_mod))
shapiro.test(resid(LUX_yield_mod))
#W = 0.9469, p-value = 0.6799
anova(LUX_yield_mod)
#FertStatus  7.263   7.263     1 2.9992  0.9687 0.3976


#Soil abiotic variables 


colnames(MMPRNT_2018.metadata_mapp_LUX)

res.pca_LUX_v4 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                       "percent_soil_moisture_dry_weight",
                                                                       "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")]), scale = TRUE)

fviz_eig(res.pca_LUX_v4)

fviz_pca_ind(res.pca_LUX_v4, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Soil Nutrients")


fviz_pca_ind(res.pca_LUX_v4, habillage=
               as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                "percent_soil_moisture_dry_weight",
                                                                "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")])$collectionDate,format = "%m/%d/%Y"),
             addEllipses=TRUE, ellipse.level=0.95, palette = viridis(15)
)+ggtitle("Soil Nutrients")


fviz_pca_var(res.pca_LUX_v4,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Soil Nutrients")


fviz_pca_biplot(res.pca_LUX_v4, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)+ggtitle("Soil Nutrients")

fviz_pca_biplot(res.pca_LUX_v4, repel = TRUE,
                col.var = "black", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                   "percent_soil_moisture_dry_weight",
                                                                   "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")])$collectionDate,
                          format = "%m/%d/%Y"),
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Soil Nutrients")


fviz_pca_biplot(res.pca_LUX_v4, repel = TRUE,
                col.var = "black", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                 na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("plotRep","ph","p_ppm","k_ppm","ca_ppm","mg_ppm",
                                                                   "percent_soil_moisture_dry_weight",
                                                                   "ugN_NH4_g_dry_soil_na","ugN_NO3_g_dry_soil_na")])$plotRep,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Soil Nutrients")


MMPRNT_2018.metadata_mapp_LUX_soil=MMPRNT_2018.metadata_mapp_LUX%>%group_by(FertStatus,plotRep)%>%
  summarise(p_ppm=mean(p_ppm),ph=mean(ph))




(soil_phos_LUX_p<-ggplot(MMPRNT_2018.metadata_mapp_LUX_soil,aes(x=(FertStatus), y=p_ppm))+
    geom_boxplot(data = MMPRNT_2018.metadata_mapp_LUX_soil,aes(x=FertStatus, y=p_ppm))+
    geom_text(aes(label=plotRep, color=factor(plotRep)),size=20)+
    theme_classic(base_size = 25)+theme(legend.position = "none",axis.title.x = element_blank()))

(soil_PH_LUX_p<-ggplot(MMPRNT_2018.metadata_mapp_LUX_soil,aes(x=FertStatus, y=ph))+
    geom_boxplot(data = MMPRNT_2018.metadata_mapp_LUX_soil,aes(x=FertStatus, y=ph,))+
    geom_text(aes(label=plotRep, color=factor(plotRep)),size=20)+
    theme_classic(base_size = 25)+theme(legend.position = "none",axis.title.x = element_blank()))

plot_grid(soil_phos_LUX_p,soil_PH_LUX_p,rel_widths = c(1,1))
#1400x600
#Soil phosporus

LUX_soil_P_mod=lmer((p_ppm)~FertStatus+(1|plotRep),data = MMPRNT_2018.metadata_mapp_LUX_soil)
plot(LUX_soil_P_mod)
hist(resid(LUX_soil_P_mod))
qqPlot(resid(LUX_soil_P_mod))
shapiro.test(resid(LUX_soil_P_mod))
#W = 0.93144, p-value = 0.5293
anova(LUX_soil_P_mod)
#FertStatus   12.5    12.5     1     3  0.9259 0.4069

emmeans(LUX_soil_P_mod,pairwise~siteID,adjust="fdr")




#Soil pH


LUX_soil_pH_mod=lmer((ph)~FertStatus+(1|plotRep),data = MMPRNT_2018.metadata_mapp_LUX_soil)
plot(LUX_soil_pH_mod)
hist(resid(LUX_soil_pH_mod))
qqPlot(resid(LUX_soil_pH_mod))
shapiro.test(resid(LUX_soil_pH_mod))
#W = 0.93152, p-value = 0.53
anova(LUX_soil_pH_mod)
#FertStatus 0.03125 0.03125     1     3  2.1429 0.2394




#Microbial traits
#MMPRNT_2018.metadata_mapp_all_sites_G5_Msub=subset(MMPRNT_2018.metadata_mapp_LUX,MBN_NA<300)

colnames(MMPRNT_2018.metadata_mapp_LUX)

res.pca_LUX_v5 <- prcomp(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("ugN_fixed_gdrysoil_day","TOC","TON","MBC_NA","MBN_NA",
                                                                            "MBC_MBN_ratio","TOC_TON_ratio","potential_NH4_mineralization_rate_na",
                                                                            "potential_NO3_nitrification_rate_na"
)]), scale = TRUE)

fviz_eig(res.pca_LUX_v5)

fviz_pca_ind(res.pca_LUX_v5, 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Microbial traits")


fviz_pca_ind(res.pca_LUX_v5, habillage=
               as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","ugN_fixed_gdrysoil_day","TOC",
                                                                              "TON","MBC_NA","MBN_NA",
                                                                      "MBC_MBN_ratio","TOC_TON_ratio",
                                                                      "potential_NH4_mineralization_rate_na",
                                                                      "potential_NO3_nitrification_rate_na")])$collectionDate,
                       format="%m/%d/%Y"),
             addEllipses=TRUE, ellipse.level=0.95, palette = viridis(15)
)+ggtitle("Microbial traits")


fviz_pca_var(res.pca_LUX_v5,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("Plant data")


fviz_pca_biplot(res.pca_LUX_v5, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20)
)+ggtitle("Microbial traits")

fviz_pca_biplot(res.pca_LUX_v5, repel = TRUE,
                col.var = "black", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  as.Date(na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("collectionDate","ugN_fixed_gdrysoil_day","TOC",
                                                                   "TON","MBC_NA","MBN_NA",
                                                                   "MBC_MBN_ratio","TOC_TON_ratio",
                                                                   "potential_NH4_mineralization_rate_na",
                                                                   "potential_NO3_nitrification_rate_na")])$collectionDate,
                          format="%m/%d/%Y"),
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Microbial traits")

fviz_pca_biplot(res.pca_LUX_v5, repel = TRUE,
                col.var = "black", # Variables color
                col.ind = "#696969",  # Individuals color,
                select.var = list(contrib = 20),label = "var",
                habillage=
                  na.omit(MMPRNT_2018.metadata_mapp_LUX[,c("plotRep","ugN_fixed_gdrysoil_day","TOC",
                                                                   "TON","MBC_NA","MBN_NA",
                                                                   "MBC_MBN_ratio","TOC_TON_ratio",
                                                                   "potential_NH4_mineralization_rate_na",
                                                                   "potential_NO3_nitrification_rate_na")])$plotRep,
                addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2"
)+ggtitle("Microbial traits")



#value_max_N_min<-MMPRNT_2018.metadata_mapp_all_sites_G5_Msub %>% group_by(siteID) %>% 
#  summarize(max_value = max(potential_NH4_mineralization_rate_na,na.rm = T))
#value_max_MBN<-MMPRNT_2018.metadata_mapp_all_sites_G5_Msub %>% group_by(siteID) %>% summarize(max_value = max(MBN_NA))



(TOC_lux_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=TOC))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus),color=FertStatus),position = position_dodge(8))+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(legend.position = "none"))



(MBC_lux_p=ggplot(MMPRNT_2018.metadata_mapp_LUX, aes(x=as.Date(collectionDate, format="%m/%d/%Y"),y=MBC_NA))+
    geom_boxplot(aes(group=interaction(collectionDate,FertStatus),color=FertStatus),position = position_dodge(8))+
    theme_cowplot(font_size = 24)+scale_x_date(name = NULL)+theme(strip.text = element_blank(),strip.background = element_blank()))

plot_grid(TOC_lux_p,MBC_lux_p,rel_widths = c(1,1.45))
#1400x600


#Total organic carbon

lux_TOC_mod=lmer((TOC)^(-1/2)~FertStatus*collectionDate+(1|plotRep),data = MMPRNT_2018.metadata_mapp_LUX)
plot(all_site_N_min_mod)
hist(resid(lux_TOC_mod))
qqPlot(resid(lux_TOC_mod))
shapiro.test(resid(lux_TOC_mod))
#W = 0.98865, p-value = 0.01068
anova(lux_TOC_mod)
#FertStatus                0.000284 0.0002842     1 300.05  0.9351 0.3343    
#collectionDate            0.067515 0.0048225    14 300.17 15.8697 <2e-16 ***
#FertStatus:collectionDate 0.004478 0.0003199    14 300.07  1.0527 0.4009 



#Microbial nitrogen

LUX_microb_C_mod=lmer((MBC_NA)~FertStatus*collectionDate+(1|plotRep),data = MMPRNT_2018.metadata_mapp_LUX)
plot(LUX_microb_C_mod)
hist(resid(LUX_microb_C_mod))
qqPlot(resid(LUX_microb_C_mod))
shapiro.test(resid(LUX_microb_C_mod))
#W = 0.99268, p-value = 0.1137
anova(LUX_microb_C_mod)
#collectionDate            3810497  272178    14 291.10 13.1712 <2e-16 ***
emmeans(LUX_microb_C_mod,pairwise~FertStatus|collectionDate)


