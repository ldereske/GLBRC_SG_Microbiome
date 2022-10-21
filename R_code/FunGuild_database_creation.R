library(here)
library(stringr)

#####Construction of database for FunGuild####



#FunGuild of CONSTAX
GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus = read.delim(here::here("Publish_data",
                                                                     "constax_taxonomy.txt"),
                                                          sep = "\t",header = T,fill=T, row.names = 1)
head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)
GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus[GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus==""]="Unknown"

GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Species=
  str_replace(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Species," ","_")
GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$taxonomy=
  paste(paste("d",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Kingdom,sep="__"),
        paste("p",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Phylum,sep="__"),
        paste("c",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Class,sep="__"),
        paste("o",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Order,sep="__"),
        paste("f",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Family,sep="__"),
        paste("g",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Genus,sep="__"),
        paste("s",GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus$Species,sep="__"),sep = ";")

head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)
tail(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus)

GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw=data.frame("OTU ID"=row.names(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus),
                                                                     GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus[,c("Kingdom","taxonomy")])

head(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw)
write.table(GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw,
            here::here("Publish_data","GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw.txt"),row.names = F,sep = '\t')

#Code ran in python (https://github.com/UMNFuN/FUNGuild)

#python FUNGuild_1.py taxa -otu  GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw.txt -format tsv -column taxonomy -classifier unite
#there was an encoding error so I had to add  encoding="utf-8" to line 214 of FUNGuild.py
#python FUNGuild_1.py guild -taxa GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw.taxa.txt
#####Found 11284 records in the db.




#Load the FunGuild classifications


GLBRC_CON_FunGuild_raw=read.table(here::here("R_files","Publish_data","GLBRC018_fung_raw_CONSTAX_UNITE8.2_consensus_FunGuild_raw.taxa.guilds.txt"),
                                  sep = '\t',header = T)
head(GLBRC_CON_FunGuild_raw)
unique(GLBRC_CON_FunGuild_raw$confidenceRanking)
#Only Possible categories or higher!
GLBRC_CON_FunGuild_raw_F=subset(GLBRC_CON_FunGuild_raw,confidenceRanking!="Possible"&confidenceRanking!="na")
dim(GLBRC_CON_FunGuild_raw_F)
#2432   17

unique(GLBRC_CON_FunGuild_raw_F$confidenceRanking)
unique(GLBRC_CON_FunGuild_raw_F$guild)
#Creation of a summarized guild from taxa with multiple classications. 
GLBRC_CON_FunGuild_raw_F$simp_guild=
  with(GLBRC_CON_FunGuild_raw_F,
       ifelse(guild=="Endophyte-Litter Saprotroph-Soil Saprotroph-Undefined Saprotroph"|
                guild=="Lichenized-Undefined Saprotroph"|
                guild=="Ectomycorrhizal-Undefined Saprotroph"|
                guild=="Undefined Saprotroph-Undefined Symbiotroph"|
                guild=="Endophyte-Undefined Saprotroph"|
                guild=="Epiphyte-Undefined Saprotroph"|
                guild=="Lichenized-Wood Saprotroph"|
                guild=="Dung Saprotroph-Ectomycorrhizal", 
              "Symbiotroph-Saprotroph",
              ifelse(guild=="Endophyte-Lichen Parasite-Plant Pathogen-Undefined Saprotroph"|
                       guild=="Endophyte-Plant Pathogen-Undefined Saprotroph"|
                       guild=="Endophyte-Plant Pathogen-Wood Saprotroph"|
                       guild=="Animal Pathogen-Endophyte-Epiphyte-Fungal Parasite-Plant Pathogen-Wood Saprotroph"|
                       guild=="Ectomycorrhizal-Fungal Parasite-Plant Saprotroph-Wood Saprotroph"|
                       guild=="Animal Endosymbiont-Animal Pathogen-Endophyte-Plant Pathogen-Undefined Saprotroph"|
                       guild=="Animal Endosymbiont-Animal Pathogen-Undefined Saprotroph"|
                       guild=="Ectomycorrhizal-Fungal Parasite-Plant Pathogen-Wood Saprotroph", 
                     "Symbiotroph-Pathogen-Saprotroph",
                     ifelse(guild=="Animal Pathogen-Undefined Saprotroph"|
                              guild=="Animal Pathogen-Fungal Parasite-Undefined Saprotroph"|
                              guild=="Plant Pathogen-Wood Saprotroph"|
                              guild=="Dung Saprotroph-Plant Parasite-Soil Saprotroph-Undefined Saprotroph-Wood Saprotroph"|
                              guild=="Fungal Parasite-Plant Pathogen-Plant Saprotroph"|
                              guild=="Plant Pathogen-Plant Saprotroph"|
                              guild=="Plant Pathogen-Undefined Parasite-Undefined Saprotroph"|
                              guild=="Dung Saprotroph-Nematophagous"|
                              guild=="Bryophyte Parasite-Litter Saprotroph-Wood Saprotroph"|
                              guild=="Plant Pathogen-Undefined Saprotroph"|
                              guild=="Bryophyte Parasite-Leaf Saprotroph-Soil Saprotroph-Undefined Saprotroph-Wood Saprotroph"|
                              guild=="Fungal Parasite-Litter Saprotroph"|
                              guild=="Dung Saprotroph-Endophyte-Plant Pathogen-Undefined Saprotroph"|
                              guild=="Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"|
                              guild=="Algal Parasite-Bryophyte Parasite-Fungal Parasite-Undefined Saprotroph"|
                              guild=="Algal Parasite-Fungal Parasite-Undefined Saprotroph"|
                              guild=="Animal Pathogen-Plant Pathogen-Undefined Saprotroph",
                            "Pathogen-Saprotroph",
                            ifelse(guild=="Endophyte-Plant Pathogen"|
                                     guild=="Ectomycorrhizal-Fungal Parasite"|
                                     guild=="Endophyte-Insect Pathogen"|
                                     guild=="Lichen Parasite-Lichenized",
                                   "Symbiotroph-Pathogen",
                                   ifelse(guild=="Soil Saprotroph-Undefined Saprotroph"|
                                            guild=="Dung Saprotroph-Soil Saprotroph-Undefined Saprotroph"|
                                            guild=="Plant Saprotroph-Wood Saprotroph"|
                                            guild=="Dung Saprotroph-Plant Saprotroph-Soil Saprotroph"|
                                            guild=="Dung Saprotroph-Plant Saprotroph"|
                                            guild=="Dung Saprotroph-Undefined Saprotroph"|
                                            guild=="Dung Saprotroph-Plant Saprotroph-Wood Saprotroph"|
                                            guild=="Dung Saprotroph-Soil Saprotroph"|
                                            guild=="Dung Saprotroph-Wood Saprotroph"|
                                            guild=="Dung Saprotroph-Soil Saprotroph-Wood Saprotrop"|
                                            guild=="Litter Saprotroph-Soil Saprotroph-Wood Saprotroph"|
                                            guild=="Leaf Saprotroph-Wood Saprotroph"|
                                            guild=="Dung Saprotroph-Leaf Saprotroph", 
                                          "Multiple Saprotroph",
                                          ifelse(guild=="Animal Parasite-Fungal Parasite","Multiple Pathogen",
                                                 guild)))))))

unique(GLBRC_CON_FunGuild_raw_F$simp_guild)

write.csv(GLBRC_CON_FunGuild_raw_F, here::here("Publish_data","GLBRC_CON_FunGuild_raw_F.csv"),
          row.names = F)
