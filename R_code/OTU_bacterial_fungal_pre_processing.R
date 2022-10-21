
library(here)
library(phyloseq)
library(plyr); library(dplyr)
library(reshape2)
library(tidyr)
library(seqinr)
library(ggplot2)
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
nrow(MMPRNT_16s.metadata_mapp_raw)
MMPRNT_16s.metadata_mapp_raw$treatment_GLBRC=paste(MMPRNT_16s.metadata_mapp_raw$plotType,MMPRNT_16s.metadata_mapp_raw$plotRep,sep = "R")

#Add in the missing gps ponts

MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correct2=read.csv(here::here("ABIOTIC_data","MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correction_coef.csv"))
MLE_coord2_m2_RHN_G5_UTM_correct2=subset(MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correct2,site=="Rhinelander")
MLE_coord2_m2_HAN_G5_UTM_correct2=subset(MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correct2,site=="Hancock")


#Add in the gps points for Rhinelander 

MMPRNT_16s.metadata_mapp_raw_rhn=subset(MMPRNT_16s.metadata_mapp_raw, siteID=="RHN")
nrow(MMPRNT_16s.metadata_mapp_raw_rhn)
#216
MMPRNT_16s.metadata_mapp_raw_rhn_correction=merge(MMPRNT_16s.metadata_mapp_raw_rhn,MLE_coord2_m2_RHN_G5_UTM_correct2, by.x="treatment_GLBRC",by.y="treatment" )

MMPRNT_16s.metadata_mapp_raw_rhn_correction$UTM_Lon_Cord=MMPRNT_16s.metadata_mapp_raw_rhn_correction$lon+MMPRNT_16s.metadata_mapp_raw_rhn_correction$xCord+
  (MMPRNT_16s.metadata_mapp_raw_rhn_correction$abs_yCord*MMPRNT_16s.metadata_mapp_raw_rhn_correction$coeff_lon_m)
MMPRNT_16s.metadata_mapp_raw_rhn_correction$UTM_Lat_Cord=MMPRNT_16s.metadata_mapp_raw_rhn_correction$lat+MMPRNT_16s.metadata_mapp_raw_rhn_correction$abs_yCord+
  (MMPRNT_16s.metadata_mapp_raw_rhn_correction$xCord*MMPRNT_16s.metadata_mapp_raw_rhn_correction$coeff_lat_m)


ggplot(MMPRNT_16s.metadata_mapp_raw_rhn_correction,aes(x=UTM_Lon_Cord,y=UTM_Lat_Cord))+geom_point()+geom_point(data = MMPRNT_16s.metadata_mapp_raw_rhn_correction,aes(x=lon,y=lat),color="red",size=4)
nrow(MMPRNT_16s.metadata_mapp_raw_rhn_correction)
#72

#Add in the gps points for Hancoc 

MMPRNT_16s.metadata_mapp_raw_han=subset(MMPRNT_16s.metadata_mapp_raw, siteID=="HAN")

MMPRNT_16s.metadata_mapp_raw_han_correction=merge(MMPRNT_16s.metadata_mapp_raw_han,MLE_coord2_m2_HAN_G5_UTM_correct2, by.x="treatment_GLBRC",by.y="treatment" )

MMPRNT_16s.metadata_mapp_raw_han_correction$UTM_Lon_Cord=MMPRNT_16s.metadata_mapp_raw_han_correction$lon+MMPRNT_16s.metadata_mapp_raw_han_correction$xCord+
  (MMPRNT_16s.metadata_mapp_raw_han_correction$abs_yCord*MMPRNT_16s.metadata_mapp_raw_han_correction$coeff_lon_m)
MMPRNT_16s.metadata_mapp_raw_han_correction$UTM_Lat_Cord=MMPRNT_16s.metadata_mapp_raw_han_correction$lat+MMPRNT_16s.metadata_mapp_raw_han_correction$abs_yCord+
  (MMPRNT_16s.metadata_mapp_raw_han_correction$xCord*MMPRNT_16s.metadata_mapp_raw_han_correction$coeff_lat_m)


ggplot(MMPRNT_16s.metadata_mapp_raw_han_correction,aes(x=UTM_Lon_Cord,y=UTM_Lat_Cord))+geom_point()+geom_point(data = MMPRNT_16s.metadata_mapp_raw_han_correction,aes(x=lon,y=lat),color="red",size=4)

MMPRNT_16s.metadata_mapp_raw_han_correction[,c("X","avg_corr_lon","avg_corr_lat","lon","lat","site","coeff_lon_m","coeff_lat_m")]=NULL
MMPRNT_16s.metadata_mapp_raw_rhn_correction[,c("X","avg_corr_lon","avg_corr_lat","lon","lat","site","coeff_lon_m","coeff_lat_m")]=NULL
MMPRNT_16s.metadata_mapp_raw2=rbind(subset(MMPRNT_16s.metadata_mapp_raw,siteID!="RHN"&siteID!="HAN"),MMPRNT_16s.metadata_mapp_raw_han_correction,MMPRNT_16s.metadata_mapp_raw_rhn_correction)
nrow(MMPRNT_16s.metadata_mapp_raw2)
#2430

MMPRNT.metadata_mapp_trunc=MMPRNT_16s.metadata_mapp_raw2[,c("sampleID_long","collectionDate","siteID","plotType","FertStatus","plotRep","pseudo_no","year","dis_Mid","UTM_Lat_Cord","UTM_Lon_Cord")]

MMPRNT.metadata_mapp_trunc[1:10,1:10]
nrow(MMPRNT.metadata_mapp_trunc)
#2430
#colnames(GLBRC018_bact_map)[colnames(GLBRC018_bact_map)=="field_sampleID"]="sampleID_long"

GLBRC018_bact_map_metadata=merge(GLBRC018_bact_map,MMPRNT.metadata_mapp_trunc, by="sampleID_long", all.x = T)
head(GLBRC018_bact_map_metadata)
summary(GLBRC018_bact_map_metadata)
colnames(GLBRC018_bact_map_metadata)
row.names(GLBRC018_bact_map_metadata)=GLBRC018_bact_map_metadata$SampleID 
GLBRC018_bact_map_metadata[1:10,1:10]
unique(GLBRC018_bact_map_metadata$Project)

write.csv(subset(GLBRC018_bact_map_metadata, Root_soil=="Root"&Project=="MMPRNT"),here::here("Bact_HPCC_out","NCBI_submission","Root_Raw_mapping.csv"))


#OTU table
otu_GLBRC018_bact=otu_table(read.table(here::here("Bact_HPCC_out","combined_soil_roots_merged_16S_GLBRC018_OTU_table.txt"),sep = "\t", header = T, row.names = 1), taxa_are_rows = T)
otu_GLBRC018_bact[1:10,1:10]


#Load Taxa table

#I am going to use the sintax classification for now 


GLBRC018_bact_raw_SINTAX_SILVA123 = read.delim(here::here("Bact_HPCC_out","combined_soil_roots_merged_16S_GLBRC018_otus_silva123_taxonomy.SINTAX"),
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


#Let's choose the best samples from the Re PCR of poorly amplifying samples 
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


#Pre-filitering abundances
GLBRC018_OTU_bact_MMPRNT_fin=subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_fin,seq_samp_type!="mock")
GLBRC018_OTU_bact_MMPRNT_fin=prune_taxa(taxa_sums(GLBRC018_OTU_bact_MMPRNT_fin) > 0, GLBRC018_OTU_bact_MMPRNT_fin)
ntaxa(GLBRC018_OTU_bact_MMPRNT_fin)
# 46923
sum(taxa_sums(GLBRC018_OTU_bact_MMPRNT_fin))
#31944520

mean(sample_sums(GLBRC018_OTU_bact_MMPRNT_fin))
#28935.25

min(sample_sums(GLBRC018_OTU_bact_MMPRNT_fin))
#162
max(sample_sums(GLBRC018_OTU_bact_MMPRNT_fin))
#112729

nsamples(GLBRC018_OTU_bact_MMPRNT_fin)
#1104

#Roots 
sum(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Root")))
#11536151
nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Root"))
#234
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Root")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Root")))
#27005

#Soils 
sum(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Soil")))
#20408369
nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Soil"))
#870
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Soil")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_fin,Root_soil=="Soil")))
#45633

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

#####Remove non bacterial reads####
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

#save(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar, file = here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))
load(here::here("R_files","GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData"))


write.table(data.frame(otu_table(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)),
            here::here("R_files","Publish_data","OTU_GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.txt"),
            sep = "\t")

write.csv(data.frame(tax_table(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)),
          here::here("R_files","Publish_data","Tax_table_GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.csv"))

write.csv(data.frame(sample_data(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar)),
          here::here("R_files","Publish_data","metadata_GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.csv"))


ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,seq_samp_type!="mock")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,seq_samp_type!="mock")))
#31837
sum(otu_table(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,seq_samp_type!="mock")))
# 10920000
nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,seq_samp_type!="mock"))
#1092

sort(sample_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,seq_samp_type!="mock")))[1:20]

#roots
nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"))
#233
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root")))
#17499
sum(taxa_sums(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Root"))))
#2330000

#Soils
nsamples(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Soil"))
#859
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Soil")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Soil")))
#31361
sum(taxa_sums(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Soil")) > 0, 
                 subset_samples(GLBRC018_OTU_bact_MMPRNT_mock_bact_rar,Root_soil=="Soil"))))
#8590000

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
MMPRNT_16s.metadata_mapp_raw=read.table(here::here("D:/MMPRNT_16S_016-018/ABIOTIC_data/MMPRNT.metadata_mapp_raw_comp_OUTLIER_filter.txt"))
nrow(MMPRNT_16s.metadata_mapp_raw)
MMPRNT_16s.metadata_mapp_raw$treatment_GLBRC=paste(MMPRNT_16s.metadata_mapp_raw$plotType,MMPRNT_16s.metadata_mapp_raw$plotRep,sep = "R")

#Add in the missing gps ponts

MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correct2=read.csv(here::here("ABIOTIC_data","MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correction_coef.csv"))
MLE_coord2_m2_RHN_G5_UTM_correct2=subset(MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correct2,site=="Rhinelander")
MLE_coord2_m2_HAN_G5_UTM_correct2=subset(MLE_coord2_m2_ESC_RHN_HAN_G5_UTM_correct2,site=="Hancock")


#Add in the gps points for Rhinelander 

MMPRNT_16s.metadata_mapp_raw_rhn=subset(MMPRNT_16s.metadata_mapp_raw, siteID=="RHN")
nrow(MMPRNT_16s.metadata_mapp_raw_rhn)
#216
MMPRNT_16s.metadata_mapp_raw_rhn_correction=merge(MMPRNT_16s.metadata_mapp_raw_rhn,MLE_coord2_m2_RHN_G5_UTM_correct2, by.x="treatment_GLBRC",by.y="treatment" )

MMPRNT_16s.metadata_mapp_raw_rhn_correction$UTM_Lon_Cord=MMPRNT_16s.metadata_mapp_raw_rhn_correction$lon+MMPRNT_16s.metadata_mapp_raw_rhn_correction$xCord+
  (MMPRNT_16s.metadata_mapp_raw_rhn_correction$abs_yCord*MMPRNT_16s.metadata_mapp_raw_rhn_correction$coeff_lon_m)
MMPRNT_16s.metadata_mapp_raw_rhn_correction$UTM_Lat_Cord=MMPRNT_16s.metadata_mapp_raw_rhn_correction$lat+MMPRNT_16s.metadata_mapp_raw_rhn_correction$abs_yCord+
  (MMPRNT_16s.metadata_mapp_raw_rhn_correction$xCord*MMPRNT_16s.metadata_mapp_raw_rhn_correction$coeff_lat_m)


ggplot(MMPRNT_16s.metadata_mapp_raw_rhn_correction,aes(x=UTM_Lon_Cord,y=UTM_Lat_Cord))+geom_point()+geom_point(data = MMPRNT_16s.metadata_mapp_raw_rhn_correction,aes(x=lon,y=lat),color="red",size=4)
nrow(MMPRNT_16s.metadata_mapp_raw_rhn_correction)
#72

#Add in the gps points for Hancoc 

MMPRNT_16s.metadata_mapp_raw_han=subset(MMPRNT_16s.metadata_mapp_raw, siteID=="HAN")

MMPRNT_16s.metadata_mapp_raw_han_correction=merge(MMPRNT_16s.metadata_mapp_raw_han,MLE_coord2_m2_HAN_G5_UTM_correct2, by.x="treatment_GLBRC",by.y="treatment" )

MMPRNT_16s.metadata_mapp_raw_han_correction$UTM_Lon_Cord=MMPRNT_16s.metadata_mapp_raw_han_correction$lon+MMPRNT_16s.metadata_mapp_raw_han_correction$xCord+
  (MMPRNT_16s.metadata_mapp_raw_han_correction$abs_yCord*MMPRNT_16s.metadata_mapp_raw_han_correction$coeff_lon_m)
MMPRNT_16s.metadata_mapp_raw_han_correction$UTM_Lat_Cord=MMPRNT_16s.metadata_mapp_raw_han_correction$lat+MMPRNT_16s.metadata_mapp_raw_han_correction$abs_yCord+
  (MMPRNT_16s.metadata_mapp_raw_han_correction$xCord*MMPRNT_16s.metadata_mapp_raw_han_correction$coeff_lat_m)


ggplot(MMPRNT_16s.metadata_mapp_raw_han_correction,aes(x=UTM_Lon_Cord,y=UTM_Lat_Cord))+geom_point()+geom_point(data = MMPRNT_16s.metadata_mapp_raw_han_correction,aes(x=lon,y=lat),color="red",size=4)

MMPRNT_16s.metadata_mapp_raw_han_correction[,c("X","avg_corr_lon","avg_corr_lat","lon","lat","site","coeff_lon_m","coeff_lat_m")]=NULL
MMPRNT_16s.metadata_mapp_raw_rhn_correction[,c("X","avg_corr_lon","avg_corr_lat","lon","lat","site","coeff_lon_m","coeff_lat_m")]=NULL
MMPRNT_16s.metadata_mapp_raw2=rbind(subset(MMPRNT_16s.metadata_mapp_raw,siteID!="RHN"&siteID!="HAN"),MMPRNT_16s.metadata_mapp_raw_han_correction,MMPRNT_16s.metadata_mapp_raw_rhn_correction)
nrow(MMPRNT_16s.metadata_mapp_raw2)
#2430

MMPRNT.metadata_mapp_trunc=MMPRNT_16s.metadata_mapp_raw2[,c("sampleID_long","collectionDate","siteID","plotType","FertStatus","plotRep","year","dis_Mid","UTM_Lat_Cord","UTM_Lon_Cord","pseudo_no")]

MMPRNT.metadata_mapp_trunc[1:10,1:10]
nrow(MMPRNT.metadata_mapp_trunc)
#2430
colnames(GLBRC018_fung_map)[colnames(GLBRC018_fung_map)=="field_sampleID"]="sampleID_long"

GLBRC018_fung_map_metadata=merge(GLBRC018_fung_map,MMPRNT.metadata_mapp_trunc, by="sampleID_long", all.x = T)
head(GLBRC018_fung_map_metadata)
summary(GLBRC018_fung_map_metadata)
colnames(GLBRC018_fung_map_metadata)
row.names(GLBRC018_fung_map_metadata)=GLBRC018_fung_map_metadata$sampleID_seq 
GLBRC018_fung_map_metadata[1:10,1:10]

unique(GLBRC018_fung_map_metadata$project)
#NCBI SUMBMISSION FILES 

write.csv(subset(GLBRC018_fung_map_metadata, Root_soil=="Soil"&project=="Field2018"),here::here("Fung_HPCC_out","NCBI_submission","Soil_Raw_mapping.csv"))
write.csv(subset(GLBRC018_fung_map_metadata, Root_soil=="Root"&project=="Field2018"),here::here("Fung_HPCC_out","NCBI_submission","Root_Raw_mapping.csv"))


#OTU table
otu_GLBRC018_fung=otu_table(read.table(here::here("Fung_HPCC_out","combined_merged_ITS1_2_GLBRC_OTU_table.txt"),sep = "\t", header = T, row.names = 1), taxa_are_rows = T)
otu_GLBRC018_fung[1:10,1:10]


#Load Taxa table


#I am going to use the sintax classification for now 


GLBRC018_fung_raw_SINTAX_UNITE8.2 = read.delim( here::here("Fung_HPCC_out","taxonomy_combined_merged_ITS1_2_GLBRC_otus.SINTAX"),
                                                sep = "\t",header = F,fill=T)

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


#Let's choose the the best samples from the Re PCR of poorly amplifying samples
unique(sample_data(GLBRC018_OTU_fung_MMPRNT_mock)$rep_pcr) 


GLBRC018_OTU_fung_MMPRNT_mock_re_PCR=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock,rep_pcr!="none")
sample_data(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)
nsamples(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)
#76

unique(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR)$sampleID_long) 


GLBRC018_OTU_fung_MMPRNT_mock_re_PCR_diversity=merge(estimate_richness(GLBRC018_OTU_fung_MMPRNT_mock_re_PCR,
                                                                       measures = c("Observed", "InvSimpson", "Shannon", "Chao1")),
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


#Pre-filitering abundances
GLBRC018_OTU_fung_MMPRNT_fin=subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fin,seq_samp_type!="mock")
GLBRC018_OTU_fung_MMPRNT_fin=prune_taxa(taxa_sums(GLBRC018_OTU_fung_MMPRNT_fin) > 0, GLBRC018_OTU_fung_MMPRNT_fin)
ntaxa(GLBRC018_OTU_fung_MMPRNT_fin)
# 6871
sum(taxa_sums(GLBRC018_OTU_fung_MMPRNT_fin))
#105206210

mean(sample_sums(GLBRC018_OTU_fung_MMPRNT_fin))
#95729.04

min(sample_sums(GLBRC018_OTU_fung_MMPRNT_fin))
#1
max(sample_sums(GLBRC018_OTU_fung_MMPRNT_fin))
#342074

nsamples(GLBRC018_OTU_fung_MMPRNT_fin)
#1099

#Roots 
sum(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Root")))
#14818135
nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Root"))
#234
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Root")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Root")))
#5963

#Soils 
sum(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Soil")))
#90388075
nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Soil"))
#865
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Soil")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_fin,Root_soil=="Soil")))
#6858
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

write.table(data.frame(otu_table(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)),
            here::here("R_files","Publish_data","OTU_GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.txt"),
            sep = "\t")

write.csv(data.frame(tax_table(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)),
          here::here("R_files","Publish_data","Tax_table_GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.csv"))

write.csv(data.frame(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)),
          here::here("R_files","Publish_data","metadata_GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.csv"))



ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,seq_samp_type!="mock")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,seq_samp_type!="mock")))
#3088
sum(otu_table(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,seq_samp_type!="mock")))
#10740000
nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,seq_samp_type!="mock"))
#1074

unique(sample_data(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar)$Root_soil)

sort(sample_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,seq_samp_type!="mock")))[1:20]

#roots
nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&
                          seq_samp_type!="mock"))
#228
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&
                                            seq_samp_type!="mock")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&
                                  seq_samp_type!="mock")))
#2604

sum(taxa_sums(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&
                                            seq_samp_type!="mock")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Root"&
                                  seq_samp_type!="mock"))))

#Soils
nsamples(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Soil"&
                          seq_samp_type!="mock"))
#846
ntaxa(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Soil"&
                                            seq_samp_type!="mock")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Soil"&
                                  seq_samp_type!="mock")))
#3076

sum(taxa_sums(prune_taxa(taxa_sums(subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Soil"&
                                            seq_samp_type!="mock")) > 0, 
                 subset_samples(GLBRC018_OTU_fung_MMPRNT_mock_fung_rar,Root_soil=="Soil"&
                                  seq_samp_type!="mock"))))
#8460000
#####END: Pre-processing Fungal community####