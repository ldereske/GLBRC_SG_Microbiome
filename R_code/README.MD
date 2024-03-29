# Code for the analysis of the switchgrass microbiome in the R environment

### OTU_bacterial_fungal_pre_processing.R
Code for the pre processing and assessment of microbial communities post bioinformatic processing detailed in the Bioinformatic_code directory of this repository. The files are not provided for this processing due the inclusion of samples that are not published. 

Please contact [Lukas Bell-Dereske](lukas.dereske@gmail.com) for the raw data or find the raw sequences in the NCBI Short Read Archive under BioProject PRJNA733764 for bacterial communities and PRJNA799201 for fungal communities.

This code was used to produce 
+ Tables S3 - S5

### NMDS_GDM_MLE_graph_PERMANOVA.R and NMDS_GDM_Lux_Arbor_graph_PERMANOVA.R  
The analyses and plotting of community composition, along with abiotic and biotic drivers of community composition, from both bacterial and fungal communities. 
This code was used to produce 
+ Tables 1 & 2 with the PERMANOVA analyzed using [PRIMER v6](https://www.researchgate.net/publication/235425881_Primer_v6_User_ManualTutorial)
+ Table 3 
+ Tables S1, S2, & S11
+ Figures 1-3, S10, & S13 - S15


### Diversity_MLE_graph_analyses.R and Diversity_Lux_Arbor_graph_analyses.R 
The analyses and plotting of community diversity from both bacterial and fungal communities. 

This code was used to produce 
+ Tables S8 & S10
+ Figures S10 & S12

### Similarity_soilDate_root_Lux_Arbor_graph_analyses.R 
The analyses of the similarity of root and soil communities based on the temporal sampling of the two niches.

This code was used to produce 
+ Table S12
+ Figure 4 
### Procrustes_MLE_graph_residual_analyses.R and Procrutes_Lux_Arbor_graph_residual_analyses.R
The analyses of the concordance between bacterial and fungal communities in roots and soils.

This code was used to produce 
+ Tables S7 & S9
+ Figures S9 & S11


### FunGuild_database_creation.R                        
Code for the creation and processing of the [FunGuild](https://github.com/UMNFuN/FUNGuild) database used for the charaterization of the traits in the core fungal community.

This code was used to produce 
+ Tables S6
                            
### Core_Taxa_metacoder_MLE_graph_analyses.R and Core_Taxa_metacoder_Lux_Arbor_graph_analyses.R 
Code for the charaterization and plotting of the core communities.

This code was used to produce 
+ Table S14
+ Figures 5 & 6
+ Figures S1 - S8, & S17 - S20

### Betadisp_Lux_Arbor_graph_analyses.R 
The analysis of the betadispersion of the microbial communities of Lux Arbor.

This code was used to produce 
+ Table S13
+ Figure S16