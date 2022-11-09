# Code for the analysis of spatiotemporal dynamics in switchgrass root and soil microbiomes

This is the code for a manuscript that is currently under review. The project explored the spatiotemporal structure of the switchgrass microbiome of the [GLBRC Marginal Lands Experiment](https://lter.kbs.msu.edu/research/long-term-experiments/marginal-land-experiment/). 


## Data necessary for the analyses can be found at [DOI 10.17605/OSF.IO/5VW9C](https://osf.io/5vw9c/)
You can download the data into your R environment by running the below code

```
remotes::install_github("ropensci/osfr")
library(osfr)
library(here)

#Create the directory
dir.create("Publish_data")

#Retrieve the files
osf_retrieve_file("63513ed993d3520c4891c7da") |>
  osf_download(path = here::here("Publish_data"))
  
  
```

## Directories in the repository

### R_code
Code necessary for the analysis of the microbiome data. 

### R_functions
Functions called by the code in  R_code. Code is modified from R package [gdm](https://www.rdocumentation.org/packages/gdm/versions/1.5.0-3) and code originally written by [Gian Maria Niccolò Benucci](https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R). Please cite these sources. 

### Bioinformatic_code
Code used to process the raw bacterial and fungal reads from MiSeq sequencing. We used [usearch](https://www.drive5.com/usearch/), [QIIME2](https://qiime2.org/), and [CONSTAXv2](https://constax.readthedocs.io/en/latest/) within our modified pipeline. The raw reads can be found in the NCBI Short Read Archive under BioProject PRJNA733764 for bacterial communities and PRJNA799201 for fungal communities.

