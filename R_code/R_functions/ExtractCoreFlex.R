#Modified from https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R
# ExtractCore function ------------------------------------------------------------------------------------
# Modofied form Stopnisek and Shade 2019 - Current Opinion in Microbiology
# https://www.sciencedirect.com/science/article/pii/S1369527419300426?via%3Dihub


# Extract Core OTUs function ------------------------------------------------------------------------------

# Example: ExtractCore(phyloseq_object, "Date", "elbow", Treatment", "Fertilized", "")

# This function will filter the dataset by Fertilized samples only and calculate
# the core according to the different sampling dates.
#method is "total", "elbow", or number which sets the Bray-Curtis cut of a percent (whole number)

ExtractCoreFlex <- function(physeq, Var, method, Group=NULL, Level=NULL){
  require(magrittr)
  require(vegan)
  require(dplyr)
  require(tidyverse)
  library(lazyeval)
  #set.seed(100)
  # input dataset needs to be rarified and minimum depth included
  if (min(sample_sums(physeq)) == max(sample_sums(physeq))) {
    print("Dataset is rarefied at a depth of:")
    nReads = min(sample_sums(physeq))
    nReads %T>% print()
    rare <- physeq
    rare %T>% print()
    taxon <- as(tax_table(rare), "matrix")
    taxon <- as.data.frame(taxon)
    dim(taxon)  %T>% print()
  } else {
    print("Performing Rarefaction. Minimum depth:")
    nReads = min(sample_sums(physeq))
    nReads %T>% print()
    rare <-
      rarefy_even_depth(
        physeq, sample.size = nReads,
        trimOTUs = TRUE,
        rngseed = TRUE,
        replace = TRUE)
    taxon <- as(tax_table(rare), "matrix")
    taxon <- as.data.frame(taxon)
  }
  # choosing a subset or using the whole phyloseq object as is
  if (is.null(Group)) {
    print("Dataset not subsetted")
    otu <- rare@otu_table %>% as("matrix")
    map <- rare@sam_data %>% as("data.frame")
  } else{
    print("Dataset subsetted using the Group variable")
    sub_group <- substitute(Group)
    sub_set <- subset(sample_data(rare), eval(parse(text=sub_group)) %in% Level)
    physeq1 <- merge_phyloseq(otu_table(rare),
                              tax_table(rare),
                              refseq(rare),
                              sub_set)
    otu_table(physeq1) <- otu_table(physeq1)[which(rowSums(otu_table(physeq1)) > 0), ]
    otu <- physeq1@otu_table %>% as("matrix")
    map <- physeq1@sam_data %>% as("data.frame")
    print("Grouping Factor")
    map[,Group] %T>% print()
    taxon <- as(tax_table(physeq1), "matrix")
    taxon <- as.data.frame(taxon)
  }
  map$SampleID <- rownames(map)
  print("Check: dimension of datframe and metadata")
  dim(otu) %T>% print() # funcitons form magrittr package
  dim(map) %T>% print() 
  # calculating occupancy and abundance
  otu_PA <-
    1 * ((otu > 0) == 1)                                             # presence-absence data
  otu_occ <-
    rowSums(otu_PA) / ncol(otu_PA)                                   # occupancy calculation
  otu_rel <-
    apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # mean relative abundance
  occ_abun <-
    tibble::rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)), "otu")     # combining occupancy and abundance data frame
  # NOTE! tibble::rownames_to_column is deprecated and generates a warning, a bug of tidyverse, 
  # alternative you can use:
  # occ_abun <- tibble::rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),"otu")
  # Ranking OTUs based on their occupancy
  # For caluclating raking index we included following conditions:
  # - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
  # - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
  Var <- enquo(Var) # lazy evaluation
  PresenceSum <-
    data.frame(otu = as.factor(row.names(otu)), otu) %>%
    gather(SampleID, abun,-otu) %>%
    left_join(map, by = "SampleID") %>%
    group_by(otu, !!Var) %>%
    dplyr::summarise(
      time_freq = sum(abun > 0) / length(abun), # frequency of detection between time points
      coreTime = ifelse(time_freq == 1, 1, 0)) %>% # 1 only if occupancy 1 with specific time, 0 if not
    group_by(otu) %>%
    dplyr::summarise(
      sumF = sum(time_freq),
      sumG = sum(coreTime),
      nS = length(!!Var)* 2,
      Index = (sumF + sumG) / nS) # calculating weighting Index based on number of points detected
  # PresenceSum %T>% print()
  # ranking otus
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by = "otu") %>%
    transmute(otu = otu,
              rank = Index) %>%
    arrange(desc(rank))
  #otu_ranked %T>% print()
  # calculating BC dissimilarity based on the 1st ranked OTU
  otu_start = otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start, ])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x)
    sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
  x_names <-
    apply(combn(ncol(start_matrix), 2), 2, function(x)
      paste(colnames(start_matrix)[x], collapse = "-"))
  # creating a data.frame and adding first OTU name as 1
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  # Initialize your data structures:  calculating the contribution of ranked OTUs to the BC similarity
  BCaddition <- NULL
  BCaddition <- rbind(BCaddition,df_s)
  # calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. 
  # Can be set to the entire length of OTUs in the dataset.
  # it might take some time if more than 5000 OTUs are included.
  for(i in 2:1000){ #nrow(otu_ranked)
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    y <-
      apply(combn(ncol(start_matrix), 2), 2, function(y)
        sum(abs(start_matrix[, y[1]] - start_matrix[, y[2]])) / (2 * nReads))
    df_a <- data.frame(x_names, y)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  # Calculating the BC dissimilarity of the whole dataset (not needed if the second loop 
  # is already including all OTUs)
  z <-
    apply(combn(ncol(otu), 2), 2, function(z)
      sum(abs(otu[, z[1]] - otu[, z[2]])) / (2 * nReads))
  # overwrite the names here
  x_names <-
    apply(combn(ncol(otu), 2), 2, function(x)
      paste(colnames(otu)[x], collapse = "-"))
  df_full <- data.frame(x_names, z)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition, df_full, by='x_names')
  # ranking the obtained BC 
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    # Calculate mean Bray-Curtis dissimilarity
    summarise(MeanBC=mean(BC)) %>%            
    arrange(desc(-MeanBC)) %>%
    # Calculate proportion of the dissimilarity explained by the n number of ranked OTUs 
    mutate(proportionBC=MeanBC/max(MeanBC))
  #BC_ranked %T>% print()
  # Calculating the increase BC
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  if (method=="elbow"){
    fo_difference <- function(pos){
      left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
      right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
      return(left - right)
    }
    BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
    elbow <- which.max(BC_ranked$fo_diffs)
    core_otus <- otu_ranked$otu[1:elbow]
    core_otus %T>% print()
  }
  # Creating threshold for core inclusion - last call method using a
  # final increase in BC similarity of equal or greater than number set by method
  else{
    lastCall <-
      as.numeric(as.character(dplyr::last(
        subset(BC_ranked, IncreaseBC >= 1+(method*0.01))$rank)))
    # lastCall <- last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))
    core_otus <- otu_ranked$otu[1:lastCall]
    core_otus %T>% print()
  }
  # Adding Core otus for reating occupancy abundance plot
  occ_abun$fill <- 'no'
  occ_abun$fill[occ_abun$otu %in% core_otus] <- "core"
  return_list <-
    list(core_otus, BC_ranked, otu_ranked, occ_abun, otu, map, taxon)
  return(return_list)
}
