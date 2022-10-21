# FitNeutral function ------------------------------------------------------------

# Example: FitNeutral(core_otus) 

# Fitting neutral models through the sncm.fit.R function modified form 
# Burns et al. - 2016 ISME JOURNAL

FitNeutral <- function(otu_core) {
  source(here::here("R_code","R_functions","sncm.fit.R"))
  #generating df
  taxon <- otu_core[[7]]
  spp <- t(otu_core[[5]])
  occ_abun <- otu_core[[4]]
  names(occ_abun)[names(occ_abun) == "otu"] <- "OTU_ID"
  occ_abun %T>% print()
  # source community pool
  meta <- otu_core[[6]]
  # using root as the pool gave weird results
  # soil_pool <- spp[ rownames(subset(meta, Origin=="Soil")), ]
  # root_pool <- spp[ rownames(subset(meta, Origin=="Root")), ]
  #fitting model
  obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)
  obs.np$OTU_ID=row.names(obs.np)
  sta.np <- sncm.fit(spp, taxon, stats = TRUE, pool = NULL)
  #
  obs.np$fit_class <- "As predicted"
  obs.np[which(obs.np$freq < obs.np$pred.lwr), "fit_class"] <-
    "Below prediction"
  obs.np[which(obs.np$freq > obs.np$pred.upr), "fit_class"] <-
    "Above prediction"
  obs.np[which(is.na(obs.np$freq)), "fit_class"] <- "NA"
  as.data.frame(left_join(occ_abun, obs.np, by = 'OTU_ID')) -> fit_table
  #
  sta.np$above.pred <-
    sum(obs.np$freq > (obs.np$pred.upr), na.rm = TRUE) / sta.np$Richness
  sta.np$below.pred <-
    sum(obs.np$freq < (obs.np$pred.lwr), na.rm = TRUE) / sta.np$Richness
  fit_res <- as.data.frame(sta.np)
  rownames(fit_res) <- "Value"
  fit_res
  list_tab <- list(fit_res, fit_table)
  return(list_tab)
}