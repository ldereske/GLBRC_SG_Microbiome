# Functions called by the R code in the R_code folder of this directory

### ExtractCoreFlex.R                    
Code for the identification and extraction of core taxa within an OTU table. Original code is from [Gian Maria Niccolò Benucci]{https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R}. I modified the code to allow for the ease of setting the Bray-Curtis threshold. 

### PlotBCincrease_Richness_Abundance.R                             
Code for the graphing of the richness of OTUs that were included in the core community at a given Bray-Curtis threshold. The PlotBCincrease() function is from [Gian Maria Niccolò Benucci]{https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R}


### gdm.varImp_MOD.R
Code for calculating the importance of factors in a gdm. This is a modified form of [gdm.varImp]{https://www.rdocumentation.org/packages/gdm/versions/1.5.0-3/topics/gdm.varImp} from the R package [gdm]{https://www.rdocumentation.org/packages/gdm/versions/1.5.0-3}. I modified the code to allow for the assessment of models with 3 or less factors. This code does not allow for the model selection that was available in the original function.

### FitNeutral.R    
Code for fitting and assessment of neutral models defined in [Burns et al. - 2016 ISME JOURNAL]{https://www.nature.com/articles/ismej2015142} with code modified and implemented by [Gian Maria Niccolò Benucci]{https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R}

### PlotNeutral.R                        
Code for graphing of neutral models defined in [Burns et al. - 2016 ISME JOURNAL]{https://www.nature.com/articles/ismej2015142} with code modified and implemented by [Gian Maria Niccolò Benucci]{https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R}

### dataUncertainty.R 
Code for the calculation of uncertainty intervals from [gdm]{https://www.rdocumentation.org/packages/gdm/versions/1.5.0-3} models. This modified from the original function [plotUncertainty]{https://rdrr.io/cran/gdm/man/plotUncertainty.html} to allow for the use of this model uncertainy in other plotting frameworks.

### sncm.fit.R
Code for fitting and assessment of neutral models defined in [Burns et al. - 2016 ISME JOURNAL]{https://www.nature.com/articles/ismej2015142} with code modified and implemented by [Gian Maria Niccolò Benucci]{https://github.com/Gian77/Scientific-Papers-R-Code/blob/master/VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome/ExtractCore.R}
