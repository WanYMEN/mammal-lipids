# mammal-lipids

This depository contains code to generate the main results and plots in the paper that analyzes coevolution between mammalian brain and milk.

The scripts include:

   data_normalize.R - normalizes the raw intensity data matrixes of fatty acids (FAs) based on quantiles.
   
   milk_FA.mds.R - performs multidimensional scaling (MDS) for milk samples.
   
   brain_FA.mds.age.R - performs MDS for brain samples and plots the results with age being considered.
   
   milk_FA.percentage_bar.species.R - calculates concentration percentages of milk samples across species. Similar ways were used to calculate concentration percentages of milk samples across groups, primates and two human populations.
   
   brain_FA.percentage_bar.R - calculates concentration percentages of brain samples across species and primates.

   milk_FA.specific.R - calculate FAs specific to species, groups and primates based on milk samples.

   milk_FA.LS22_40.lmAdjustedNoPopulation.R - calculate FAs specific to human populations based on milk samples.

   brain_FA.specific.R - calculate FAs specific to species, groups and primates based on brain samples for each brain region.

   brain_FA.factorAdjusted.R - adjust impacts of age and sex before concentration percentage calculation.

   brain_FA.age_slope.primate.R - calculate the slope of an FA’s intensity changes with age across primate brains.

   brain_milk_FA.species_FC.R - calculate fold changes (FCs) of intensity between two species for milk and brain, respectively.

   brain_milk_FA.species_FC_LS_correlation.R - calculate correlation of intensity FCs of human to other species between milk and brain for each group of lactation stages.

   brain_milk_FA.speciesIntensityCorrelation.R - calculate correlation of intensity between milk and brain for a same species.
