# NutNetLitterDecomp
This project organizes data files and code associated with a seven year oak (*Quercus ellipsoidalis*) leaf litter decomposition across Nutrient Network sites in North America, Australia, and Europe. Information about the Nutrient Network experiment can be found [here](https://nutnet.org/). 

## Methods associated with litter decomposition experiment
Freshly fallen leaf litter (*Quercus ellipsoidalis*) was collected at Cedar Creek Ecosystem Science Researve (45.4020 °N, 93.1994 °W) in October 2008. Litter bags were constructed of 1-mm mesh fiberglass window screen and contained approximately 10 g (dry weight) of autoclave-steriled leaf litter. Bags were strung together in groups of seven to allow for seven annual harvests and sent to individual Nutrient Network site investigators. Site investigators deployed litterbags by pinning a single string to the ground in each experimental plot between December 2009 and October 2010. Here we evaluate the effect of factorial fertilization treatments on litter decomposition dynamics.

One litter bag was collected rom each plot at approximately annual intervals from 2010 to 2016. Harvested bags were cleaned of contaminating material, dried at 65°C to constant mass, and sent to the University of Minnesota for further processing. Litter was ground and analyzed for total carbon (C) and nitrogen (N) on an elemental analyzer (ESC 4010, Costech Analytical Technologies, Valencia, CA, USA), and a subsample was ashed (600°C for 6 hours) to determine ash-free dry mass (AFDM). Proportion initial mass remaining was converted to proportion initial C remaining for further analyses to account for contamination of litter by soil. Outliers were handled as follows. Quantile plots were used to screen for possible outlier values that were then inspected individually. Values were excluded if investigators noted signs of major damage to litter bags that likely caused loss of material from bag. In addition, for outliers caused by anomalously high or low values of %C, we adjusted their %C by predicting it from %AFDM using the regression between %C and %AFDM among all samples (%C = 0.75 + 0.50*%AFDM, R2 = 0.77, P<0.0001). When there was insufficient mass remaining to analyze for %C, we used the average %C from other samples harvested at the same time in the same treatment. Sites differed in their frequency of litter bag harvests, and some sites were unable to complete the study. In our anaalysis, we summarize data across 20 sites that contibuted data from at least three and up to seven annual harvests. 

The code and datasets provided here allow you to:
1. **Fit Decomposition Models by Site x Treatment combination**. This code allows you to fit four decomposition models (single exponential, double exponential, asymptotic exponential, and Weibull) to C mass loss data pooled across site x treatment combinations. We used this approach to identify the best-fit models across sites and treatments. Follwing the site x treatment model fitting, we excluded the double exponential model from all subsequent analyses, as it best-described the decomposition series in a small minority of site x treatments. 

 - **Data**: NutNet_OakLitterDecay_HarvestData.csv
 - **Metadata**: NutNet_OakLitter_Decomp_Metadata.csv
 - **Code**: DecompModelFitting_SiteTrt.R
 
 Single exponential, asymptotic exponential, and Weibull decomposition models were then fit to decomposition series within individual plots.
 - **Code**: DecompModelFitting_Plot.R
 
2. **Evaluate Treatment Responses**. Evaluate effect of fertilization treatments on model decomposition model parameters. 
 - **Data**: 
 - **Metadata**: 
 - **Code**: NutNetOakLitter_TreatmentResponses.R
 
The first line of each code file needs to be updated to match the data directory on the user's computer.  This project was built under R version3.6.1 and depends on the following packages: stats4, bbmle, readr, plyr, dplyr, Matrix, Hmisc, gglplot2, nlme, multcomp, MuMIn, boot.