# CausalChallenge2016
Code and notes for the "Causal Challenge" - July 2016

## Meeting June 22, 2016
- Gabrielle: Try pooling predicted values from individual MICE-imputed datasets to get a final set of predicted values.  Test this by deleting observations from another variable, imputing, and comparing predicted values to true values.  For now, do without litter variable, but will introduce later when we figure out how to handle it.
- Kevin:  Look into machine learning methods for clustering data.  The goal is to reduce the 190 litters into a variable with a much smaller number of categories by clustering similar litters together based on phenotype.
- Will hold off on causal part for now, since we need to get good predicted values before considering the problems in causal inference.

## Meeting June 30, 2016
- Gabrielle: Run simulations with MICE (simulations as described above), determine whether linear regression vs. predictive mean matching performs better. Eventually include litter effect variable(s).
- Kevin: Construct variable(s) summarizing the litter effect via PC. Idea: group litters using the median of the 1st PC per litter. Including the first two PC directly in the prediction model did not work because of collinearity. 
- For next week's meeting: look into the causal interpretation part. Start thinking about the poster.

## Meeting July 8, 2016
Prediction of missing data 
* GAM works better for one of the variable. We may predict that variable with GAM, and then run MICE for the 4 remaining variables. Or code our own MICE algorithm, integrating GAM for that variable.
* Think of a way to present the prediction results on the poster (beautiful plot).

Causal interpretation 
* Literature review for methods to build a DAG
- [Joint estimation of causal effects from observational and intervention gene expression data] (ftp://ftp.sam.math.ethz.ch/sfs/pub/Manuscripts/buhlmann/pcalg-jss.pdf): really not sure if that would work for us. Top of page 2: this is the method applied in the `load_mous_data.R` example R code provided with the challenge statement.
- This is the [paper] (ftp://ftp.sam.math.ethz.ch/sfs/pub/Manuscripts/buhlmann/pcalg-jss.pdf) introducing the `pcalg` package used in the example R code.
- It seems like [this](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12071/epdf) is what we are trying to do.









