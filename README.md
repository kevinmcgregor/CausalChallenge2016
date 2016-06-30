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
