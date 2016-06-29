# CausalChallenge2016
Code and notes for the "Causal Challenge" - July 2016

## Meeting June 22, 2016
- Gabrielle: Try pooling predicted values from individual MICE-imputed datasets to get a final set of predicted values.  Test this by deleting observations from another variable, imputing, and comparing predicted values to true values.  For now, do without litter variable, but will introduce later when we figure out how to handle it.
- Kevin:  Look into machine learning methods for clustering data.  The goal is to reduce the 190 litters into a variable with a much smaller number of categories by clustering similar litters together based on phenotype.
- Will hold off on causal part for now, since we need to get good predicted values before considering the problems in causal inference.
