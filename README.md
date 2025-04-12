# Conditional Logistic Quantile Regression

Motivation:

Biomarkers are often continuous and bounded variables that are expensive to measure. This can be especially problematic in clinical trials where financial restraints could result in not enough data being collected, or having to skip some patients's data collection. Instead of doing a complete case analysis or scrapping the experiment altogether, Conditional Logistic Quantile Imputation (CLQI) can be used so that we are able to perform the desired analysis without the loss of data. 

After introducing CLQI in this simple case, we will show how CLQI can handle limit-of-detection / left-censored / limited-range data with a few modifications. This proves the method to be flexible in imputing different types of missingness, but also meets the demands of researchers. 
