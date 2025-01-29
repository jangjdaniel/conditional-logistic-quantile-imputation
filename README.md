# truncated_missing_data
Spring 2025 (and potentially beyond) research on handling truncated missing data

Motivation:

The quality of instruments to measure biomarkers between different research facilities is variable due to various factors. For example, let’s say facility A has a lot of funding, and therefore has a machine that can fully detect a biomarker (called CB: a continuous variable with a limited range; in our case, it cannot be negative), while facility B gets less funding, so they have a machine that can only detect a concentration starting at 1 mg/L. Facility C on the other hand does not have a machine to obtain this biomarker value. 

Say each of these facilities performs a very similar analysis looking at some outcome Y with a treatment A, controlling for a set of confounders C, which contains CB. However, recall the different situations we had for each facility with CB. We would say that CB is Sporadically Missing for Facility B, meaning that a decent chunk of data is missing. Also note that this is data Missing Not at Random (MNAR), since we know the reason why we have missing data is not random. We would say that CB is Systematically Missing for Facility C, meaning the variable is fully missing (100% missing) from Facility C’s analysis. Let’s also say that we cannot pool the individual level data together from these three studies, potentially due to privacy issues or laws that prevent sharing this type of sensitive data. 

What is the proper course of action here if we want to perform a meta analysis? Do we fully remove CB from our analysis? Do we try to find similar studies to facility A so we can have access to this variable in its complete form? Either of these options is unrealistic. Instead, we propose an extension to (insert paper here) that uses existing methods for imputing missing truncated data using facility A’s analysis as a basis for imputation so we can keep facility B and C’s analysis, saving money and resources while increasing statistical power. 

The general process would be facility A sends over aggregated statistics to facility B and C for them to impute CB. Details for this process, as well a simulation study motivated by real-world scenarios will be the bulk of my work for this project.

I am not expecting to fully finish this work this semester, since it is quite dense and there is a slight learning curve with quantile regression, missing data methods, and meta analysis (the multivariate normal distribution hurts my brain). My goal is to make a significant dent in this question. My overall goal, however, is a publication (dream big). 
