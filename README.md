The current folder includes R code for reproducing all of the tables and figures in the article "Sample size and power calculation for testing treatment effect heterogeneity in cluster randomized crossover designs" by Wang et al. 

For questions or comments about the code please contact Xueqi Wang at xueqi.wang@yale.edu. 

I. List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript and supporting information.
1. functions_calc_ss.R = function to calculate the predicted number of clusters to get at least pre-specified power under a cluster randomized crossover (CRXO) design;
2. functions_gendata.R = function to generate the data;
3. functions_empirical.R = function to compute the empirical power or empirical type I error;
4. functions_calc_power.R = function to calculate the predicted power given number of clusters and mean/CV of cluster-period size;

II. List of Main Files: These main files are used to reproduce the results in the submitted manuscript and supporting information.
5. numStudy.R = reproduce results for the numerical study – roles of intracluster correlation coefficient (ICC) parameters;
6. VIFplot_CS.R = reproduce figures of variance inflation factor (VIF) for unequal to equal cluster-period sizes under the cross-sectional design;
7. VIFplot_CC.R = reproduce figures of VIF for unequal to equal cluster-period sizes under the closed-cohort design;
8. VIFplot_PCS.R = reproduce figures of relative efficiency (RE) for parallel cluster randomized trial (CRT) to cross-sectional CRXO;
9. VIFplot_PCC.R = reproduce figures of RE for parallel CRT to closed-cohort CRXO;
10. sim_CS_HTE_1.R = reproduce results for Simulations under the cross-sectional design;
11. sim_CC_HTE_1.R = reproduce results for Simulations under the closed-cohort design;
12. applications.R = reproduce results for Applications.
