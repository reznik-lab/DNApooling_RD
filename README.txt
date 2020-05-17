
#############################################################################
#############################################################################
##                                                                         ##
##       MULTIREGIONAL DNA POOLING OVERCOMES INTRATUMORAL HETEROGENEITY    ##
##                   IN CLINICAL SEQUENCING OF TUMORS                      ##
##                                                                         ##
##                            DiNatale RG et al.                           ##
##                                                                         ##
#############################################################################
#############################################################################


To run the analyses performed in this study you must first download the original directory which includes all the mutation-annotation files (MAFs) needed as well as two R scripts with the necessary code. You must then source the R scripts from the original directory to run the analyses.

Each script will create a set of subdirectories where the results will be placed. Results include graphs (in PDF formats) as well as tables and Rdata objects (SM3 only).  The code in SM3 can be further customized to use different simulation parameters. Additional details and usage below.

Contact: reznike@mskcc.org (for any questions related to this code)


###############################################################
#							      #
# 	SUPPLEMENTARY MATERIAL 2 - POOLING EXPERIMENT         #
#							      #
###############################################################


This R script uses multi-regional tumor mutation data (in MAF format) from the single-region and tumor DNA pools sequenced as part of the study and compares the outcomes of these approaches against a conventional 'multi-regional assessment'

The following outputs will be produced (panels of supplementary figure 2):

	1- An oncoprint showing the concordance in clonality calls for each sample (correct, incorrect, missed) (SF2A)
	2- The comparison of observed versus expected CCF and purity values in the DNA pools (SF2B)
	3- A summary of the outcomes for each type of sequencing approach (SF2C)
	4- The PR and ROC curves of the CCF in the pool to classify clonality (SF2D)
	5- The performance of each CCF cutoff in these pools (containing 6 regions) to classify clonality (SF2E)





#######################################################################################
#							          		      #
#       SUPPLEMENTARY MATERIAL 3 - BOOTSTRAPPING PROCEDURE ON TRACERx COHORTS         #
#									              #
#######################################################################################


This R script uses multi-regional tumor mutation data (in MAF format) from the TRACERx cohorts and simulates the confection of a tumor DNA pool using the formulas described in the manuscript (SM1) this is done by iteratively obtaining a subset of tumors from the cohort and sampling random regions from each one.

A process of bootstrapping with replacement is performed using increasing number of 
tumor regions to simulate these pools. The code then performs a comparison of several outcome measures between conventional multiregional and 'pseudo-bulk' pooled sequencing. 
The accuracy of the CCF to classify variants in clonal or subclonal is then assessed using several classifier accuracy measures such as the F1 score and the Matthew's Correlation Coefficient.
Finally, the code allows the user to specify a number of regions of interest and returns detailed results. 

Results are provided in the form of graphs (PDF format) as well as a table (CSV format) and an Rdata object that contains all the data from the bootstrapping iterations.





### USAGE ###


# INPUTS #

1) The mutation-annotation files (MAF) of the TRACERx RCC and NSCLC cohorts. 
    This data consists of multiregional somatic variant data from next-generation sequencing assays
     - RCC cohort - consists of targeted panel sequencing, 
     - NSCLC cohort - whole-exome sequencing data

     Mandatory columns: Hugo_Symbol, Tumor_Sample_Barcode, CCF (cancer-cell fraction), Purity, patient_ID, cohort, OncoKB_evidencelvl, mutation_ID (in the format GENE_Ch:StartPosition:Reference_Allele;PATIENTID), mutation_ID_region (in the format GENE_Ch:StartPosition:Reference_Allele;SAMPLEID)


2) A set of parameters to specify to performed the desired simulations
    - cohort = the cohort to use, choose from 'renal' and 'lung' (default:renal)

    - iter_num = number of samples (iterations) to use (default: 2)

    - max_region_num = maximum number of regions to include in the hypothetical pools (if greater than those available defaults to 5)

    - spec_region_num = specific number of regions in the pools, to be displayed more detailed results (if greater than those available defaults to 4)

    - sample_size = sample size of each of the iterations of the simulation, as a proportion of the total (default: 0.7)

    - Clonality definitions: a set of thresholds in the cancer-cell fraction (CCF) estimated in the pool (POOL_threshold, default: 0.75) and in the multiregional assessment (MR_threshold, default: 0.5) that will be used to classify mutations as clonal (>= threshold) or subclonal (< threshold).

    - detection_limit = specifies the CCF below which a mutation is likely to be missed in the pool (default: 0.02)

   - evidence_levels = Levels of evidence from the OncoKB platform that will be used to select variant for the targetable mutation analysis

   - costs = a list with the costs of sample profiling using NGS platforms
              (default: 'cost_seq' = 500 [cost of sample sequencing], 
                        'cost_persample_prep' = 50, [cost of DNA sample preparation],
                        'cost_overhead' = 1000, [overhead costs from processing a DNA sample],
                        'cost_mix' = 50 [hypothetical cost of mixing several DNA samples])


# OUTPUTS #

After editing the script and running the function inside the same repository as the MAF files
the script will produce the following outcomes:
  Creates a subdirectory in identified with a timestamp and the specified simulation parameters
  Contents:
      simresults_outcomes.csv - a table with the numerical outcomes and other details of each sample
      simresults_detailed.Rdata - a list that contains the merged detailed results from each simulation
      plots_outcomes - a directory with 14 plots showing the results of the analysis
      plots_classifiers - a directory with 20 plots showing the accuracy estimates of the pools in each simulation
      plots_detailed_poolXregions - a directory with the 4 plots showing the results of all the simulated pools with a number of regions equal to spec_region_num (4 plots)


