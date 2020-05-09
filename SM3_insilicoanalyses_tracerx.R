rm(list=ls())
setwd("/Users/dinatalr/Documents/MSKCC/rcc_translational/rcc_ith/analyses/scripts/final_code/")


#############################################################################
#       MULTIREGIONAL DNA POOLING OVERCOMES INTRATUMORAL HETEROGENEITY      #
#                   IN CLINICAL SEQUENCING OF TUMORS                        #
#                                                                           #
#                            DiNatale RG et al.                             #
#############################################################################

# The following code in R uses multiregional tumor mutation data and simulates
# the confection of a tumor DNA pool the formulas described in the manuscript (ST1)
# this is done by obtaining a series of random samples, each one smaller than the total cohort.

# A bootstrapping with replacement process is then performed using increasing number of 
# tumor regions to simulate these pools. The code then performs a comparison of several outcome measures
# between multiregional assessment and tumor pool assessment, and calculates several accuracy measures
# for the CCF in the pooled sample to classify mutations into clonal or subclonal.

# Finally, the code allows the user to specify a number of regions of interest and returns detailed results.
# Results are provided visually in the form of graphs (showing the different simulations)
# as well as a table and an Rdata object that contains all the detailed results


# USAGE

### INPUTS ###

# 1) The mutation-annotation files (MAF) of the TRACERx RCC and NSCLC cohorts. 
    # This data consists of multiregional somatic variant data from next-generation sequencing assays
        # RCC cohort - consists of targeted panel sequencing (XXX genes), 
        # NSCLC cohort - whole-exome sequencing data (XXXXX genes)

    # Mandatory columns: Hugo_Symbol, Tumor_Sample_Barcode, CCF (cancer-cell fraction), Purity, patient_ID
    #                    cohort, OncoKB_evidencelvl
    #                    mutation_ID (in the format GENE_Ch:StartPosition:Reference_Allele;PATIENTID), 
    #                    mutation_ID_region (in the format GENE_Ch:StartPosition:Reference_Allele;SAMPLEID)


# 2) A set of parameters to specify to performed the desired simulations
    # cohort = the cohort to use, choose from 'renal' and 'lung' (default:renal)

    # iter_num = number of samples (iterations) to use
                 # (default: 2)

    # max_region_num = maximum number of regions to include in the hypothetical pools 
                      # (if greater than those available defaults to 5)

    # spec_region_num = specific number of regions in the pools, to be displayed more detailed results
                      # (if greater than those available defaults to 4)

    # sample_size = sample size of each of the iterations of the simulation, as a proportion of the total 
                    # (default: 0.7)

    # CCF definitions: a set of thresholds in the cancer-cell fraction (CCF) estimated in the 
                          # pool (POOL_threshold, default: 0.75) and in the multiregional
                          # assessment (MR_threshold, default: 0.5) that will be used to classify mutations
                          # as clonal (>= threshold) or subclonal (< threshold).
                          # detection_limit = specifies the CCF below which a mutation is likely to be 
                                            # missed in the pool (default: 0.02)
    # evidence_levels = Levels of evidence from the OncoKB platform that will be used to select variant for the 
                        # targetable mutation analysis
    # costs = a list with the costs of sample profiling using NGS platforms
              # (default: 'cost_seq' = 500 [cost of sample sequencing], 
              #           'cost_persample_prep' = 50, [cost of DNA sample preparation],
              #           'cost_overhead' = 1000, [overhead costs from processing a DNA sample],
              #           'cost_mix' = 50 [hypothetical cost of mixing several DNA samples])



### OUTPUTS ###

# After editing the script and running the function inside the same repository as the MAF files
# the script will produce the following outcomes:
  # Creates a subdirectory in identified with a timestamp and the specified simulation parameters
  # Contents:
      # simresults_outcomes.csv - a table with the numerical outcomes and other details of each sample
      # simresults_detailed.Rdata - a list that contains the merged detailed results from each simulation
      # plots_outcomes - a directory with 14 plots showing the results of the analysis  (in PDF format)
      # plots_classifiers - a directory with 20 plots showing the accuracy estimates of the pools in each simulation
      # plots_detailed_poolXregions - a directory with the 4 plots showing the results of all the simulated pools 
                                    # with a number of regions equal to spec_region_num (4 plots)







#Install/load packages and basic functions
pkg_touse = c("readxl","dplyr","reshape2","ggplot2","cowplot",
              "PRROC","ggforce","grDevices","RColorBrewer") ; {
    pkg_toinstall = pkg_touse[!pkg_touse %in% rownames(installed.packages())]
    if(length(pkg_toinstall)>0){install.packages(pkg_toinstall)}
    lapply(pkg_touse, require, character.only = TRUE)}



#################################################################
###### SET UP DEFINITIONS AND PARAMETERS FOR BOOTSTRAPPING ######
#################################################################

#Cohort selection: 'renal' for RCC and 'lung' for NSCLC TRACERx cohorts
cohort_select = "lung" ; {
    if(tolower(cohort_select) == "renal"){maf = read.csv("ST4_MAF_TRACERx_RCC.csv",stringsAsFactors = F) ; cohort = "RCC"}
    if(tolower(cohort_select) == "lung"){maf = read.csv("ST5_MAF_TRACERx_NSCLC.csv",stringsAsFactors = F) ; cohort = "NSCLC"}

    #Split MAF into 'nested' list of tumors with multiple regions     
    list.pats = split(maf, maf$patient_ID) # Tumor-level list
    
    # Region-level: each element is a list of data frames with region-specific variants
    list.regs = sapply(list.pats, function(x) split(x, x$Tumor_Sample_Barcode))} 


#Detection limit in a tumor DNA pool CCF threshold below which a variant might be missed in a tumor DNA pool
detection_limit = 0.02

#Clonality definitions: CCF to define a variant as 'clonal' (>= threshold, which has to be >0)
MR_threshold = 0.5 #Multiregional threshold (must meet condition in all regions)
POOL_threshold = 0.75 # Tumor DNA pool threshold



# Evidence levels to consider when evaluating targetable alterations
evidence_levels = 1:3 
    
#Estimated costs
costs = list('cost_seq' = 500, # cost of sample sequencing
             'cost_prep_persample' = 50, # cost of DNA sample preparation
             'cost_overhead' = 1000, # overhead costs from processing a DNA sample (include bioinformatics analysis)
             'cost_mix' = 50) # added cost of tumor DNA pool confection (i.e. DNA aliquot mixing)   

#Relative outcomes: outcomes to be displayed relative to single-region profiling
relative_outcomes = c("cost","cea_detection","cea_clonality","targets")  

#Bootstrapping parameters
#Sample size: proportion of tumors to be selected (out of the total with that amount of regions available)
sample_size = 0.7   
iter_num = 100 #Number of iterations
max_region_num = 5 # Number of regions to include in the analysis (all possibilities from 1:max_region_num will be included)
spec_region_num = 4 # Specifies the specific simulation where more detailed results will be shown (must be less than max_region_num)

    # If the number of regions selected is not available, use 5 regions as maximum, and 4 for detailed results
    if(max_region_num > max(sapply(list.regs, length))){max_region_num = 5} 
    if(spec_region_num > max(sapply(list.regs, length))){spec_region_num = 4} 


    #Create output directory to store simulation results
    time_of_simulation = gsub("[[:punct:]]","",Sys.time())
        time_of_simulation = gsub(" ","_",time_of_simulation) #Timestamp in the format: YYYYMMDD_hhmmss
        output_path = paste0(time_of_simulation,"_",cohort,"_ccf.MR",MR_threshold*100,".POOL",POOL_threshold*100,"_i",iter_num,"r",max_region_num,"/")
        system(paste("mkdir",output_path))
        
        
        #Basic functions used in the body of the simulation and true results (i.e. risk, subtypes): 
        #callClonal_MR,callClonal_POOL,genomicRiskRCC, subtypes_txRCC_MR, subtypes_txRCC_POOL, compute_mean95ci
        {
          
          # Clonality functions - to use in a 'nested' list of tumor regions.
          # Given a specific variant, the functions take a numeric vector of CCF values in each region (0 or NA for variants not detected)
          # The POOL function takes an additional argument 'pur', a vector of purity values in each region  
          # as well as 'detection limit', the CCF value below which a variant would likely not be detected
          
          callClonal_MR = function(ccfs, cutpoint = MR_threshold){
            if(all(is.na(ccfs))){result = "miss"}else{
              if(any(is.na(ccfs) | ccfs < cutpoint)){result = "subclonal"
              }else{result = "clonal"}}
            return(result)}
          
          callClonal_POOL = function(ccfs, pur, cutpoint = POOL_threshold, detect_lim = detection_limit){
            if(length(ccfs) != length(pur)){break("Both vectors must be of the same length")}
            if(all(is.na(ccfs))){result = 0 ; names(result) = "miss"}else{
              if(any(is.na(ccfs))){ccfs[is.na(ccfs)] = 0}
              
              ccfpool = sum(ccfs*pur)/sum(pur) # Formula of ccf in the tumor DNA pool (sum(CCF*p)/sum(p))
              if(ccfpool < detect_lim){result = 0 ; names(result) = "miss"}else{
                result = ccfpool
                names(result) = ifelse(ccfpool >= cutpoint,"clonal","subclonal")
              }
            }
            return(result)
          }
          
          
          #Functions for translational outcomes, RCC specific 
          #RCC risk stratification (Voss MH, et al. Lancet Oncol 2018)  
          genomicRiskRCC = function(muts){ 
            #takes a character vector with gene names 
            criterion_a = ifelse(any(muts %in% c("BAP1","TP53")),1,0)
            criterion_b = ifelse(all(muts %in% c("BAP1","TP53","PBRM1")),1,0)
            criterion_c = ifelse(!any(muts %in% "PBRM1"),1,0)
            points_bc = ifelse(criterion_b == 1 | criterion_c == 1,1,0)
            
            points = criterion_a + points_bc
            atrisk = ifelse(points>=1,1,0)
            res = c(criterion_a, criterion_b,criterion_c,
                    points, atrisk)
            names(res) = c("bap1.tp53","multidriver","pbrm1wt",
                           "risk_points","atrisk")
            return(res)
          }
              true_risk = sapply(split(maf, maf$patient_ID), function(x) genomicRiskRCC(unique(x$Hugo_Symbol)) )
              true_risk = as.data.frame(t(true_risk),stringsAsFactors = F)
          
          
          #RCC subtype assignment (Turajlic S, et al. Cell 2018)
          driver_events_rcc = c("VHL", "PBRM1", "SETD2", "PIK3CA", "MTOR", "PTEN", "KDM5C", "CSMD3", "BAP1", "TP53", "TSC1", "TSC2", "ARID1A", "TCEB1")
          
          subtypes_txRCC_MR = function(multiregion_list, clonal_threshold = MR_threshold, 
                                       gene.col = "Hugo_Symbol", pur.col = "Purity", ccf.col = "CCF"){
            
            #Make sure data is in the right format
            if(!is.list(multiregion_list)){stop("Please provide a list with the mutations per region")}
            if(!is.numeric(clonal_threshold) | !(clonal_threshold>0 & clonal_threshold<1)){stop("please provide a clonal")}
            
            #Genes
            driver_events = c("VHL", "PBRM1", "SETD2", "PIK3CA", "MTOR", "PTEN", "KDM5C", "CSMD3", "BAP1", "TP53", "TSC1", "TSC2", "ARID1A", "TCEB1")
            pi3k_pathway = c("PIK3CA","MTOR","PTEN","TSC1","TSC2")
            core_drivers = c("VHL", "PBRM1", "SETD2","BAP1","PTEN") #-------------------------- change 1, core events includes PTEN as the only representative of PI3K pathway
            #core_drivers = c("VHL", "PBRM1", "SETD2","BAP1",pi3k_pathway) 
            
            #Vector of drivers (clonal/subclonal/miss) - if tandem muts, selects the one with highest CCF in that region
            
            #matrix of ccf per region for each driver event, 
            #if the gene is not present in a region, return NA, 
            #if there is >1 mutation  in the gene, pick the one with the highest CCF (in that region)
            ccfs = sapply(driver_events, function(g){
              value_pergene =  sapply(multiregion_list, function(x){
                ifelse(g %in% x[,gene.col], 
                       max(x[,ccf.col][x[,gene.col] == g]), NA)})
              if(length(multiregion_list) == 1){names(value_pergene) = NULL}
              return(value_pergene)
            })
            ccfs = matrix(ccfs, nrow = length(multiregion_list), ncol = length(driver_events), 
                          dimnames = list(names(multiregion_list),driver_events))
            
            
            #obtain purity for each region
            pur = sapply(multiregion_list, function(x) unique(x[,pur.col]))
            
            #Make an assertion about the clonality of each variant
            # clonality_calls = a named character vector containing a 'clonal status' for each of the genes (clonal/subclonal/miss) 
            clonality_calls = apply(ccfs,2,callClonal_MR, cutpoint = clonal_threshold)
            
            #events observed
            observed = clonality_calls[clonality_calls != "miss"]
            clonal_events = names(observed[observed == "clonal"])
            subclonal_events = names(observed[observed == "subclonal"])
            
            
            #Rules in hierarchical order -> Swanton's definitions
            subtype = NA
            
            #1 - Multiple clonal driver (>= clonal muts in BAP1, PBRM1, SETD2 or PTEN)
            num_clonal_drivers = sum(c("BAP1","PBRM1","SETD2","PTEN") %in% clonal_events)
            if(num_clonal_drivers >=2){subtype = "Multiple_Clonal_Drivers"}
            
            #2 - BAP1-driven tumors
            bap1_driven = ("BAP1" %in% names(observed) & is.na(subtype))
            if(bap1_driven){
              #core_events_novhl.bap1 = core_drivers[which(core_drivers %in% names(observed) & !core_drivers %in% c("VHL"))]
              core_events_novhl.bap1 = core_drivers[which(core_drivers %in% names(observed) & !core_drivers %in% c("VHL","PBRM1"))]#------------ change 2, BAP1 driven if BAP1 sole driver other than VHL and PBRM1 (change refers to addition of the latter)
              
              
              #If BAP1 is single clonal driver other than VHL, subtype is BAP1-driven
              if(length(core_events_novhl.bap1) == 1){
                subtype = "BAP1_driven"
              }else{
                subclone_bap1status = apply(ccfs,1,function(y){
                  
                  y = y[core_events_novhl.bap1]
                  y = y[!is.na(y)]
                  bap1 = y["BAP1"] 
                  
                  
                  if(is.na(bap1) | is.null(bap1) | length(bap1) == 0){
                    res="other"
                  }else{
                    rest = y[!names(y) %in% "BAP1"]
                    if(bap1 >= clonal_threshold & all(rest < clonal_threshold)){res = "bap1_driven_clone"}else{
                      res = "other"}
                  }
                  return(res)})  
                subclone_bap1status = subclone_bap1status[!is.na(subclone_bap1status)]
                subtype = ifelse(any(subclone_bap1status == "bap1_driven_clone"), "BAP1_driven",NA)
              }
            }
            
            #3-5 PBRM1-driven tumors
            #Includes groups 3(PBRM1>SETD2), 4(PBRM1>PIK3CA), 5(PBRM1>SCNA, indeterminate here),
            pbrm1_driven = ("PBRM1" %in% names(observed) & is.na(subtype))
            if(pbrm1_driven){
              
              subclone_pbrm1status = apply(ccfs,1,function(y){
                relevant_genes = c("PBRM1","SETD2",pi3k_pathway)
                #core_events = core_events
                
                y = y[relevant_genes]
                y = y[!is.na(y)]
                pbrm1 = ifelse(all(names(y) == "PBRM1"),y,y["PBRM1"]) 
                
                res = NA
                
                #No PRM1 mutation in that subclone makes the result = NA
                if(is.na(pbrm1) | is.null(pbrm1) | length(pbrm1) == 0 | length(y) == 0){
                  res=NA
                }else{
                  
                  # 3 - PBRM1 --> SETD2
                  pbrm1_setd2 = ('SETD2' %in% names(y) & is.na(res))
                  if(pbrm1_setd2){
                    setd2 = y['SETD2']
                    if(pbrm1 > setd2){res = "PBRM1-->SETD2"}
                    if(pbrm1 > clonal_threshold & setd2 <= clonal_threshold){res = "PBRM1-->SETD2"}
                    if(pbrm1 > clonal_threshold & setd2 > clonal_threshold){res = "PBRM1-->SETD2"}
                    if("SETD2" %in% clonal_events){res = NA}
                    
                  }
                  
                  #4 - PBRM1 --> PI3K
                  pbrm1_pi3k = (any(pi3k_pathway %in% names(y)) & is.na(res))
                  if(pbrm1_pi3k){
                    pi3k = y[names(y) %in% pi3k_pathway]
                    if(pbrm1 > clonal_threshold & max(pi3k) < clonal_threshold){res = "PBRM1-->PI3K"}
                    if(pbrm1 < clonal_threshold & max(pi3k) < clonal_threshold & pbrm1 > max(pi3k)){res = "PBRM1-->PI3K"}
                    
                    if(pbrm1 > clonal_threshold & max(pi3k) > clonal_threshold){res = "PBRM1-->PI3K"}
                    if(any(names(pi3k) %in% clonal_events)){res = NA}
                    
                  }
                  
                  
                  #5 - PBRM1 --> SCNA
                  if(is.na(res)){
                    res = "PBRM1-->sCNA"
                  }
                }
                return(res)
              })
              subclone_pbrm1status = unique(subclone_pbrm1status[which(!is.na(subclone_pbrm1status) & subclone_pbrm1status != "other")])
              subclone_pbrm1status = factor(subclone_pbrm1status, levels = c("PBRM1-->SETD2","PBRM1-->PI3K","PBRM1-->sCNA"))
              
              if(length(subclone_pbrm1status) <= 1){
                subtype = subclone_pbrm1status
                subtype = ifelse(length(subclone_pbrm1status) == 0, NA,as.character(subclone_pbrm1status))
              }else{
                subtype = as.character(subclone_pbrm1status[which.min(subclone_pbrm1status)])
              }
            }
            
            #6 - VHL wild-type
            vhlwt = (!"VHL" %in% names(observed) & (is.na(subtype) | subtype == "PBRM1-->sCNA")) #------------ change 3, if PBRM1-CN/indeterminate, still evaluate for VHL-wt group and prioritize VHL-wt over PBRM1-->CN
            if(vhlwt){subtype = "VHL_wt"}
            
            #7 - VHL monodriver
            vhlmono = ("VHL"%in% names(observed) & length(observed) == 1 & is.na(subtype))
            if(vhlmono){subtype = "Monodriver_VHL"}
            
            #8 - Unknown
            if(is.na(subtype)){subtype = "non_driver_subtype"}
            
            subtype = factor(subtype, levels = c("Multiple_Clonal_Drivers","BAP1_driven","PBRM1-->SETD2","PBRM1-->PI3K","PBRM1-->sCNA","VHL_wt","Monodriver_VHL","non_driver_subtype"))
            return(subtype)
          }
          subtypes_txRCC_POOL = function(multiregion_list, clonal_threshold = POOL_threshold, ccf_clonality_margin = 0.02, 
                                         gene.col = "Hugo_Symbol", pur.col = "Purity", ccf.col = "CCF"){
            
            #Make sure data is in the right format
            if(!is.list(multiregion_list)){stop("Please provide a list with the mutations per region")}
            if(!is.numeric(clonal_threshold) | !(clonal_threshold>0 & clonal_threshold<1)){stop("please provide a clonal")}
            
            
            #Genes
            driver_events = c("VHL", "PBRM1", "SETD2", "PIK3CA", "MTOR", "PTEN", "KDM5C", "CSMD3", "BAP1", "TP53", "TSC1", "TSC2", "ARID1A", "TCEB1")
            pi3k_pathway = c("PIK3CA","MTOR","PTEN","TSC1","TSC2")
            core_drivers = c("VHL", "PBRM1", "SETD2","BAP1","PTEN") #-------------------------- change 1, core events includes PTEN as the only representative of PI3K pathway
            #core_drivers = c("VHL", "PBRM1", "SETD2","BAP1",pi3k_pathway) 
            
            #Vector of drivers (clonal/subclonal/miss) - if tandem muts, selects the one with highest CCF in that region
            
            #matrix of ccf per region for each driver event, 
            #if the gene is not present in a region, return NA, 
            #if there is >1 mutation  in the gene, pick the one with the highest CCF (in that region)
            ccfs = sapply(driver_events, function(g){
              value_pergene =  sapply(multiregion_list, function(x){
                ifelse(g %in% x[,gene.col], 
                       max(x[,ccf.col][x[,gene.col] == g]), NA)})
              if(length(multiregion_list) == 1){names(value_pergene) = NULL}
              return(value_pergene)
            })
            ccfs = matrix(ccfs, nrow = length(multiregion_list), ncol = length(driver_events), 
                          dimnames = list(names(multiregion_list),driver_events))
            
            #obtain purity for each region
            pur = sapply(multiregion_list, function(x) unique(x[,pur.col]))
            
            #Expected CCF in the pooled sample
            ccf_pool = apply(ccfs,2, callClonal_POOL,pur = pur, cutpoint = clonal_threshold)
            ccf_pool = ccf_pool[ccf_pool != 0]
            
            
            #Make an assertion about the clonality of each variant
            # clonality_calls = a named character vector containing a 'clonal status' for each of the genes (clonal/subclonal/miss) 
            clonality_calls = apply(ccfs,2, function(x)
              names(callClonal_POOL(x,pur = pur, cutpoint = clonal_threshold)))
            
            #events observed
            observed = clonality_calls[clonality_calls != "miss"]
            clonal_events = names(observed[observed == "clonal"])
            subclonal_events = names(observed[observed == "subclonal"])
            
            #Rules in hierarchical order -> Swanton's definitions
            subtype = NA
            
            #1 - Multiple clonal driver (>= clonal muts in BAP1, PBRM1, SETD2 or PTEN)
            num_clonal_drivers = sum(c("BAP1","PBRM1","SETD2","PTEN") %in% clonal_events)
            if(num_clonal_drivers >=2){subtype = "Multiple_Clonal_Drivers"}
            
            
            #2 - BAP1-driven tumors
            bap1_driven = ("BAP1" %in% names(observed) & is.na(subtype))
            if(bap1_driven){
              #core_events_novhl.bap1 = core_drivers[which(core_drivers %in% names(observed) & !core_drivers %in% c("BAP1","VHL"))]
              core_events_novhl.bap1 = core_drivers[which(core_drivers %in% names(observed) & !core_drivers %in% c("BAP1","VHL","PBRM1"))]#------------ change 2, BAP1 driven if BAP1 sole driver other than VHL and PBRM1 (change refers to addition of the latter)
              core_ccf = ccf_pool[core_events_novhl.bap1]
              bap1 = ccf_pool["BAP1"] 
              
              #If BAP1 is single clonal driver other than VHL, subtype is BAP1-driven
              if(length(core_ccf) == 0 | "BAP1" %in% clonal_events){
                subtype = "BAP1_driven"
              }else{
                if(bap1 > max(core_ccf)+ccf_clonality_margin){subtype = "BAP1_driven"}    
              }
            }
            
            
            # PBRM1-driven tumors
            #Includes groups 3(PBRM1>SETD2), 4(PBRM1>PIK3CA), 5(PBRM1>SCNA, indeterminate here),
            pbrm1_driven = ("PBRM1" %in% names(observed) & is.na(subtype))
            if(pbrm1_driven){
              
              relevant_genes = c("PBRM1","SETD2",pi3k_pathway)
              
              relevant = ccf_pool[relevant_genes]
              pbrm1 = ccf_pool['PBRM1']
              setd2 = ccf_pool['SETD2']
              pi3k = ccf_pool[names(ccf_pool) %in% pi3k_pathway]
              
              # 3 - PBRM1 --> SETD2
              pbrm1_setd2 = ('SETD2' %in% names(ccf_pool) & is.na(subtype))
              if(pbrm1_setd2){
                if("PBRM1" %in% clonal_events & !'SETD2' %in% clonal_events){subtype = "PBRM1-->SETD2"}
                if(!"PBRM1" %in% clonal_events & !"SETD2" %in% clonal_events & 
                   pbrm1 > setd2+ccf_clonality_margin){subtype = "PBRM1-->SETD2"}
              }
              
              #4 - PBRM1 --> PI3K
              pbrm1_pi3k = (any(pi3k_pathway %in% names(ccf_pool)) & is.na(subtype))
              if(pbrm1_pi3k){
                if("PBRM1" %in% clonal_events & !any(pi3k_pathway %in% clonal_events)){subtype = "PBRM1-->PI3K"}
                if(!"PBRM1" %in% clonal_events & !any(pi3k_pathway %in% clonal_events) & 
                   pbrm1 > max(pi3k)+ccf_clonality_margin){subtype = "PBRM1-->PI3K"}
              }
              
              #5 - PBRM1 --> SCNA
              if(is.na(subtype) & 'PBRM1' %in% names(observed)){
                subtype = "PBRM1-->sCNA"
              }
            }
            
            
            #6 - VHL wild-type
            vhlwt = (!"VHL" %in% names(observed) & (is.na(subtype) | subtype == "PBRM1-->sCNA")) #------------ change 3, if PBRM1-CN/indeterminate, still evaluate for VHL-wt group and prioritize VHL-wt over PBRM1-->CN
            if(vhlwt){subtype = "VHL_wt"}
            
            
            #7 - VHL monodriver
            vhlmono = ("VHL"%in% names(observed) & length(observed) == 1 & is.na(subtype))
            if(vhlmono){subtype = "Monodriver_VHL"}
            
            
            #8 - Unknown
            if(is.na(subtype)){subtype = "non_driver_subtype"}
            
            
            subtype = factor(subtype, levels = c("Multiple_Clonal_Drivers","BAP1_driven","PBRM1-->SETD2","PBRM1-->PI3K","PBRM1-->sCNA","VHL_wt","Monodriver_VHL","non_driver_subtype"))
            return(subtype)
          }
              true_subtypes = sapply(list.regs, subtypes_txRCC_MR, clonal_threshold = MR_threshold)
          
          #Function to calculate classifiers - takes two vectors of the same length with the same labels
          calculateClassifiers =  function(test, true, res = "classifiers", lvls = c("clonal","subclonal")){
            
                test = factor(test, levels = lvls)
                true = factor(true, levels = lvls)
                
                test_labels = names(table(test))
                true_labels = names(table(true))
                
                conf_matrix = as.data.frame(table(interaction(test, true,drop = F)),stringsAsFactors = F)
                confusion = matrix(0, nrow=2, ncol = 2, 
                                   dimnames = list('test' = test_labels, 'truth' = true_labels))
                #Fill the empty matrix
                for(cat in 1:nrow(conf_matrix)){
                  val = conf_matrix$Freq[cat]
                  testdata = strsplit(conf_matrix$Var1[cat],split = "\\.")[[1]][1]
                  truth = strsplit(conf_matrix$Var1[cat],split = "\\.")[[1]][2]
                  confusion[testdata,truth] = val}
                
                #labels = "clonal"    "subclonal"
                tp = confusion[test_labels[1],true_labels[1]] 
                fp = confusion[test_labels[1],true_labels[2]]
                tn = confusion[test_labels[2],true_labels[2]] 
                fn = confusion[test_labels[2],true_labels[1]]
                
                #CLONALITY CASSIFIER MEASURES 
                #Sources: Fawcett (2006), Powers (2011), Ting (2011), and CAWCR[4] Chicco & Jurman (2020)
                classifiers = list()
                
                    #Curves estimates - PRROC package
                    classifiers$auc_roc = roc$auc
                    classifiers$auc_pr.integral = pr$auc.integral            
                    classifiers$auc_pr.dg = pr$auc.davis.goadrich
                    
                    # Pre-test measures
                    classifiers$sensitivity_tpr = tp/(tp+fn) # proportion of clonal events correctly classified
                    classifiers$specificity_tnr = tn/(tn+fp) # proportion of subclonal events correctly classified
                    classifiers$missrate_fnr = fn/(fn+tp) # (1-sens) proportion of clonal missclasified as subclonal 
                    classifiers$fallout_fpr = fp/(fp+tn) # (1-spec) proportion of false missclasified as clonal
                    
                    # Post-test measures
                    classifiers$ppv_precision = tp/(tp+fp) # proportion of real clonal events out of all clonal calls
                    classifiers$npv_prevention = tn/(tn+fn) # proportion of real subclonal events out of all subclonal calls
                    classifiers$fdr = fp/(fp+tp) # (1-PPV) proportion of incorrect out of all the clonal calls
                    classifiers$for_omission = fn/(fn+tn) # (1-NPV) proportion of incorrect out of all the subclonal calls
                    
                    #Other measures of classifier accuracy
                    classifiers$informedness_bm = classifiers$sensitivity + classifiers$specificity -1
                    classifiers$markedness = classifiers$ppv_precision + classifiers$npv - 1
                    classifiers$accuracy = (tp+tn)/sum(confusion)
                    classifiers$accuracy_bal = (classifiers$sensitivity + classifiers$specificity)/2
                    
                    classifiers$f1 = 2 * (classifiers$ppv_precision * classifiers$sensitivity)/
                                          (classifiers$ppv_precision + classifiers$sensitivity)
                    classifiers$clonalityillusion_fprate = fp/sum(confusion)
                    
                    #Cohen's Kappa coefficient ( Po - Pe / 1 - Pe)
                    prandom_clonal = ( (tp+fp)/sum(confusion) ) * ( (tp+fn)/sum(confusion) )
                    prandom_subclonal = ( (tn+fn)/sum(confusion) ) * ( (tn+fp)/sum(confusion) )
                    pexpect = prandom_clonal + prandom_subclonal
                    classifiers$kappa_cohen = (classifiers$accuracy - pexpect) / (1-pexpect)
                    
                    #Fowlkes–Mallows index
                    classifiers$fowlkes_mallows = sqrt(classifiers$ppv_precision * classifiers$npv_prevention)
                    
                    #Jaccard similarity coefficient and its distance measure
                    classifiers$jacsimilarity = tp/(tp+fp+fn) # Jaccard Index (J), Threat Score (TS), Critical Success Index (CSI) or Gilbert Index (G)
                    classifiers$jacdistance = (fp+fn)/(tp+fp+fn) # 1-J, a disagreement measure (higher = more disagreement)
                    
                    #Matthew's correlation coefficient (MCC) 
                    #PMID 22905111 - PLOS One A Comparison of MCC and CEN Error Measures in Multi-Class Prediction
                    classifiers$mcc =          ( (tp*tn) - (fp*fn) )/
                                      sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn))
                
                
                classifiers = unlist(classifiers)
                
                if(res == "classifiers"){results = classifiers}
                if(res == "confusion"){results = confusion}
                return(results)
                
              }
              
              
          #Function to calculate the 95%CI of a mean
          compute_mean95ci = function(values){
            if(any(is.na(values))){values = values[!is.na(values)]}
            
            n = length(values)
            mean = mean(values, na.rm = T)
            se = sd(values)/sqrt(n)
            error = qnorm(0.975)*se
            
            res = c(mean, mean-error, mean+error)
            names(res) = c("mean","low95","up95")
            return(res)
            
          }
          
          #Function to make piecharts, takes a vector of discrete values
          plot_pie = function(vec, plot_title = NULL, colpal = NULL, path = NULL, output = 'plot', plot_type = "pie"){
            
            if(!is.character(vec) & !is.factor(vec)){stop('A character or factor must be provided')}
            
            #Setting up the data
            dat = as.data.frame(table(vec),stringsAsFactors = F)
            if(is.factor(vec)){
              dat$vec = factor(dat$vec, levels = levels(vec))
            }else{
              dat$vec = as.factor(dat$vec)
            }
            dat = dat[order(dat$vec),]
            dat$prop = dat$Freq/sum(dat$Freq)
            dat$cumprop = cumsum(dat$prop)
            dat$label = paste0(dat$vec,"\n", round(dat$prop*100,1),"%")
            for(r in 1:nrow(dat)){
              basepoint = ifelse(r == 1, 0, dat$cumprop[r-1])
              #if(ind == 1){
              #  basepoint = 0}else{basepoint = dat[ind-1,'cumprop']}
              position = basepoint + as.numeric(dat$prop[r])/2
              position = 1 - position
              dat[r,'placeholder_label'] = position
            }
            
            #Aesthetics
            plot_title = ifelse(is.null(plot_title),"",plot_title)
            if(is.null(colpal)){
              colpal = colorRampPalette(brewer.pal(8,"Set3"))(length(levels(vec)))}
            
            #Plotting
            plot =  ggplot(dat) +
              geom_bar(aes(fill = vec, x = colnames(dat)[1], y = Freq),stat = "identity", position = "fill", color = "#000000") +
              geom_text(aes(x=colnames(dat)[1], y= placeholder_label, label = label), size = 6, color = "#000000") +
              scale_fill_manual("",values = colpal) +
              ggtitle(plot_title) +
              xlab("") +  ylab("") +
              guides(fill = guide_legend(reverse = TRUE)) +
              theme_classic() +
              theme(plot.title = element_text(size=18, color = "black", face="bold", hjust = 0.5), # controls the aesthetics of the plot
                    axis.text.x = element_blank(), 
                    axis.title.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.title.y = element_blank(),
                    legend.text = element_text(size=16, color = "black"),
                    legend.title = element_text(size = 16, face = "bold", color = "black"),
                    legend.position = "",
                    axis.ticks = element_blank(),
                    axis.line = element_blank())
            if(plot_type == "pie"){
              plot = plot + coord_polar(theta="y",start=0, direction=1)
            }
            
            if(!is.null(path)){
              pdf(path,7,7)
              print(plot)
              dev.off()
            }
            
            if(output == "plot"){return(plot)}
            if(output == "data"){return(dat)}
          }
          
          
          
        }        
        



#########################
###### SIMULATIONS ######
#########################
#Create dataframe to store bootstrapping results
bootstrap = data.frame()
details = list()

    # Loop through all iterations of the bootstrapping 
    # (each simulation is done with a specific number of regions)
    for(region_num in 1:max_region_num){ 
    for(i in 1:iter_num){
        
        #Simulation id to store results
        iteration = paste0("i",i)
        sim_id = paste0(iteration,"_r",region_num) 
        
        #Bootstrapping: obbtain a subset of tumors (with replacement) out of those with enough regions profiled   
        tum_enough_samps = which(sapply(list.regs, length) >= region_num)
            tum_enough_samps = names(list.regs)[tum_enough_samps]    
            set.seed(i)    
            n_tumors = round(sample_size*length(tum_enough_samps),0)
            tumors_selected = sample(tum_enough_samps, n_tumors, replace = T)
    
            
        #EVENT-LEVEL DATA
        #Store event-level results in an empty dataframe, loops through every tumor
        res_eventlvl = lapply(tumors_selected, function(tumor){
          
          #Sample a random subset of regions from each tumor
          t = which(tumors_selected == tumor) 
          regions_all = list.regs[[tumor]]
          set.seed(i)
          region_select = sample(names(regions_all))[1:region_num]
          regions = regions_all[region_select] # Nested list of tumor regions (region-specific mutations)
          regions_chosen = paste(region_select,collapse=";")
          
          #All somatic mutations present in the tumor
          all_muts = sapply(regions_all, function(x) x$mutation_ID)
          all_muts = unique(as.character(unlist(all_muts)))
          eventids = paste0(all_muts,"_t",t)
          tumorid = paste0(tumor,"_t",t)
          gene_events = sapply(strsplit(all_muts,split="_"), function(x)x[1])
          
          
          #Extracts maf ids and OncoKB info
          mafindex = which(maf$mutation_ID %in% all_muts)
          oncokb = maf[mafindex,c("mutation_ID","OncoKB_evidencelvl")]
          oncokb = oncokb[!duplicated(oncokb$mutation_ID),]
          evidence_targets = oncokb$OncoKB_evidencelvl ; 
          names(evidence_targets) = oncokb$mutation_ID
          
          
          ### MUTATION DETECTION ###
          detection_events_mr = sapply(all_muts,function(event){
            detect = sapply(regions, function(x){
              res = ifelse(event %in% x$mutation_ID,1,0)
              return(res)})
            detect = ifelse(any(detect == 1),1,0)
            return(detect)})
          
          detection_events_pool = sapply(all_muts, function(event){
            ccf = sapply(regions, function(x){
              if(event %in% x$mutation_ID){
                res = x$CCF[x$mutation_ID == event]
              }else{res = NA}
              return(res)})
            pur = sapply(regions, function(x){unique(x$Purity)})
            
            pool_status = callClonal_POOL(ccf, pur, detect_lim = detection_limit, cutpoint = POOL_threshold)
            if(region_num == 1){ # If single region, use the MR definitions instead
              pool_status = callClonal_POOL(ccf, pur, detect_lim = detection_limit, cutpoint = MR_threshold)}
            detect = ifelse(names(pool_status) != "miss",1,0) ; return(detect)})
          
          
          ### CORRECT CLONALITY CALLS ###
          #Determine the true clonality (with all regions available)
          clonality_truth = sapply(all_muts, function(event){
            ccf = sapply(regions_all, function(x){
              if(event %in% x$mutation_ID){
                res = x$CCF[x$mutation_ID == event]}else{res = NA}
              return(res)})
            clonality = callClonal_MR(ccf, cutpoint = MR_threshold)})    
          
          clonality_mr = sapply(all_muts, function(event){
            ccf = sapply(regions, function(x){
              if(event %in% x$mutation_ID){
                res = x$CCF[x$mutation_ID == event]}else{res = NA}
              return(res)})
            clonality = callClonal_MR(ccf, cutpoint = MR_threshold)})   
          
          
          clonality_pool = sapply(all_muts, function(event){
            ccf = sapply(regions, function(x){
              if(event %in% x$mutation_ID){
                res = x$CCF[x$mutation_ID == event]}else{res = NA}
              return(res)})
            pur = sapply(regions, function(x){unique(x$Purity)})
            pool_status = callClonal_POOL(ccf, pur, detect_lim = detection_limit, cutpoint = POOL_threshold)
            if(region_num == 1){pool_status = callClonal_POOL(ccf, pur, detect_lim = detection_limit, cutpoint = MR_threshold)}
            return(pool_status)})
          ccf_pool = clonality_pool
          pur_pool = sapply(all_muts, function(event){
            pur = sapply(regions, function(x){unique(x$Purity)})
            res_eventlvl = mean(pur, na.rm = T)
            return(res_eventlvl)})
          clonality_pool = sapply(strsplit(names(clonality_pool),split="\\.", fixed = F), function(x) x[2])
          
          
          
          resall = list('mutation_ID' = all_muts, 'tumor_ID'= rep(tumorid,length(all_muts)),
                        'tumor_ID_raw' = rep(tumor,length(all_muts)),
                       'Gene' = gene_events, 'regions' = rep(regions_chosen,length(all_muts)), 'OncoKB_evidencelvl' = evidence_targets[all_muts],
                       'mr_detection' = detection_events_mr, 'pool_detection' = detection_events_pool, 
                       'true_clonality' = clonality_truth, 'mr_clonality' = clonality_mr,'pool_clonality' = clonality_pool,
                       'pool_ccf' = ccf_pool, 'pool_pur' = pur_pool)
                          
          
          
          #RCC SUBTYPE ASSIGNMENT
          if(cohort == "RCC"){    
            pool_ccf_threshold = ifelse(region_num == 1, MR_threshold,POOL_threshold)
            subtype_mr = subtypes_txRCC_MR(multiregion_list = regions,clonal_threshold = MR_threshold)
            subtype_pool = subtypes_txRCC_POOL(multiregion_list = regions, clonal_threshold = pool_ccf_threshold)
            resall[['mr_subtype']] = rep(as.character(subtype_mr),length(all_muts))
            resall[['pool_subtype']] = rep(as.character(subtype_pool),length(all_muts))
            }
          
          names_resall = names(resall)
          resall = Reduce(cbind,resall)
          colnames(resall) = names_resall
          return(resall)
        })
            colnames_eventlvl = colnames(res_eventlvl[[1]])
            res_eventlvl = as.data.frame(Reduce(rbind,res_eventlvl),stringsAsFactors = F)
            colnames(res_eventlvl) = colnames_eventlvl
            rownames(res_eventlvl) = 1:nrow(res_eventlvl)
            numeric_vars_reseventlvl = grep('_detection|_ccf|_pur', colnames_eventlvl)
            res_eventlvl[,numeric_vars_reseventlvl] = sapply(res_eventlvl[,numeric_vars_reseventlvl], as.numeric)
        
            
            
            
        #PATIENT/TUMOR-LEVEL DATA
        res_list = split(res_eventlvl, res_eventlvl$tumor_ID)
        res_patientlvl = data.frame(row.names = names(res_list)) ; {
            
            #Save the claean patient ID
            res_patientlvl[names(res_list),"tumor_ID_raw"] = sapply(res_list, function(x) unique(x$tumor_ID_raw))
            res_patientlvl[names(res_list),"regions_ID"] = sapply(res_list, function(x) unique(x$regions))
            
    
            #DETECTION
            #All variants detected
            detection_patient_mr = sapply(res_list, function(x) ifelse(all(x$mr_detection == 1),1,0))
            detection_patient_pool = sapply(res_list, function(x) ifelse(all(x$pool_detection == 1),1,0))
                res_patientlvl[names(res_list),"mr_detection"] = detection_patient_mr
                res_patientlvl[names(res_list),"pool_detection"] = detection_patient_pool
                
            #All drivers detected
            detectiondriver_patient_mr = sapply(res_list, function(x) ifelse(all(x$mr_detection[which(x$Gene %in% driver_events_rcc)] == 1),1,0))
            detectiondriver_patient_pool = sapply(res_list, function(x) ifelse(all(x$pool_detection[which(x$Gene %in% driver_events_rcc)] == 1),1,0))
                res_patientlvl[names(res_list),"mr_detectiondriver"] = detectiondriver_patient_mr
                res_patientlvl[names(res_list),"pool_detectiondriver"] = detectiondriver_patient_pool
                
            #Proportion of variants detected per tumor  
            detectionpertum_patient_mr = sapply(res_list, function(x){sum(x$mr_detection)/length(x$mr_detection)})
            detectionpertum_patient_pool = sapply(res_list, function(x){sum(x$pool_detection)/length(x$pool_detection)})
                res_patientlvl[names(res_list),"mr_detectionpertum"] = detectionpertum_patient_mr
                res_patientlvl[names(res_list),"pool_detectionpertum"] = detectionpertum_patient_pool

                
            #CLONALITY
            #At least one variant misclassified
            clonality_patient_mr = sapply(res_list, function(x){
                muts_sampled = which(x$mr_detection == 1)
                correct_clonality = x$mr_clonality[muts_sampled] == x$true_clonality[muts_sampled]
                ifelse(all(correct_clonality),1,0)})
            clonality_patient_pool = sapply(res_list, function(x){
                muts_sampled = which(x$pool_detection == 1)
                correct_clonality = x$pool_clonality[muts_sampled] == x$true_clonality[muts_sampled]
                ifelse(all(correct_clonality),1,0)})
                res_patientlvl[names(res_list),"mr_clonality"] = clonality_patient_mr
                res_patientlvl[names(res_list),"pool_clonality"] = clonality_patient_pool
                
                #At least one driver misclassified
                clonalitydriver_patient_mr = sapply(res_list, function(x){
                  muts_sampled = which(x$mr_detection == 1 & x$Gene %in% driver_events_rcc)
                  correct_clonalitydriver = x$mr_clonality[muts_sampled] == x$true_clonality[muts_sampled]
                  ifelse(all(correct_clonalitydriver),1,0)})
                clonalitydriver_patient_pool = sapply(res_list, function(x){
                  muts_sampled = which(x$pool_detection == 1 & x$Gene %in% driver_events_rcc)
                  correct_clonalitydriver = x$pool_clonality[muts_sampled] == x$true_clonality[muts_sampled]
                  ifelse(all(correct_clonalitydriver),1,0)})
                res_patientlvl[names(res_list),"mr_clonalitydriver"] = clonalitydriver_patient_mr
                res_patientlvl[names(res_list),"pool_clonalitydriver"] = clonalitydriver_patient_pool
                
                #Proportion of variants misclassified per tumor
                clonalitypertum_patient_mr = sapply(res_list, function(x){
                  muts_sampled = which(x$mr_detection == 1)
                  correct_prop = sum(x$mr_clonality[muts_sampled] == x$true_clonality[muts_sampled])/length(muts_sampled)})
                clonalitypertum_patient_pool = sapply(res_list, function(x){
                  muts_sampled = which(x$pool_detection == 1)
                  correct_prop = sum(x$mr_clonality[muts_sampled] == x$true_clonality[muts_sampled])/length(muts_sampled)})
                  res_patientlvl[names(res_list),"mr_clonalitypertum"] = clonalitypertum_patient_mr
                  res_patientlvl[names(res_list),"pool_clonalitypertum"] = clonalitypertum_patient_pool
                
            #TARGETABLE MUTATION DETECTION
            evidence_levels = paste(evidence_levels, collapse="|") #Levels of evidence included: 1 through 3
            targets_patient_mr = sapply(res_list, function(x){
                muts_sampled = which(x$mr_detection == 1)
                all_targets = grepl(evidence_levels,x$OncoKB_evidencelvl[muts_sampled])
                ifelse(any(all_targets),1,0)})
            targets_patient_pool = sapply(res_list, function(x){
                muts_sampled = which(x$pool_detection == 1)
                all_targets = grepl(evidence_levels,x$OncoKB_evidencelvl[muts_sampled])
                ifelse(any(all_targets),1,0)})
                res_patientlvl[names(res_list),"mr_targets"] = targets_patient_mr
                res_patientlvl[names(res_list),"pool_targets"] = targets_patient_pool

            #MUTATIONAL LOAD (TMB)
            tmb_truth = sapply(res_list,nrow)   
            tmb_mr = sapply(res_list,function(x) length(x$mutation_ID[x$mr_detection == 1]))   
            tmb_pool = sapply(res_list,function(x) length(x$mutation_ID[x$pool_detection == 1]))   
                res_patientlvl[names(res_list),"true_tmb"] = tmb_truth
                res_patientlvl[names(res_list),"mr_tmb"] = tmb_mr
                res_patientlvl[names(res_list),"pool_tmb"] = tmb_pool
            
                
            #ADDITIONAL OUTCOMES. RENAL CELL CARCINOMA ONLY        
            if(cohort == "RCC"){
                #RCC RISK CATEGORIZATION
                riskrcc_truth = true_risk[res_patientlvl$tumor_ID_raw,]
                riskrcc_mr = as.data.frame(t(sapply(res_list, function(x){genomicRiskRCC(x$Gene[x$mr_detection == 1])})),stringsAsFactors = F)
                riskrcc_pool = as.data.frame(t(sapply(res_list, function(x){genomicRiskRCC(x$Gene[x$pool_detection == 1])})),stringsAsFactors = F)
                    res_patientlvl[,paste0("true_riskrcc_",colnames(riskrcc_truth))] = riskrcc_truth
                    res_patientlvl[names(res_list),"mr_riskrcc"] = ifelse(riskrcc_mr$risk_points == riskrcc_truth$risk_points,1,0)
                    res_patientlvl[names(res_list),"pool_riskrcc"] = ifelse(riskrcc_pool$risk_points == riskrcc_truth$risk_points,1,0)

                #MOLECULAR SUBTYPE ASSIGNMENT
                res_patientlvl$true_subtypercc = true_subtypes[res_patientlvl$tumor_ID_raw]
                    res_patientlvl[names(res_list),"mr_subtypercc"] = sapply(res_list, function(x) unique(x$mr_subtype))
                    res_patientlvl[names(res_list),"pool_subtypercc"] = sapply(res_list, function(x) unique(x$pool_subtype))
                    res_patientlvl$mr_subtypercc = ifelse(res_patientlvl$mr_subtypercc == res_patientlvl$true_subtypercc,1,0)
                    res_patientlvl$pool_subtypercc = ifelse(res_patientlvl$pool_subtypercc == res_patientlvl$true_subtypercc,1,0)
            }
        }
    

        
        #OUTCOME CALCULATION
        detected_mr = which(res_eventlvl$mr_detection == 1) ; 
        detected_pool = which(res_eventlvl$pool_detection == 1)
        
        
  outcomes = list(
          
        #DROPOUT RATE: proportion of mutations missed
        dropout_eventlvl_mr =  1 - mean(res_eventlvl$mr_detection),
        dropout_eventlvl_pool = 1 - mean(res_eventlvl$pool_detection),
        
        #Proportion of patients with at least 1 mutation missed
        dropout_patientlvl_mr =  1 - mean(res_patientlvl$mr_detection),
        dropout_patientlvl_pool = 1 - mean(res_patientlvl$pool_detection),
        
        #Proportion of patients with at least 1 mutation missed
        dropoutdriver_patientlvl_mr =  1 - mean(res_patientlvl$mr_detectiondriver),
        dropoutdriver_patientlvl_pool = 1 - mean(res_patientlvl$pool_detectiondriver),
        
        #Average dropout across tumors
        dropoutmean_patientlvl_mr =  1 - mean(res_patientlvl$mr_detectionpertum),
        dropoutmean_patientlvl_pool = 1 - mean(res_patientlvl$pool_detectionpertum),
            
            
        #CLONALITY ERROR RATE: proportion of incorrect clonality calls 
        #Event level: out of the events detected (all calls - correct calls)
        clonalityerror_eventlvl_mr = 1 - mean(ifelse(res_eventlvl$mr_clonality[detected_mr] == res_eventlvl$true_clonality[detected_mr],1,0)),
        clonalityerror_eventlvl_pool = 1 - mean(ifelse(res_eventlvl$pool_clonality[detected_pool] == res_eventlvl$true_clonality[detected_pool],1,0)),
            
        #Proportion of patients with at least 1 mutation misclassified in terms of clonality
        clonalityerror_patientlvl_mr = 1 - mean(res_patientlvl$mr_clonality) ,
        clonalityerror_patientlvl_pool = 1 - mean(res_patientlvl$pool_clonality),
        
        #Proportion of patients with at least 1 driver misclassified in terms of clonality
        clonalityerrordriver_patientlvl_mr = 1 - mean(res_patientlvl$mr_clonalitydriver),
        clonalityerrordriver_patientlvl_pool = 1 - mean(res_patientlvl$pool_clonalitydriver) ,
        
        #Average clonality error per tumor
        clonalityerrormean_patientlvl_mr = 1 - mean(res_patientlvl$mr_clonalitypertum),
        clonalityerrormean_patientlvl_pool = 1 - mean(res_patientlvl$pool_clonalitypertum),
            
       
        #DETECTION OF TARGETABLE MUTATIONS: number of patients with at least one targetable mutation detected (with evidence equal to evidence_levels)
        targets_mr = mean(res_patientlvl$mr_targets),
        targets_pool = mean(res_patientlvl$pool_targets),
            
        
        #TMB UNDERESTIMATION (true - MR/POOL)/true
        tmbunder_mr = mean((res_patientlvl$true_tmb - res_patientlvl$mr_tmb)/res_patientlvl$true_tmb),
        tmbunder_pool = mean((res_patientlvl$true_tmb - res_patientlvl$pool_tmb)/res_patientlvl$true_tmb)
        


        )  
  
        #COSTS
  cost_sampleprep = costs$cost_prep_persample*(1+region_num) # sample prepping cost is the same for MR and POOL approaches
        outcomes$cost_mr =  costs$cost_overhead + cost_sampleprep + costs$cost_seq*(1+region_num) # 1 normal sample + number of tumor samples
        outcomes$cost_pool = costs$cost_overhead + cost_sampleprep + costs$cost_seq*2# 1 normal sample + 1 pool
        #Cost-effectiveness measures (CE = outcome/cost)
        outcomes$cea_detection_mr = (1-outcomes$dropout_eventlvl_mr)/outcomes$cost_mr
        outcomes$cea_detection_pool = (1-outcomes$dropout_eventlvl_pool)/outcomes$cost_pool
        outcomes$cea_clonality_mr = (1-outcomes$clonalityerror_eventlvl_mr)/outcomes$cost_mr
        outcomes$cea_clonality_pool = (1-outcomes$clonalityerror_eventlvl_pool)/outcomes$cost_pool
  

        # ADDITIONAL RCC-SPECIFIC OUTCOMES
        if(cohort == "RCC"){
            # RISK ERROR
            outcomes$rccriskerror_mr = 1 - mean(res_patientlvl$mr_riskrcc,na.rm = T) 
            outcomes$rccriskerror_pool = 1 - mean(res_patientlvl$pool_riskrcc,na.rm = T) 
            
            #CORRECT SUBTYPE ASSIGNMENT
            outcomes$rccsubtype_mr = mean(res_patientlvl$mr_subtypercc)
            outcomes$rccsubtype_pool = mean(res_patientlvl$pool_subtypercc)}

        outcomes = unlist(outcomes)
        
      
      
    

        

        #CLONALITY CLASSIFIER MEASURES FOR THE DNA POOL
        #Select events that were detected in the pool for classifier  analysis
        events_sampled = res_eventlvl[detected_pool,]
            events_clonal = events_sampled$pool_ccf[events_sampled$true_clonality == "clonal"]
            events_subclonal = events_sampled$pool_ccf[events_sampled$true_clonality == "subclonal"]
        
        #Curves estimates
        #ROC curve, area under the curve
        roc = roc.curve(scores.class0 = events_clonal, scores.class1 = events_subclonal, curve = T)
        #Precision-Recall curve, area under the curve
        pr = pr.curve(scores.class0 = events_clonal, scores.class1 = events_subclonal, curve = T)
    
        #Confusion matrix (numerical) and classifier performance measures 
        classifiers = calculateClassifiers(test = events_sampled$pool_clonality, true = events_sampled$true_clonality, res = "classifiers")
        confusion = calculateClassifiers(test = events_sampled$pool_clonality, true = events_sampled$true_clonality, res = "confusion")
        
        

        
    
        #SAVE RESULTS IN A DATAFRAME (bootstrap) and a LIST (details)
        all_results = list('outcomes' = outcomes,
                           'res_eventlvl' = res_eventlvl,
                           'res_patientlvl' = res_patientlvl,
                           'res_list' = res_list,
                           'classifiers' = classifiers,
                           'confusion_matrix' = confusion)
            details[[sim_id]] = all_results
            
        bootstrap[sim_id,"ID_simulation"] = sim_id
            bootstrap[sim_id,"iteration"] = iteration
            bootstrap[sim_id,"num_regions"] = region_num
            bootstrap[sim_id,"n_tumors"] = n_tumors
            bootstrap[sim_id,"n_events"] = nrow(res_eventlvl)
            bootstrap[sim_id,names(outcomes)] = outcomes
    
        print(sim_id)
        
    
    }
    }




##########################################
###### RESULTS - SUMMARY AND EXPORT ######
##########################################

# Save a dataframe of outcomes (per bootstrapping iteration)
df_path = paste0(output_path,"simresults_outcomes.csv")
write.csv(bootstrap,file = df_path, row.names = F)        

# Rdata object with detailed results of outcomes (per bootstrapping iteration)
details_path = paste0(output_path,"simresults_detailed.Rdata")
save(details, file = details_path, compress = F)




#################################
###### GRAPHICS AESTHETICS ######
#################################

#New labels for outcome and classifier measures
#Names of outcomes  
clean_outcome_names = c("Dropout rate (event-level)","Dropout rate (patient-level)",
                        "Driver dropout rate (patient-level)","Average dropout (patient-level)",
                        "Clonality error (event-level)","Clonality error rate (patient-level)",
                        "Driver clonality error rate (patient-level)","Average clonality error (patient-level)",
                        "Targetable alteration detection (%)","Underestimation of tumor mutational burden (%)",
                        "Absolute cost (in USD)", "Cost-effectiveness (detection)", "Cost-effectiveness (clonality)")
    if(cohort == "RCC"){clean_outcome_names = c(clean_outcome_names, "Risk misattribution rate","Subtype assignment")}
    outcome_measures = unique(gsub("(_mr|_pool)","",names(outcomes)))
    names(outcome_measures) = clean_outcome_names

#Names of classifier measures    
clean_classifier_names = c("AUC of the ROC curve",
                           "AUC of the Precision-Recall (PR) curve\n(piecewise function integration)",
                           "AUC of the Precision-Recall (PR) curve\n(Davis & Goadrich interpolation)",
                           "Sensitivity","Specificity",
                           "False-negative rate (miss-rate)", "False-positive rate (fallout)",
                           "Positive-predictive value (precision)","Negative-predictive value",
                           "False-discovery rate (FDR)","False-omission rate",
                           "Informedness (Youden's J statistic)", "Markedness (deltaP)",
                           "Accuracy","Balanced accuracy",
                           "F1 score", "Illusion of clonality",
                           "Cohen's Kappa statistic", "Fowlkes-Mallows",
                           "Jaccard index", "Jaccard distance",
                           "MCC")
    names(clean_classifier_names) = names(classifiers)
    
# Set up graphic aesthetics and clean labels
color_palette = c("#CF5743","#3985D9","#FFFBC9","#78B96E","#00A00E","#0098BD","#AB355C","#FF8346")
colors_clonality = c('clonal' = '#006E89', 'subclonal' = '#F4A554')
themefull = theme_classic() +
            theme(axis.text.x = element_text(size=14, color = "black"),
                  axis.text.y = element_text(size=14,color="black"),
                  axis.title = element_text(size=14,color="black", face="bold"),
                  plot.title = element_text(size=15.5,color="black", face="bold",hjust=0.5),
                  plot.subtitle = element_text(size=12,color="black", face="bold",hjust=0.5),
                  
                  legend.position = "bottom",
                  legend.title = element_text(size=14,color="black", face="bold"),
                  legend.text = element_text(size=12, color = "black"))
    
    



#############################################
###### PLOTTING RESULTS OF SIMULATIONS ######
#############################################   

#Results per number of regions profiled
results_perregion = lapply(split(bootstrap, bootstrap$num_regions), function(x){
  iter_data = x[,5:ncol(x)]
  mean_95ci = sapply(x[,5:ncol(x)], compute_mean95ci)
  mean_95ci = as.data.frame(t(mean_95ci),stringsAsFactors = F)
  return(mean_95ci)})

# OUTCOMES: Plotting outcomes of multiregional and pooled approaches to NGS, barplots
          # Selects some outcomes to be displayed relative to single-region assessment

relative_outcomes = outcome_measures[outcome_measures %in% relative_outcomes]
    outcomes_path = paste0(output_path,"plots_outcomes/")
    outcome_plots = list()
    
    system(paste("mkdir",outcomes_path))
    for(o in outcome_measures){

      
      #Outcome summary per regions
      res_plot = data.frame(row.names = names(results_perregion))
      res_out = lapply(results_perregion, function(x){
          profiling_results = x[grepl(o,rownames(x)),]
          profiling_results$assay = toupper(gsub(paste0(o,"_"),"",rownames(profiling_results)))
          profiling_results = melt(profiling_results, id.vars = "assay")})
          means = as.data.frame(t(sapply(res_out, function(x){
              out = x$value[x$variable == "mean"]
              names(out) = x$assay[x$variable == "mean"]
              return(out)})), stringsAsFactors = F)
              means$region_num = factor(rownames(means),levels = 1:max_region_num)
              means = melt(means,id.vars = "region_num")
              
          singleregion_mean = which(means$region_num == 1 & means$variable == "MR")
              means$value_rela = means$value/means$value[singleregion_mean]
              
              low95 = as.data.frame(t(sapply(res_out, function(x){
                  out = x$value[x$variable == "low95"]
                  names(out) = x$assay[x$variable == "low95"]
                  return(out)})), stringsAsFactors = F)
              low95$region_num = factor(rownames(low95),levels = 1:max_region_num)
              low95 = melt(low95,id.vars = "region_num")
              low95$value_relative = low95$value/low95$value[singleregion_mean]
              
              
              up95 = as.data.frame(t(sapply(res_out, function(x){
                out = x$value[x$variable == "up95"]
                names(out) = x$assay[x$variable == "up95"]
                return(out)})), stringsAsFactors = F)
              up95$region_num = factor(rownames(up95),levels = 1:max_region_num)
              up95 = melt(up95,id.vars = "region_num")
              up95$value_relative = up95$value/up95$value[singleregion_mean]
              
    
          cis = as.data.frame(cbind(up95, low95[,3:4]),stringsAsFactors = F)
              colnames(cis) = c("region_num","assay","up95","up95_rela","low95","low95_rela")
          
      #GGPLOT OBJECTS AND PRINTING THE PLOT AS PDFs
      outcome_plots[[o]] = ggplot() + 
          geom_errorbar(data = cis, aes(x = region_num, ymax =up95, ymin = low95, group = assay),
                        color = "#000000", position = "dodge", size = 0.5,width = 0.7) +
        
          geom_bar(data = means,aes(x = region_num, y = value, group = variable, fill = variable), stat = "identity", 
                   position = "dodge", color = "black", width = 0.75) +
          ggtitle(names(outcome_measures)[outcome_measures == o]) +
          scale_fill_manual("Type of profiling",values = color_palette) +
          ylab(o) + xlab("Number of regions profiled") +
          theme_classic() +
          themefull

      
          path_plot = paste0(outcomes_path,"barplot_",o,".pdf") #Path for plot
          pdf(path_plot,7,7)
          print(outcome_plots[[o]])
          dev.off()
      
     
          
      #EXPRESS SOME OUTCOMES RELATIVE TO SINGLE-REGION RESULTS    
      if(o %in% relative_outcomes){
        ro = paste0("r_",o)
        outcome_plots[[ro]] = ggplot() + 
            geom_errorbar(data = cis, aes(x = region_num, ymax =up95_rela, ymin = low95_rela, group = assay),
                          color = "#000000", position = "dodge", size = 0.5,width = 0.7) +
            
            geom_bar(data = means,aes(x = region_num, y = value_rela, group = variable, fill = variable), stat = "identity", 
                     position = "dodge", color = "black", width = 0.75) +
            ggtitle(names(outcome_measures)[outcome_measures == o],
                    subtitle = "(relative to single region profiling)") +
            ylab(paste0(o," (vs 1-region)")) +
            theme_classic() +
            themefull +
            scale_fill_manual("Type of profiling",values = color_palette)
          
            path_plot = paste0(outcomes_path,"barplot_",o,"_relat.pdf") #Path for plot
            pdf(path_plot,7,7)
            print(outcome_plots[[ro]])
            dev.off()
      }
      
    }
    #Combined plot
    combplot_outcome = plot_grid(plotlist = outcome_plots, nrow = 4, ncol = 5, labels = LETTERS[1:15])
        path_plot = paste0(output_path,"outcomes_combined.pdf") #Path for plot
        pdf(path_plot,25,20)
        print(combplot_outcome)
        dev.off()
    
        

# CLASSIFIER ACCURACY: Plotting classifier accuracy measures, boxplots
# USES MR_threshold and POOL_threshold  
#Create a subdirectory to save the region_specific results
classifiers_path = paste0(output_path,"plots_classifiers/")
    system(paste("mkdir",classifiers_path))
    classifier_plots = list()
    
    classifier_data = sapply(details, function(x) x[['classifiers']])
        classifier_data = as.data.frame(t(classifier_data),stringsAsFactors = F)
        classifier_data$regions_pooled = sapply(strsplit(rownames(classifier_data),split="_"),function(x)x[2])
            classifier_data$regions_pooled = factor(gsub("r","",classifier_data$regions_pooled), levels = 1:max_region_num)
        classifier_data$iteration = sapply(strsplit(rownames(classifier_data),split="_"),function(x)x[1])
        classifier_data$sim_id = rownames(classifier_data)
        
    
    for(cl in names(classifiers)){
        classifier_toplot = classifier_data[,cl]
        classifier_label = clean_classifier_names[cl]
        
        summarized_classifier_data = sapply(split(classifier_data[,cl], classifier_data$regions_pooled), compute_mean95ci)
            summarized_classifier_data = as.data.frame(t(summarized_classifier_data),stringsAsFactors = F)
            summarized_classifier_data$regions_pooled = factor(rownames(summarized_classifier_data), levels = rownames(summarized_classifier_data))

        classifier_plots[[cl]] = ggplot(summarized_classifier_data, 
                                        aes(x=regions_pooled, y=mean, ymax = up95, ymin = low95)) +
            geom_errorbar(color = "#3985D9", size = 1.05, width = 0.3) +
            geom_point(size = 6, color = "#3985D9", shape = 18) +
            xlab("Regions pooled") + ylab("Value") +
            ggtitle(classifier_label) +
            themefull
            path_plot = paste0(classifiers_path,"boxplot_",cl,".pdf") #Path for plot
            pdf(path_plot,7,7)
            (print(classifier_plots[[cl]])) #### WARNING_PLOTS
            dev.off()}
        
    combplot_classifier = plot_grid(plotlist = classifier_plots, nrow = 5, ncol = 5, labels = LETTERS[1:length(classifiers)])
        path_plot = paste0(output_path,"classifiers_combined.pdf") #Path for plot
        pdf(path_plot,25,25)
        print(combplot_classifier) 
        dev.off()
        
  
        
        
        
        
        
        
#######################################
###### PLOTTING DETAILED RESULTS ######
#######################################     
          
# Plotting detailed results in a separate subdirectory
# This is, estimates for all the simulated DNA pools of X number of regions (where X is 'spec_region_num')        
details_path = paste0(output_path,"plots_detailed_pool",spec_region_num,"regions/")
    system(paste("mkdir",details_path))
    
    #Obtain the event-level results over the iterations with the specified number of regions (defaults to 4)
    # Defaults to all iterations but allows any number
    detail_res = paste0("i",1:iter_num,"_r",spec_region_num)
        detail_res = details[detail_res]
        detail_allevents = data.frame() ; detail_alltumors = data.frame()
        for(i in 1:length(detail_res)){
            eventdata = detail_res[[i]]$res_eventlvl
                detail_allevents = as.data.frame(rbind(detail_allevents,eventdata),stringsAsFactors = F)
            tumordata = detail_res[[i]]$res_patientlvl
                detail_alltumors = as.data.frame(rbind(detail_alltumors,tumordata),stringsAsFactors = F)
        }
        
    #Detection - genes dropped in pool
    genes_dropped = which(detail_allevents$pool_detection == 0) ; {  
        genes_dropped = detail_allevents[genes_dropped,"Gene"] 
        genes_dropped = as.data.frame(sort(table(genes_dropped),dec=T), stringsAsFactors = F)
        genes_dropped$Prop = genes_dropped$Freq/sum(genes_dropped$Freq)
        genes_dropped_top10 = genes_dropped[1:10,]
        genes_dropped_top10$genes_dropped = factor(genes_dropped_top10$genes_dropped, levels = rev(genes_dropped_top10$genes_dropped))
        
      plotdet_dropout = ggplot(genes_dropped_top10, aes(x=genes_dropped,y = Prop)) +
        geom_bar(stat="identity", fill = "#015980", color = "#000000", width = 0.7) +
        xlab("") + ylab("Proportion (out of n-dropped)") + ggtitle("Dropout, top 10 genes") +
        themefull +
        coord_flip() +
        theme(axis.text.y = element_text(face = "bold"))}

    #Clonal error - genes misclassified     
    genes_misclassified = which(detail_allevents$pool_clonality != detail_allevents$true_clonality) ; {
        genes_misclassified = detail_allevents[genes_misclassified,"Gene"] 
        genes_misclassified = as.data.frame(sort(table(genes_misclassified),dec=T), stringsAsFactors = F)
        genes_misclassified$Prop = genes_misclassified$Freq/sum(genes_misclassified$Freq)
        genes_misclassified_top10 = genes_misclassified[1:10,]
        genes_misclassified_top10$genes_misclassified = factor(genes_misclassified_top10$genes_misclassified, levels = rev(genes_misclassified_top10$genes_misclassified))
      
    plotdet_clonalerror = ggplot(genes_misclassified_top10, aes(x=genes_misclassified,y = Prop)) +
        geom_bar(stat="identity", fill = "#015980", color = "#000000", width = 0.7) +
        xlab("") + ylab("Proportion (out of n-misclassified)") + ggtitle("Clonal error, top 10 genes") +
        themefull +
        coord_flip() +
        theme(axis.text.y = element_text(face = "bold"))}
    
        #Combine both plots and print in a single PDF
        plots_combined_eventlvl = plot_grid(plotdet_dropout, plotdet_clonalerror, nrow = 1, ncol = 2, labels = c("A","B"))
            plotname = paste0(details_path,"outcomes_eventlvl_ccfpool.pdf")
            pdf(plotname,10,8)
            print(plots_combined_eventlvl)
            dev.off()

    #Targetable mutation detection
    targets_dropped = which(detail_allevents$pool_detection == 0 & grepl(evidence_levels,detail_allevents$OncoKB_evidencelvl)) ; {  
      targets_dropped = detail_allevents[targets_dropped,] 
      pieplot_targets = plot_pie(targets_dropped$Gene, colpal = c("#55C0D5","#E9972A",color_palette),plot_title = "Targetable mutation dropout")
      plotname = paste0(details_path,"outcomes_translational.pdf")}
    
    if(cohort == "RCC"){
      
        #Risk misclassification
        riskcond_missed = which(detail_alltumors$pool_riskrcc == 0)
            riskcond_missed = detail_alltumors[riskcond_missed,c("true_riskrcc_bap1.tp53","true_riskrcc_multidriver","true_riskrcc_pbrm1wt")] 
            riskcond_missed$id = rownames(riskcond_missed)
            riskcond_missed = melt(riskcond_missed,id.vars = "id")
            riskcond_missed = riskcond_missed[which(riskcond_missed$value != 0),]
            riskcond_missed$variable = toupper(gsub("^true_riskrcc_","",riskcond_missed$variable))
            pieplot_riskcond = plot_pie(riskcond_missed$variable, colpal = c("#00A969","#38499D","#474747",color_palette),plot_title = "Risk misclassification features")
        
        #Subtype assignment
        subtypes_missed = which(detail_alltumors$pool_subtypercc == 0)   
            subtypes_missed = detail_alltumors[subtypes_missed,] 
            pieplot_subtypes = plot_pie(subtypes_missed$true_subtypercc, colpal = color_palette,plot_title = "Subtype misclassification")
        
        #Combining all plots and printing    
        pieplots_comb = plot_grid(pieplot_targets, pieplot_riskcond, pieplot_subtypes, labels = c("A","B","C"), nrow = 1, ncol = 3)
            pdf(plotname,15,5)
            print(pieplots_comb)
            dev.off()
        
    }else{
        pdf(plotname,5,5)
        print(pieplot_targets)
        dev.off()
        }
            
    #CLASSIFIER ACCURACY MEASURES: selects only detected events for the additional classifier results
    detail_allevents = detail_allevents[which(detail_allevents$pool_detection == 1),] ; {
        clonal_muts = detail_allevents$pool_ccf[detail_allevents$pool_clonality == "clonal"]
        subclonal_muts = detail_allevents$pool_ccf[detail_allevents$pool_clonality == "subclonal"]
    
        # ROC and Precision-recall (PR) curves - PRROC package
        roc = roc.curve(scores.class0 = clonal_muts, scores.class1 = subclonal_muts, curve = T)
        pr = pr.curve(scores.class0 = clonal_muts, scores.class1 = subclonal_muts, curve = T)
        plotname = paste0(details_path,"curves_ROC_PR.pdf")
            pdf(plotname,5,4)
            par(mfrow = c(2,1))
            plot(roc) ; plot(pr)
            dev.off()
            
            
        #Additional plots to display the accuracy of the CCF in the pool to classify clonality 
        otherplots = list()
          
          #Histogram of CCF values in the simulated pool, broken down by the true clonality status
          otherplots$histo_clonal = ggplot(detail_allevents, aes(x = pool_ccf, y=..density.., fill = true_clonality, group = true_clonality)) + 
              geom_histogram(color = "#000000", binwidth = .05) +
              themefull + xlab("CCF tumor pool") + ylab("Density") +
              scale_fill_manual('True clonality', values = colors_clonality)
          
          #Violin plot of CCF values in the simulated pool, broken down by the true clonality status
          otherplots$violin_clonal = ggplot(detail_allevents, aes(y = pool_ccf, x=true_clonality, color = true_clonality)) + 
              geom_sina(size = 3, alpha = 0.5) +
              themefull + ylab("CCF tumor pool") + xlab("True clonality") +
              scale_color_manual('True clonality', values = colors_clonality)
    
          #Optimal CCF cutoff in the simulated pool, given the specified number of regions and MR_threshold
          allcutoffs = seq(0.52,0.99,0.0001)
          accuracy_res = sapply(allcutoffs, function(ccf_threshold){
                  #New clonality calls with specified cutoff
                  clonality_calls = ifelse(detail_allevents$pool_ccf >= ccf_threshold,"clonal","subclonal") #Makes new clonality calls to 
                  classifiers_res = calculateClassifiers(test = clonality_calls, true = detail_allevents$true_clonality, res = "classifiers")
                  res = classifiers_res[c('accuracy',"f1","mcc")] ; names(res) = c("Accuracy","F1 score","MCC")
                  return(res)})
              accuracy_res = as.data.frame(t(accuracy_res),stringsAsFactors = F)
              accuracy_res$ccf_threshold = allcutoffs
              accuracy_res = melt(accuracy_res, id.vars = "ccf_threshold")
    
          otherplots$plot_classifiers = ggplot(accuracy_res, aes(x=ccf_threshold, y=value, group = variable, color = variable)) +
              geom_line(size = 1.2, alpha = 0.9) +
              scale_color_manual("Measure",values = color_palette[c(1,2,4)]) +
              xlab(paste("CCF cutoff performance in pools of",spec_region_num,"regions")) + ylab("Value") +
              themefull
          
          #Printing single plots
          sapply(1:length(otherplots), function(p){
              label = names(otherplots)[p]
              plotname = paste0(details_path,label,"_ccfpool.pdf")
              pdf(plotname,7,7)
              print(otherplots[[p]])
              dev.off()})
          
        
        }
                
        
        










