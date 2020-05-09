#rm(list=ls())
#setwd("/Users/dinatalr/Documents/MSKCC/rcc_translational/rcc_ith/analyses/scripts/final_code/")


#############################################################################
#       MULTIREGIONAL DNA POOLING OVERCOMES INTRATUMORAL HETEROGENEITY      #
#                   IN CLINICAL SEQUENCING OF TUMORS                        #
#                                                                           #
#                            DiNatale RG et al.                             #
#############################################################################

# This R script uses multi-regional tumor mutation data (in MAF format) from the single-region 
# and tumor DNA pools sequenced as part of the study and compares the outcomes of these approaches 
# against a conventional 'multi-regional assessment'

#The following outputs will be produced (panels of supplementary figure 2):
#1- An oncoprint showing the concordance in clonality calls for each sample (correct, incorrect, missed) (SF2A)
#2- The comparison of observed versus expected CCF and purity values in the DNA pools (SF2B)
#3- A summary of the outcomes for each type of sequencing approach (SF2C)
#4- The PR and ROC curves of the CCF in the pool to classify clonality (SF2D)
#5- The performance of each CCF cutoff in these pools (containing 6 regions) to classify clonality (SF2E)



#Install/load packages and basic functions
pkg_touse = c("readxl","dplyr","reshape2","ggplot2","cowplot",
              "PRROC","ggforce","grDevices","RColorBrewer") ; {
                pkg_toinstall = pkg_touse[!pkg_touse %in% rownames(installed.packages())]
                if(length(pkg_toinstall)>0){install.packages(pkg_toinstall)}
                lapply(pkg_touse, require, character.only = TRUE)}


MR_threshold = 0.5
POOL_threshold = 0.75
detection_limit = 0.02

# Basic functions 
{
  #Call an event clonal in multiregional assessment
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
    
    #Sensitivity, specificity and derived measures
    classifiers$sensitivity = tp/(tp+fn) # proportion of clonal events correctly classified
    classifiers$specificity = tn/(tn+fp) # proportion of subclonal events correctly classified
    classifiers$missrate_fnr = fn/(fn+tp) # (1-sens) proportion of clonal missclasified as subclonal 
    classifiers$fallout_fpr = fp/(fp+tn) # (1-spec) proportion of subclonal missclasified as clonal
    classifiers$informedness_bm = classifiers$sensitivity + classifiers$specificity -1
    
    # PPV, NPV and derived measures
    classifiers$ppv_precison = tp/(tp+fp) # proportion of real clonal events out of all clonal calls
    classifiers$npv = tn/(tn+fn) # proportion of real subclonal events out of all subclonal calls
    classifiers$fdr = fp/(fp+tp) # (1-PPV) proportion of incorrect out of all the clonal calls
    classifiers$fomissionr = fn/(fn+tn) # (1-NPV) proportion of incorrect out of all the subclonal calls
    classifiers$markedness = classifiers$ppv_precison + classifiers$npv - 1
    
    #Conventional measures
    classifiers$csi_threatscore =    tp/(tp+fp+fn) # true clonal out of all possible positives (real and not real)
    classifiers$accuracy = (tp+tn)/sum(confusion)
    classifiers$accuracy_bal = (classifiers$sensitivity+classifiers$specificity)/2
    classifiers$f1 = 2 * (classifiers$ppv_precison * classifiers$sensitivity)/
      (classifiers$ppv_precison + classifiers$sensitivity)
    
    #Cohen's Kappa coefficient ( Po - Pe / 1 - Pe)
    prandom_clonal = ( (tp+fp)/sum(confusion) ) * ( (tp+fn)/sum(confusion) )
    prandom_subclonal = ( (tn+fn)/sum(confusion) ) * ( (tn+fp)/sum(confusion) )
    pexpect = prandom_clonal + prandom_subclonal
    classifiers$kappa_cohen = (classifiers$accuracy - pexpect) / (1-pexpect)
    
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

#Output directory
output_path = paste0("results_experiment","/")
system(paste("mkdir",output_path))



#Read in the data and separate the multiregional from the pool samples
maf = read.csv("ST2_MAF_DNApooling_experiment.csv",stringsAsFactors = F)

mr = maf[maf$Seq_approach == "MR",]
    list.pats = split(mr, mr$patient_ID)   
    list.regs = lapply(list.pats, function(x) split(x, x$Tumor_Sample_Barcode))

pool = maf[maf$Seq_approach == "POOL",]
    pool.pats = split(pool, pool$patient_ID)  


#### DETERMINE EVENT CLONALITY AND EXTRACT INFORMATION FOR EACH EVENT #####
allevents = unique(maf$Mutation_ID)
muts = sapply(1:length(allevents), function(e){
        event = allevents[e] ; names(event) = "event"
        patient = sapply(strsplit(event,split = "_", fixed = T), function(x) x[1]) ; names(patient) = "patient"
        gene = sapply(strsplit(event,split = "_", fixed = T), function(x) x[2]) ; names(gene) = 'gene'
        
        #True clonality and expected CCF values
        mrdata = list.regs[[patient]]
        pooldata = pool.pats[[patient]]
        
            #Obtain the event information in the separate regions and in the pool
            infoevent_mr = sapply(mrdata, function(x){
                eventindex = which(x$Mutation_ID == event)
                    #if event is not present in a region, obtain the purity and set CCF = 0
                    if(length(eventindex) != 0){
                        # In the multiregional assessment, we use the upper bound of the CCF estimate to define clonality 
                        # This is because after manual inspection of the data, there is a VHL mutation (RCC003_T6) where the estimated CCF is 42%, 
                        # we considered this a likely artifact since the purity of the sample is borderline adequate (0.23)
                        info = as.numeric(x[eventindex, c("Purity","CCF_95up")])
                    }else{info = c(unique(x$Purity), 0)} 
                    names(info) = c("Purity","CCF")
                return(info)
            })
                infoevent_mr = as.data.frame(t(infoevent_mr), stringsAsFactors = F)
                #If purity is missing (i.e. <20%), set it to 0.15 by default
                infoevent_mr$Purity[is.na(infoevent_mr$Purity)] = 0.15
                #If CCF is missing (due to low purity), we assume that only the clonal mutations were detected (sets CCF=1)
                infoevent_mr$CCF[is.na(infoevent_mr$CCF)] = 1
                
            info_pool = function(x){
                eventindex = which(pooldata$Mutation_ID == x)
                    #if event is not present in a region, obtain the purity and set CCF = 0
                    if(length(eventindex) != 0){
                      info = as.numeric(pooldata[eventindex, c("Purity","CCF")])
                    }else{info = c(unique(pooldata$Purity), 0)} 
                names(info) = c("Purity","CCF")
                return(info)}
            infoevent_pool = info_pool(event)    
            
        #True clonality assertion (using all 6 regions)
        clonal_true = callClonal_MR(ccfs = infoevent_mr$CCF,cutpoint = MR_threshold) ; names(clonal_true) = "clonal_true"
        
        #Clonality assertion in each region
        clonal_sr = ifelse(infoevent_mr$CCF >= MR_threshold, "clonal","subclonal")
            clonal_sr[which(infoevent_mr$CCF == 0)] = "miss"
            names(clonal_sr) = paste0("clonal_sr",1:length(mrdata))
        ccf_sr = infoevent_mr$CCF ; names(ccf_sr) = paste0("ccf_sr",1:length(mrdata))
        pur_sr = infoevent_mr$CCF ; names(pur_sr) = paste0("pur_sr",1:length(mrdata))
        
        
        # Clonality assertion in the simulated DNA pool
        clonal_simpool = callClonal_POOL(ccfs = infoevent_mr$CCF,pur = infoevent_mr$Purity, 
                                         cutpoint = POOL_threshold, detect_lim = detection_limit)
            pur_simpool = mean(infoevent_mr$Purity) ; names(pur_simpool) = "pur_simpool"
            ccf_simpool = as.numeric(clonal_simpool) ; names(ccf_simpool) = "ccf_simpool"
            clonal_simpool = names(clonal_simpool) ; names(clonal_simpool) = "clonal_simpool"
            
        # Clonality assertion in the real DNA pool
        ccf_pool = infoevent_pool['CCF'] ; names(ccf_pool) = "ccf_realpool"    
        clonal_pool = ifelse(ccf_pool >= POOL_threshold, "clonal","subclonal") 
            clonal_pool = ifelse(ccf_pool == 0,"miss",clonal_pool) ; names(clonal_pool) = "clonal_realpool"
        pur_pool = infoevent_pool['Purity'] ; names(pur_pool) = "pur_realpool"         

        #Combine all data for the analysis
        eventdata = c(event, patient, gene, 
                      ccf_sr, pur_sr, clonal_sr,
                      clonal_true,clonal_pool, clonal_simpool,
                      pur_simpool, pur_pool, ccf_simpool, ccf_pool)
        return(eventdata)})
    muts = as.data.frame(t(muts),stringsAsFactors = F)
    rownames(muts) = muts$event
    muts[,grep("ccf|pur",colnames(muts))] = sapply(muts[,grep("ccf|pur",colnames(muts))], as.numeric)
    muts[,grep("clonal",colnames(muts))] = sapply(muts[,grep("clonal",colnames(muts))], factor, levels = c("miss","subclonal","clonal"))

    
    
    
######################            
###### PLOTTING ######
######################    
    
#Ordering genes for plotting
  order.genes = c("VHL","PBRM1","EP300","ERCC4","KDM6A","PDGFRA","TP53",
                  "TSC1","BAP1","ERBB2","MYCN","NFE2L2","NOTCH1","KDM5C")

#Create a ggplot2 theme and a color palette
colpal = c("#B22D29","#199535","#8396CE","#474543","#C57400","#7900DA")
colors_clonality = c('clonal' = '#006E89', 'subclonal' = '#F4A554')

  
mytheme = theme_classic() +
    theme(axis.text.x = element_text(size=14,color="black"),
          axis.text.y = element_text(size=14,color="black"),
          title = element_text(size=16,color="black",face="bold",hjust = 0.5),
          legend.title = element_text(size=16,color="black", face = "bold"),
          legend.text = element_text(size=14,color="black"),
          strip.text = element_text(size=14,color="black", face = "bold"),
          legend.position = "right")



##### Plot expected vs observed #####
cor.ccf = cor.test(muts$ccf_simpool, muts$ccf_realpool, method = "pearson", exact = F)   
cor.ccf_lab = paste0("Pearson, p<0.001\nrho=",round(cor.ccf$estimate,2))
ccf_obsexp =ggplot(muts) + 
    geom_abline(slope = 1, intercept = 0, linetype="dashed", alpha = 0.6) +
    geom_point(aes(x=ccf_simpool, y=ccf_realpool,color = patient), size = 3, alpha = 0.9) +
    xlab("Simulated pool") + ylab("Real pool") +
    scale_color_manual("Tumor",values = colpal) +
    scale_x_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.2)) +
    scale_y_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.2)) +
    ggtitle("CCF") +
    annotate('text', x = 0.75, y=0.1, label = cor.ccf_lab, size = 5) +
    mytheme
    plotname = paste0(output_path,"observed_expected_ccf.pdf")
    pdf(plotname,6,5)
    print(ccf_obsexp)
    dev.off()


cor.pur = cor.test(muts$pur_simpool, muts$pur_realpool, method = "pearson", exact = F)    
cor.pur_lab = paste0("Pearson, p<0.001\nrho=",round(cor.pur$estimate,2))
pur_obsexp = ggplot(muts) + 
    geom_abline(slope = 1, intercept = 0, linetype="dashed", alpha = 0.6) +
    geom_point(aes(x=pur_simpool, y=pur_realpool,color = patient), size = 3, alpha = .9) +
    scale_color_manual("Tumor",values = colpal) +
    xlab("Simulated pool") + ylab("Real pool") +
    scale_x_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.2)) +
    scale_y_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.2)) +
    ggtitle("Purity") +
    annotate('text', x = 0.75, y=0.1, label = cor.pur_lab, size = 5) +
    mytheme
  
    plotname = paste0(output_path,"observed_expected_purity.pdf")
    pdf(plotname,6,5)
    print(pur_obsexp)
    dev.off()

#####
    
    
    
##### Plot matches clonal/subclonal - gold standard is full multiregional assessment #####
    
conc = muts[,grep("^clonal_",colnames(muts))]
    conc = sapply(conc, function(x){
       res = ifelse(x == conc$clonal_true,"correct","incorrect")
       res[x == "miss"] = "missed"
       return(res)})    
    conc = as.data.frame(conc, stringsAsFactors = F)
    conc = conc[,-grep("true|simpool",colnames(conc))]
    conc$event = muts$event

    #Melt the data frame
    concmelt = melt(conc, id.vars = "event",factorsAsStrings = T)
        concmelt$patient = sapply(strsplit(concmelt$event,split = "_"), function(x)x[1])
        concmelt$gene = sapply(strsplit(concmelt$event,split = "_"), function(x)x[2])
        #Add duplets    
        dups = c("RCC002_VHL_10183838_CCTGGCA_C","RCC003_EP300_41562637_C_T",
                 "RCC004_TP53_7578535_T_A","RCC004_TSC1_135798733_A_G")
            dups = which(concmelt$event %in% dups)
            concmelt$newgene = concmelt$gene
            concmelt$newgene[dups] = paste0(concmelt$newgene[dups],"_2")
            concmelt$gene = factor(concmelt$gene, levels = order.genes)
            concmelt = concmelt[order(concmelt$gene,decreasing = F),]   
            order.newgene = unique(concmelt$newgene)
            concmelt$newgene = factor(concmelt$newgene, levels = rev(order.newgene))
        concmelt$value = factor(concmelt$value, levels = c("missed","incorrect","correct"))
        concmelt$variable = sapply(strsplit(as.character(concmelt$variable),split = "_"), function(x)x[2])
            concmelt$variable = toupper(gsub("^s|real","",concmelt$variable))
            concmelt$variable = factor(concmelt$variable, levels = c(paste0("R",1:6),"POOL"))

    #Oncoprint - note that patients with >1 mutation in each gene will appear as overlapped tiles            
    match_oncoprint = ggplot(concmelt) + 
        geom_tile(aes(x=variable, y=newgene, fill = value), color = "black",size=0.5, width = 0.75, height = 0.75) +
        theme_minimal() +
        xlab("Sample") + ylab("Gene") +
        scale_fill_manual("Clonality call",values = c("#D1D1D1","#CF5743","#00A969")) +
        guides(fill = guide_legend(reverse = TRUE)) +
        ylab("") + xlab("") +
        facet_wrap(~patient,nrow=1) +
        mytheme +
        theme(axis.text.x = element_text(angle=90))
        plotname = paste0(output_path, "oncoprint_clonalcalls.pdf")
        pdf(plotname,11,7)
        print(match_oncoprint)
        dev.off()

    #Top barplot
    match_barplot = ggplot(concmelt) +
        geom_bar(aes(x=variable, fill = value), 
                 position = "fill", color= "black", size = 0.6, width = 0.75) +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_manual("Clonality call",values = c("#D1D1D1","#CF5743","#00A969")) +
        facet_wrap(~patient,nrow=1) +
        ylab("") +
        mytheme
        plotname = paste0(output_path, "barplot_clonalcalls.pdf")
        pdf(plotname,11,2.5)
        print(match_barplot)
        dev.off()
#########        
        
        
##### Outcomes - single region vs pooled sequencing ######
outcomes = melt(conc, id.vars = "event",factorsAsStrings = F)
    outcomes$patient = sapply(strsplit(outcomes$event,split = "_"), function(x) x[1])
    outcomes$sample = paste(outcomes$patient,outcomes$variable,sep="_")
    outcomes$sampletype = ifelse(grepl("pool",outcomes$variable),"Tumor pool","Single region")
    outcomes$detection = ifelse(outcomes$value != "missed",1,0)
    outcomes$clonality = ifelse(outcomes$value == "correct",1,0)
    outcomes$clonality[outcomes$value == "missed"] = NA
    
    #Calculate dropout and clonal error per sample
    outcomes_samps = sapply(split(outcomes, outcomes$sample), 
                            function(x){ 
                              dropout = 1 - (sum(x$detection)/length(x$detection)) 
                              x = x[!is.na(x$clonality),]
                              clonalerror = 1 - (sum(x$clonality)/length(x$clonality))
                              res = c(dropout, clonalerror) ; names(res) = c("dropout","clonalerror")
                              return(res)})
    outcomes_samps = as.data.frame(t(outcomes_samps), stringsAsFactors = F)
        outcomes_samps$sampletype = ifelse(grepl("pool",rownames(outcomes_samps)),"Tumor pool","Single region")
    
    dropout = split(outcomes_samps$dropout,outcomes_samps$sampletype)
        ttest_dropout = t.test(dropout[[1]], dropout[[2]], paired = F, exact = F, correct = F, conf.int = T)
        ttest_dropout_lab = paste0("p=",signif(ttest_dropout$p.value,1))
        dropout = sapply(dropout,compute_mean95ci)    
        dropout =  as.data.frame(t(dropout), stringsAsFactors = F)
        dropout$sampletype = rownames(dropout)
        dropout$low95[dropout$low95<0]=0
    
    clonalerror = split(outcomes_samps$clonalerror,outcomes_samps$sampletype)
        ttest_clonalerror = t.test(clonalerror[[1]], clonalerror[[2]], exact = F, correct = F, conf.int = T)
        ttest_clonalerror_lab = paste0("p=",signif(ttest_clonalerror$p.value,1))
        clonalerror = sapply(clonalerror,compute_mean95ci)    
        clonalerror =  as.data.frame(t(clonalerror), stringsAsFactors = F)
        clonalerror$sampletype = rownames(clonalerror)
        clonalerror$low95[clonalerror$low95<0]=0
        
#Plots
plot_dropout = ggplot(dropout, aes(x = sampletype, fill = sampletype)) +
    geom_errorbar(aes(ymin = low95, ymax = up95), size = 0.7, width = 0.4) +
    geom_bar(stat = "identity", aes(y=mean), color = "#000000") +
    ylab('Average dropout') + xlab("") + 
    scale_fill_manual("", values = c("#3985D9","#CF5743")) +
    scale_y_continuous(limits = c(0,0.3), breaks = seq(0,0.3,.05)) +
    mytheme
    plotname = paste0(output_path, "barplot_dropout.pdf")
    pdf(plotname,6,6)
    print(plot_dropout)
    dev.off()


plot_clonalerror = ggplot(clonalerror, aes(x = sampletype, fill = sampletype)) +
    geom_errorbar(aes(ymin = low95, ymax = up95), size = 0.7, width = 0.4) +
    geom_bar(stat = "identity", aes(y=mean), color = "#000000") +
    ylab('Average clonalerror') + xlab("") + 
    scale_fill_manual("", values = c("#3985D9","#CF5743")) +
    scale_y_continuous(limits = c(0,0.3), breaks = seq(0,0.3,.05)) +
    mytheme
    plotname = paste0(output_path, "barplot_clonalerror.pdf")
    pdf(plotname,6,6)
    print(plot_clonalerror)
    dev.off()

    
#######    
    
    
    
##### Performance of the CCF in the real pool to classify mutations as clonal/subclonal #####
eventspool = muts[which(muts$clonal_realpool != "miss"),]    
    clonalpool =  eventspool$ccf_realpool[which(eventspool$clonal_true == "clonal")]   
    subclonalpool =  eventspool$ccf_realpool[which(eventspool$clonal_true == "subclonal")]   
    
    # ROC and Precision-recall (PR) curves - PRROC package
    roc = roc.curve(scores.class0 = clonalpool, scores.class1 = subclonalpool, curve = T)
    pr = pr.curve(scores.class0 = clonalpool, scores.class1 = subclonalpool, curve = T)
        plotname = paste0(output_path,"curves_ROC_PR.pdf")
        pdf(plotname,5,4)
        par(mfrow = c(2,1))
        plot(roc) ; plot(pr)
        dev.off()    

    #Optimal CCF cutoff in the simulated pool, given the specified number of regions and MR_threshold
    allcutoffs = seq(0.1,1,0.1)
    accuracy_res = sapply(allcutoffs, function(ccf_threshold){
        #New clonality calls with specified cutoff
        clonality_calls = ifelse(eventspool$ccf_realpool >= ccf_threshold,"clonal","subclonal") #Makes new clonality calls to 
        classifiers_res = calculateClassifiers(test = clonality_calls, true = eventspool$clonal_true, res = "classifiers")
        res = classifiers_res[c('accuracy',"mcc")] ; names(res) = c("Accuracy","MCC")
        return(res)})
        accuracy_res = as.data.frame(t(accuracy_res),stringsAsFactors = F)
        accuracy_res$ccf_threshold = allcutoffs
    accuracy_res = melt(accuracy_res, id.vars = "ccf_threshold")
    
    #Plotting the results
    plot_classifiers = ggplot(accuracy_res, aes(x=ccf_threshold, y=value, group = variable, color = variable)) +
        geom_vline(xintercept = POOL_threshold, linetype = "dashed", color = "#000000", alpha = 0.6) +
        geom_line(size = 1.1, alpha = 0.9) +
        geom_point(size = 2, alpha = 0.9) +
        scale_color_manual("Cut-off performance",values = colpal[c(1,2,4)]) +
        xlab("CCF in tumor pool") + ylab("Value") +
        scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.2)) +
        scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2)) +
        mytheme +
        theme(legend.position = "bottom")
        plotname = paste0(output_path,"classifiers_ccfpool.pdf")
        pdf(plotname,6.3,4.5)
        suppressWarnings(print(plot_classifiers)) # warning message: extreme CCF cutoffs can return NA (e.g. 10%, 100%)
        dev.off()
######
    
    
        
        
        
            
