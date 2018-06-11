##Code to remove all rows w/ zero from 5-hmC beta values##

#Load in oxBS-MLE beta data:
load("/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/oxBSMLE_betas.RData")

#Additionally, load in significant probes 5-hmC output files:
load("/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/5hmC/fivehmc_glialadjusted_nopmi_cingulate.RData")
load("/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/5hmC/fivehmc_glialadjusted_nopmi_parietal.RData")

#Data frame names (for reference):
#hydroxymethylcyto.all.cing.late.nopmi
#hydroxymethylcyto.all.cing.mid.nopmi
#hydroxymethylcyto.all.pari.late.nopmi
#hydroxymethylcyto.all.pari.mid.nopmi

##Create data.frames of the 5-hmC beta data:
hydroxy_mC_C = as.data.frame(hydroxymethylcyto.c) #766212 probes
hydroxy_mC_P = as.data.frame(hydroxymethylcyto.p) #774831 probes

##Go through each row and determine if a value is zero:
row_sub = apply(hydroxy_mC_C, 1, function(row) all(row !=0 ))
row_sub1 = apply(hydroxy_mC_P, 1, function(row) all(row !=0 ))

##Subset the data to remove rows with zeros
hydroxy_C_nonzero <- hydroxy_mC_C[row_sub,] #Removes all rows w/ zero as a value.
#512094 observations.
hydroxy_P_nonzero <- hydroxy_mC_P[row_sub1,] #Removes all rows w/ zero as a value.
#527655 observations.

#Number of probes with zero for any replicate:
766212-512094 #254118 probes with at least one zero for 5-hmC in cingulate
774831-527655 #247176 probes with at least one zero for 5-hmC in parietal

#New Goals (6/11/18):
#4) Filter 5hmC data (we can filter gamlss output by new more stringent criteria for non-zero beta values as a shortcut): top sig probes and volcano plots to see if this clears up the issues.
#5) New sig probe counts based on filtered data. This chunk of code is on GitHub in P1_03_SigProbes

#Move 5-hmC beta data probe row names to new column for matching w/ significant probes:
hydroxy_C_nonzero$Probe <- rownames(hydroxy_C_nonzero) #Note: "Probe" matches column name in sig probes outputs
hydroxy_P_nonzero$Probe <- rownames(hydroxy_P_nonzero)

#Merge data.frames for filtered 5hmC beta data and significant probes using "Probe" column:
hmc.cing.late.nonzero.merged <- merge(hydroxy_C_nonzero, hydroxymethylcyto.all.cing.late.nopmi, by=c("Probe"))
hmc.cing.mid.nonzero.merged <- merge(hydroxy_C_nonzero, hydroxymethylcyto.all.cing.mid.nopmi, by=c("Probe"))
hmc.pari.late.nonzero.merged <- merge(hydroxy_P_nonzero, hydroxymethylcyto.all.pari.late.nopmi, by=c("Probe"))
hmc.pari.mid.nonzero.merged <- merge(hydroxy_P_nonzero, hydroxymethylcyto.all.pari.mid.nopmi, by=c("Probe"))

#Save merged data.frames as Rdata file:
save(hmc.cing.late.nonzero.merged,hmc.cing.mid.nonzero.merged,hmc.pari.late.nonzero.merged,hmc.pari.mid.nonzero.merged, file = "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/sig.hmc.nonzero.stage.RData")

#Write .csv files of each separate merged data.frame
write.csv(hmc.cing.late.nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Cing_Late_Nonzero_Merged_SigProbes_061018.csv")
write.csv(hmc.cing.mid.nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Cing_Mid_Nonzero_Merged_SigProbes_061018.csv")
write.csv(hmc.pari.late.nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Pari_Late_Nonzero_Merged_SigProbes_061018.csv")
write.csv(hmc.pari.mid.nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Pari_Mid_Nonzero_Merged_SigProbes_061018.csv")

#Count probes that intersect between the significant probes and filtered beta data:
library(dplyr)
c.mid.hmc.nonzero <- dplyr::intersect(hydroxy_C_nonzero$Probe, hydroxymethylcyto.all.cing.mid.nopmi$Probe)
p.mid.hmc.nonzero <- dplyr::intersect(hydroxy_P_nonzero$Probe, hydroxymethylcyto.all.pari.mid.nopmi$Probe)

c.late.hmc.nonzero <- dplyr::intersect(hydroxy_C_nonzero$Probe, hydroxymethylcyto.all.cing.late.nopmi$Probe)
p.late.hmc.nonzero <- dplyr::intersect(hydroxy_P_nonzero$Probe, hydroxymethylcyto.all.pari.late.nopmi$Probe)

#hydroxymethylcyto.all.cing.late.nopmi
#hydroxymethylcyto.all.cing.mid.nopmi
#hydroxymethylcyto.all.pari.late.nopmi
#hydroxymethylcyto.all.pari.mid.nopmi

c.mid.hmc.counts <- as.numeric(c(length(c.mid.hmc.nonzero)))
c.late.hmc.counts <- as.numeric(c(length(c.late.hmc.nonzero)))
p.mid.hmc.counts <- as.numeric(c(length(p.mid.hmc.nonzero)))
p.late.hmc.counts <- as.numeric(c(length(p.late.hmc.nonzero)))
probe_counts <- data.frame(c.mid.hmc.counts,c.late.hmc.counts,p.mid.hmc.counts,p.late.hmc.counts)
names(probe_counts) <- c("Cingulate-Mid","Cingulate-Late","Parietal-Mid","Parietal-Late")
row.names(probe_counts) <- c("5hmC non-zero probes")                         
kable(probe_counts) #Creates table for Rmarkdown output -- pretty nifty!


############################# 6/10/18 ####################################
#GOAL: Filter oxBS mle estimates of 5hmc to only probes that had 3 non-zero values in 
#at least one disease stage, instead of 3 overall -- how many would we be left with modeling?

#To examine number of non-zero probes in each disease stage, split 5-hmC data.frames 
#into 5-hmC beta data for each disease stage.
hydroxy_mC_C_1 = hydroxy_mC_C[1:3] #766212 probes
hydroxy_mC_C_2 = hydroxy_mC_C[4:6] 
hydroxy_mC_C_3 = hydroxy_mC_C[7:9] 


hydroxy_mC_P_1 = hydroxy_mC_P[1:3] #774831 probes
hydroxy_mC_P_2 = hydroxy_mC_P[4:6] 
hydroxy_mC_P_3 = hydroxy_mC_P[7:9] 

##Go through each row and determine if a value is zero:
row_sub_C_1 = apply(hydroxy_mC_C_1, 1, function(row) all(row !=0 ))
row_sub_P_1 = apply(hydroxy_mC_P_1, 1, function(row) all(row !=0 ))

row_sub_C_2 = apply(hydroxy_mC_C_2, 1, function(row) all(row !=0 ))
row_sub_P_2 = apply(hydroxy_mC_P_2, 1, function(row) all(row !=0 ))

row_sub_C_3 = apply(hydroxy_mC_C_3, 1, function(row) all(row !=0 ))
row_sub_P_3 = apply(hydroxy_mC_P_3, 1, function(row) all(row !=0 ))

##Subset the data to remove rows with zeros
hydroxy_C_1_nonzero <- hydroxy_mC_C_1[row_sub_C_1,] #Removes all rows w/ zero as a value.
#613543 probes
hydroxy_C_2_nonzero <- hydroxy_mC_C_2[row_sub_C_2,] #Removes all rows w/ zero as a value.
#616237 probes
hydroxy_C_3_nonzero <- hydroxy_mC_C_3[row_sub_C_3,] #Removes all rows w/ zero as a value.
#568130 probes

hydroxy_P_1_nonzero <- hydroxy_mC_P_1[row_sub_P_1,] #Removes all rows w/ zero as a value.
#622871 probes
hydroxy_P_2_nonzero <- hydroxy_mC_P_2[row_sub_P_2,] #Removes all rows w/ zero as a value.
#595627 probes
hydroxy_P_3_nonzero <- hydroxy_mC_P_3[row_sub_P_3,] #Removes all rows w/ zero as a value.
#627211 probes

#Move 5-hmC beta data probe row names to new column for matching w/ significant probes:
hydroxy_P_1_nonzero$Probe <- rownames(hydroxy_P_1_nonzero)
hydroxy_P_2_nonzero$Probe <- rownames(hydroxy_P_2_nonzero)
hydroxy_P_3_nonzero$Probe <- rownames(hydroxy_P_3_nonzero)

hydroxy_C_1_nonzero$Probe <- rownames(hydroxy_C_1_nonzero) #Note: "Probe" matches column name in sig probes outputs
hydroxy_C_2_nonzero$Probe <- rownames(hydroxy_C_2_nonzero) 
hydroxy_C_3_nonzero$Probe <- rownames(hydroxy_C_3_nonzero) 

#Merge filtered non-zero probes across disease stages for cing and pari. 
#NOTE: all=TRUE keeps all unique rows after merging, eventually giving a complete
#set of probes with at least one disease stage with 3 non-zero values.
hydroxy_P_nonzero_merged <- merge(hydroxy_P_1_nonzero, hydroxy_P_2_nonzero, all=TRUE)
hydroxy_P_all_nonzero_merged <- merge(hydroxy_P_nonzero_merged, hydroxy_P_3_nonzero, all=TRUE)
#696124 probes

hydroxy_C_nonzero_merged <- merge(hydroxy_C_1_nonzero, hydroxy_C_2_nonzero, all=TRUE)
hydroxy_C_all_nonzero_merged <- merge(hydroxy_C_nonzero_merged, hydroxy_C_3_nonzero, all=TRUE)
#680709 probes

#Merge data.frames for filtered 5hmC beta data and significant probes using "Probe" column:
hmc.cing.late.3nonzero.merged <- merge(hydroxy_C_all_nonzero_merged, hydroxymethylcyto.all.cing.late.nopmi, by=c("Probe"))
hmc.cing.mid.3nonzero.merged <- merge(hydroxy_C_all_nonzero_merged, hydroxymethylcyto.all.cing.mid.nopmi, by=c("Probe"))
hmc.pari.late.3nonzero.merged <- merge(hydroxy_P_all_nonzero_merged, hydroxymethylcyto.all.pari.late.nopmi, by=c("Probe"))
hmc.pari.mid.3nonzero.merged <- merge(hydroxy_P_all_nonzero_merged, hydroxymethylcyto.all.pari.mid.nopmi, by=c("Probe"))

#Save merged data.frames as Rdata file:
save(hmc.cing.late.3nonzero.merged,hmc.cing.mid.3nonzero.merged,hmc.pari.late.3nonzero.merged,hmc.pari.mid.3nonzero.merged, file = "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/sig.hmc.3nonzero.in.1stage.RData")

#Write .csv files of each separate merged data.frame
write.csv(hmc.cing.late.3nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Cing_Late_3Nonzero_in_OneStage_SigProbes_061018.csv")
write.csv(hmc.cing.mid.3nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Cing_Mid_3Nonzero_in_OneStage_SigProbes_061018.csv")
write.csv(hmc.pari.late.3nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Pari_Late_3Nonzero_in_OneStage_SigProbes_061018.csv")
write.csv(hmc.pari.mid.3nonzero.merged, "/Users/jjkoch/Dropbox/Bernstein Lab/Sequencing Data Analysis/02_Differential_Methylation/P1_Pari_Mid_3Nonzero_in_OneStage_SigProbes_061018.csv")

#Count probes that intersect between the significant probes and filtered beta data:
library(dplyr)
c.mid.hmc.3nonzero <- dplyr::intersect(hydroxy_C_all_nonzero_merged$Probe, hydroxymethylcyto.all.cing.mid.nopmi$Probe)
p.mid.hmc.3nonzero <- dplyr::intersect(hydroxy_P_all_nonzero_merged$Probe, hydroxymethylcyto.all.pari.mid.nopmi$Probe)

c.late.hmc.3nonzero <- dplyr::intersect(hydroxy_C_all_nonzero_merged$Probe, hydroxymethylcyto.all.cing.late.nopmi$Probe)
p.late.hmc.3nonzero <- dplyr::intersect(hydroxy_P_all_nonzero_merged$Probe, hydroxymethylcyto.all.pari.late.nopmi$Probe)


c.mid.hmc.counts <- as.numeric(c(length(c.mid.hmc.3nonzero)))
c.late.hmc.counts <- as.numeric(c(length(c.late.hmc.3nonzero)))
p.mid.hmc.counts <- as.numeric(c(length(p.mid.hmc.3nonzero)))
p.late.hmc.counts <- as.numeric(c(length(p.late.hmc.3nonzero)))
probe_counts <- data.frame(c.mid.hmc.counts,c.late.hmc.counts,p.mid.hmc.counts,p.late.hmc.counts)
names(probe_counts) <- c("Cingulate-Mid","Cingulate-Late","Parietal-Mid","Parietal-Late")
row.names(probe_counts) <- c("5hmC non-zero (at least 3 in one stage) probes")                         
kable(probe_counts) #Creates table for Rmarkdown output -- pretty nifty!
