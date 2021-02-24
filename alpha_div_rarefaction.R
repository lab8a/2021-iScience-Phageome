# UPDATE 2020-07-14: Set a min sample size of 1000 items and added an optional manual overdrive for rarefaction depth in the parameters
# UPDATE 2019-09-18: Filter empty columns if present
# Started on 2020-03-23
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was written to calculate alpha metrics calculated from 10K rarefactions from an otu table for the specified project and generate output tables with such calculations over each metric
# Metrics simpson, inverse shannon (diversity), ace and pielou are commented but may be activated by uncommenting where # # # is found
# Input tables are contingency tables where first column has the ID, followed by each sample (column-wise). No taxonomy is included in the last column.
# It requires the vegan library 
# Run like this:
# cat 07_Raw_contingency_tables/01_tables_no_sample_singletons/01_V3_SE-OTUs_gg97-table_lvlOTU.tsv|~/bin/R-3.6.0/bin/Rscript alpha_div_rarefaction.R 01_V3_SE-OTUs_gg97-table 15
# Tested with:
# df <- read.table("01_data/4.6k_phage_contigs_table.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1,check.names=F)
# prefix = "test"
# tot_rar = 10
n=1 # UPDATE 2020-07-14: added a scaling factor to avoid losing data when many small values are present (1000 was way too slow, consider using 100) This increases computing time by 10x. Adjust it to the smallest values.

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2) { # at least, two arguments should be included: <total_rarefactions> <output_prefix>
  stop("A minimum of 2 arguments is mandatory: cat table.tsv|Rscript alpha_div_rarefaction.R <total_rarefactions> <output_prefix> [depth (optional)]", call.=FALSE)
}
tot_rar <- as.numeric(args[1]) # Get the total number of rarefaction to carry out
prefix <- as.character(args[2]) # Get a string handle to create output names
depth <- as.numeric(args[3]) # depth is optional
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names = FALSE)
library("vegan") # we will create a rarefied set with vegan using the smallest sample size as depth with a certain number of repetitions, to get the average so we can approximate the actual composition more accurately.
df <- df[,colSums(df)>0] # NEW filter 2019-09-18: remove those with empty columns (when filtered from large tables)
df <- round(n*(df)) # UPDATE 2020-11-09: We need integers
depths_vect <- sort(colSums(df))
min_limit <- depths_vect[depths_vect>=1000][1] # UPDATE 2020-07-14 Default: set the limit to the smallest sample >999 (abort if none is larger)
if (is.na(min_limit)) {stop("Cowardly aborting: No samples of 1000 items or more. Please use a larger sample", call.=FALSE)}
df <- df[,colSums(df)>=1000] # UPDATE 2020-07-14: By default, don't use samples <1000 for rarefactions, remove them
# min_limit <- min(colSums(df)) # OPTIONAL: the depth parameter for the rarefaction is set as the depth of the smallest sample #DEPRECATED, min_limit already set at start
min_limit <- ifelse(!is.na(depth), depth*n, min_limit) # optionally, set the min_limit manually (always overrides default)
df <- df[,colSums(df)>=min_limit] # UPDATE 2020-07-14: Remove items smaller than the min_limit

samples <- colnames(df)
dfr <- df # Create a copy of current dataframe to get the right dimensions
dfr[dfr!=0] <- 0 # Then empty the newly created matrix
m_shannon <- data.frame(matrix(NA,nrow=tot_rar,ncol=ncol(df))) #Initialize empty matrices for each index
colnames(m_shannon) <- samples
# m_diversity <- m_shannon # This gets the non-log form of shannon
# m_pielou <- m_shannon # This is for the distribution among totals
# m_simpson <- m_shannon # This is for the dominance
m_observed <- m_shannon # This is the raw number of unique features
m_chao1 <- m_shannon # This is for the estimated items based on singletons (adjusted if not present)
# m_ace <- m_shannon # This is for the estimated items but based on singletons, doubletons, up to 10s, similar to chao
SUM <- apply(df>0,2,sum) # Store the sum of non null observations

for (i in 1:tot_rar){ #Create n tables
	print(i) # print current iteration for tracing purposes when running large number of rarefactions
	rar <- rrarefy(t(df),min_limit)
	dfr <- dfr+t(rar) # sum the current rarefaction to the cummulative sum table
	m_shannon[i,] <- diversity(rar,index="shannon",base=exp(1)) # Store Shannon's entropy index (uncertainty of predicting next draw) # UPDATE: 2020-11-10 This was set by default to 2.71 as base but qiime uses base 2
# 	m_diversity[i,] <-exp(m_shannon[i,]) # Store the effective diversity (non-logartithmic version of shannon, true diversity)
# 	m_pielou[i,] <- m_shannon[i,]/log(SUM)# Store Pielou's evenness (distribution of observations among items (less means more abundance variation; sensitive to sample size))
# 	m_simpson[i,] <- diversity(rar,index="simpson") # Store the Simpson's measure of dominance (P or drawing two items from same features)
	Richness <- estimateR(rar) # Get the observed features, the chao1 richness, the ACE richness, estimators along with their se
	m_observed[i,] <- Richness[1,] # Store unique feature counts
	m_chao1[i,] <- Richness[2,] # Store the number of estimated features
# 	m_ace[i,] <- Richness[4,] # Store the number of estimated features based on the those having 1 through 10
}
dfr <- dfr/tot_rar # Get a table with the averaged observations for all iterations
dfr <- dfr/n # UPDATE 2020-07-14: this was a correction added for very small values, to keep them, uncomment if required
# dfr <- floor(dfr) # UPDATE 2020-06-14: OPTIONAL: Floor results so that we avoid counting those seen <1 time on average (useful for linear scaling)
# temp <- cbind(dfr,dfresp[,9]) #add the taxonomy (need a previous version of the table (dfresp)
dfr <- as.matrix(dfr[order(rowSums(dfr),decreasing=T),]) #sort by most abundant (these may have changed due to rarefactions
# output the standardized observations
write.table(dfr,paste(prefix,"table_mean",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, row.names=T, col.names =  NA) # N export the result of the mean of n iterations as a tsv file
# write.table(cbind(dfr,dfresp[rownames(dfr),][,length(dfresp)]),paste(prefix,"mean",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, row.names=T, col.names =  NA) # alternative version with taxonomy
# Output all rarefied alphas
write.table(round(m_shannon,4),paste(prefix,"shannon",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE) # and export it as a tsv file
# write.table(round(m_diversity,4),paste(prefix,"diversity",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(round(m_pielou,4),paste(prefix,"pielou",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(round(m_simpson,4),paste(prefix,"simpson",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(m_observed,paste(prefix,"observed",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(round(m_chao1,4),paste(prefix,"chao1",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(round(m_ace,4),paste(prefix,"ace",tot_rar,"resamples.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
