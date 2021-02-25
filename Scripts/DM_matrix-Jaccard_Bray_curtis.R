# Started: 2020-08-18
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script stub was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# The script is intended to create a Jaccard (composition: absence/presence) and a Bray-Curtis (abundance) dissimilarity matrices from a contingency table. No statistics calculation is carried out
# The input is preferably a previously-rarefied or normalized table (by default, no taxonomy is present in the table and the header is in row 1. Feature names are in column 1.

# Run as follows:
# cat table.tsv|Rscript DM_matrix-Jaccard_Bray_curtis.R <output_prefix>

# Tested with command:
# cat /home/rod/Desktop/16S/10_sparsity_reduction/01_Per_sample_0.001_mt_2_sam/01_transformed_tables/Filt_table_lvl5-rar-10000-depth-7000-perm.tsv|Rscript DM_matrix-Jaccard_Bray_curtis.R test
# In R:
df <- read.table("/home/rod/Desktop/16S/10_sparsity_reduction/01_Per_sample_0.001_mt_2_sam/01_transformed_tables/Filt_table_lvl5-rar-10000-depth-7000-perm.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1,check.names=F)
prefix="test"

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) { # at least, one arguments should be included: <out_name_prefix>
  stop("A minimum of 1 argument is mandatory: cat table.tsv|Rscript multi_ordination_methods.R <output_file>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1,check.names=F) #Read the table, row names expected in first column, matching colnames in first row

library("vegan") # Distance matrix calculations are included the vegan package
# Create the matrices
bray <- vegdist(t(df),method="bray",upper=TRUE,diag=TRUE)
jaccard <- vegdist(t(df),method="jaccard",upper=TRUE,diag=TRUE, binary=TRUE)
write.table(as.matrix(bray),paste(prefix,"bray_dm.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(as.matrix(jaccard),paste(prefix,"jaccard_dm.tsv", sep="-"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
