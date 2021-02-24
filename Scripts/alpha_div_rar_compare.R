# UPDATE 2020-11-16: This is a copy of script alpha_div_rar_compare.R created to specifically fix some parameters for figure 3
# UPDATE 2020-11-10: Added a small fix for shannon: changing ln (R's default) to log2 (qiime 1.9's default)
# UPDATE 2020-08-12: The script now creates a graph with the sample medians (one per sample) at the end of each pdf file. Also, I commented large tables with five " #" to save space. These correspond to per-row statistics
# Started on 2020-07-18
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This is a new version of the script to collate, graph and calculate group statistics between alpha diversity values
# The script is intended to create alpha graphs per samples and for user-defined grpups (provided as arguments)
# To deal with NAs produced by ACE and other indices, only rows with no NAs can be considered. To cope with this, NAs are left out and the user is allowed to test a cutoff of minimum non-NA values that should exist per sample for them to be included in the analyses
# There are two parts of the script:
# 1.- The whole table is considered (except empty columns). To cope with NAs, the total number of non-NA items in the smallest sample is used as the cutoff limit. This is used for group statistics and graphs. Sample boxplots and comparisons are still calculated for all obs as well.
# 2.- The table is prefiltered so that only samples having a minimum number of valid items (user-defined) are considered. The new matrix will include dim(min_nonNAs rows x passing_samples) and statistics are calculated from the trimmed table.
# The input is a matrix where each row contains the alpha metrics (any) for a table (e.g. rarefactions or jacknife recalculations) and each column is a sample
# Sample headers are expected to be in the first line
# The user mst define the following parameters (with examples)
# prefix = "test" # output prefix, this can contain a full or partial path
# sub = "ACE" # The name of the metric in use, for the graph titles
# min_nonNAs = 7000 #Minimum number of non-NA items per sample that should be used for additional tables and graphs
# # groups <- c("Z1","Z2","Z3","Z4","Z5")
# # groups <- c("2015","2016","2017","2018")
# # groups <- c("H","I","S")
# # groups <- c("HE","IE","HL","IL","HM","IM")
# # groups <- c("E","L","M")
# groups <- c("2015_Z1_H","2015_Z1_I","2015_Z2_H","2015_Z2_I","2016_Z2_H","2016_Z2_I","2015_Z3_I","2015_Z3_S","2016_Z3_S","2017_Z3_H","2017_Z3_I","2017_Z3_S","2018_Z4_I","2018_Z5_H")
# groups <- c("2015_Z1_I","2015_Z2_H","2016_Z2_H","2016_Z2_I","2017_Z3_H","2017_Z3_I","2018_Z4_I","2018_Z5_H")
# # Tested with:
# df <- read.table("02_alpha_div/10K_min_sample-round_int/4.6k_phage_contigs_alpha_shannon_10000_resamples.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names=F)
# Run like this:
# cat 13_alpha_diversity/depth_3250/Bact_lvl3_ace_10000_resamples.tsv|Rscript alpha_div_rar_compare.R 50 01_V3_SE-OTUs_gg97-table observed

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, two arguments should be included: <min_nonNAs> <prefix_output>  <name_of_alpha_metric>
  stop("A minimum of 4 arguments are mandatory: cat table.tsv|Rscript alpha_div_rar_compare.R <min_nonNAs> <prefix_output>  <name_of_alpha_metric> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>", call.=FALSE)
}
min_nonNAs <- as.numeric(args[1]) # Get the cutoff (float) for nin non-NA items (passing samples have at least this number of valid items). This should be positive. If larger than max counts use the max available.
prefix <- as.character(args[2]) # Get a string handle to create output names
sub <- as.character(args[3])  # and the name of the metric (free strings)
groups <- as.character(args[4:length(args)]) # Create a vector with the name of the groups (this should be included in the sample names and be exclusive to each group for this to work)
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names=F) # Load input table
if(sub=="shannon"){df<- log2(exp(1)^df)} # UPDATE 2020-11-10: Added to change scales for shannon into ln (default in R) to log2 (default in qiime)
# colSums(apply(df,2, is.na))/nrow(df)*100
#	C01	C02	C03	C04	C05	L01	L02	L03	L04	L05	L06 
#   0.00   0.00   0.00   0.02   0.37   0.03   0.00   0.01   0.19   2.04   1.63 
#	L07	L08	L09	L10	S02	S03	S04	S05	S06	S07	S08 
#  44.12   2.50   2.01   0.00   0.15  57.73   0.83  84.44  19.64  20.35   0.31 
#	S09	S10	S11	S12	S14	S16	V01	V02	V03	V04	V05 
#   2.60   2.14   0.03   0.63   1.21   0.15  34.61   0.00   0.07  14.66   0.00 
#	V06	V07	V08	V09	V10
#   2.42   4.28   0.10   0.17 100.00
print("Loaded parameters:")
print(paste("Prefix:",prefix))
print(paste("Metric name:",sub))
print(paste("Minimum non-NA values expected:",min_nonNAs))
print("Groups that should be present:")
print(groups)
####################################### Pre-load #######################################
 ### Define colours ### 
# Colors will only be used if the total number of groups is the same or less than the number of colors
# DEPRECATED colors
# group_cols <- c("coral1","red","cornflowerblue","blue","coral1","cornflowerblue","springgreen4","darkorchid3")
# group_cols2 <- c("bisque","pink","darkslategray2","dodgerblue2","bisque","darkslategray2","chartreuse2","orchid")
# pair or color sets used for different comparisons
# group_bord <- c("coral1","red","cornflowerblue","blue") # I've used this for random 5 group items
# group_cols <- c("bisque","pink","darkslategray2","dodgerblue2")
# group_bord <- c("black","red","green","blue","magenta") # I've used this for runs
# group_cols <- c("gray","pink","darkolivegreen1","dodgerblue2","orchid")
# group_bord <- c("coral1","cornflowerblue","turquoise1") # I've used these for ponds
# group_cols <- c("bisque","lightblue","darkslategray1")
group_bord <- c("cornflowerblue","brown2","limegreen") # I've used these for H, O, OMS
group_cols <- c("lightblue","lightpink","darkseagreen2")
# group_bord <- c("springgreen4","darkorchid3","tan1") # I've used these for organs
# group_cols <- c("chartreuse2","orchid","orangered")
# # group_bord <- c("navyblue","hotpink","forestgreen","orangered") # I've used these for years
# group_cols <- c("dodgerblue4","plum","greenyellow","tan1")
# group_bord <- c("springgreen4","darkorchid3","darkgreen","magenta","seagreen3","violetred") # I've used these for 2organs/3ponds
# group_cols <- c("chartreuse2","orchid","limegreen","orchid","darkseagreen1","pink")
# group_bord <- c("red4","navyblue","orangered4","darkorchid4","red1","cornflowerblue","darkslateblue","black","gray30","tan4","royalblue4","darkslategray","darkorchid","gold2") # I've used this for years/organs for grps: "2015_Z1_H 2015_Z1_I 
# group_cols <- rep("white",length(group_bord))
 ### Define Functions ###
# Predefine some important global objects
prepare_data <- function(df) { # create presets (sample totals and groups # MOD: ignored empty groups
	samples <- names(df) # extract names
	grp_sam <- sapply(groups,function(x) grep(x,samples)) # create index of samples in each group
	grp_less <- sapply(grp_sam,length) # and count the totals # This will now consider empty groups
	grp_all <- grp_less # Since this was modified, we still required the old version
	grp_sam <- grp_sam[(grp_less>=1)] # Remove groups that ended up with 0 samples
	grp_less <- sapply(grp_sam,length) # This will now consider only non-empty groups
	len <- length(unlist(grp_sam)) # Get total items
	grp <- rep("Other",len) # create a template for considering non-group items
	for(i in 1:length(grp_sam)){ # This dual cycle gets the actual sample distribution for all groups
		for(j in 1:length(grp_sam[[i]])){
		grp[grp_sam[[i]][j]] <- names(grp_less)[i];
		}
	}
	out <-list(grp_sam,grp_less,as.factor(grp),grp_all)
	names(out) <- c("grp_sam","grp_less","grp", "grp_all")
	return(out)
}
set_cols <- function(grp_sam, grp_less){
	# Assign colors (if no larger than the input color vectors)
	group_cols <- group_cols[1:length(groups)][grp_less>0] # Trim to grp total, then mask those not present
	group_bord <- group_bord[1:length(groups)][grp_less>0] # idem
	len <- length(unlist(grp_sam))
	bord <- cols <- unlist(grp_sam)
	cols <- rep("gray",len) # set default colors
	bord <- rep("black",len)
	for(i in 1:length(grp_sam)){ # This dual cycle gets the actual color for that sample
# 		print(i)
		for(j in 1:length(grp_sam[[i]])){
# 		print(j)
		bord[grp_sam[[i]][j]] <- group_bord[i];
		cols[grp_sam[[i]][j]] <- group_cols[i];
		}
	}
	rep <- sapply(grp_sam,length) # ALTERNATIVE VERSION: this uses a block color (for ordered samples)
	cols_ord <- rep(group_cols,times=rep)
	bord_ord <- rep(group_bord,times=rep)
	out <-list(group_cols,group_bord,cols,bord,cols_ord,bord_ord)
	names(out) <- c("group_cols","group_bord","cols","bord","cols_ord","bord_ord")
	return(out)
}
	
noNAs <-function(vector,n) { # Function to keep only non-NA samples and their sizes, requires vector and cutoff for total items (n). Can be used to sort items and leave NAs in the last rows
	noNAs <- vector[!is.na(vector)] # remove NAs
	noNAs <- noNAs[1:n] # Keep only the first n items (these are already randomized)
	new_total <- length(noNAs)
	return(noNAs)
}
pairwise_mannwhitney <- function(matrix, grp_sam, g1,g2){ # This gets the pvalue of a Mann Whitney test of two groups of values
	i <- apply(matrix, 1, function(x) wilcox.test(as.numeric(x[grp_sam[[g1]]]),as.numeric(x[grp_sam[[g2]]]),paired=F)$p.value)
# 	out <- median(i)
# 	out <- c(out, median(p.adjust(i, method="fdr")))
	out <- cbind(i, p.adjust(i, method="fdr"))
	colnames(out) <- c(paste(names(grp_sam)[g1],"vs",names(grp_sam)[g2],"pval",sep="_"),paste(names(grp_sam)[g1],"vs",names(grp_sam)[g2],"qval",sep="_"))
	return(out)
}
 # Calculate sample medians and produce single group statistics
stat_s <- function(med_trim,grp_sam){ # Requires a vector with one value per sample and the sample to group map (list object)
	temp <- c(names(med_trim),names(grp_sam))
	for(i in 1:length(grp_sam)){med_trim <- c(med_trim,median(med_trim[grp_sam[[i]]]))} # Calculate group medians
	if(length(grp_sam)>1){ # This was added as an exception if there is only one group, since we cannot calculate statistics for 1 grp
		all_comp <- combn(length(grp_sam),2) # set permutations
		for(i in 1:ncol(all_comp)){ # cycle throught the actual pairs (all permutations)
			med_trim <- c(med_trim,wilcox.test(as.numeric(med_trim[grp_sam[[all_comp[1,i]]]]),as.numeric(med_trim[grp_sam[[all_comp[2,i]]]]),paired=F)$p.value) # and statistics
			temp <- c(temp, paste(names(grp_sam)[all_comp[1,i]],"vs",names(grp_sam)[all_comp[2,i]],"pval",sep="_"))
		}
	}
	names(med_trim) <- temp
	return(med_trim)
}
# Calculate group medians and pairwise Mann-Whitney tests
stat_g <- function(matrix,grp_sam){
	# Calculate group medians for trimmed table
	grp_medians <- NULL # Create an empty object container for the medians
	for(i in 1:length(grp_sam)){ # cycle through the sample map per group
		grp_medians <- cbind(grp_medians,apply(matrix, 1, function(x){median(x[grp_sam[[i]]])})) # use the indices in grp_sam to calculate groupmedians
	}
	if(length(grp_sam)>1){ # This was added as an exception if there is only one group, since we cannot calculate statistics for 1 grp
		all_comp <- combn(length(grp_sam),2) # Create a small matrix of pairwise permutations
		MW_stat <- NULL # This object will hold the statistical comparisons (Mann-Whitney)
		for(i in 1:ncol(all_comp)){ # cycle throught the actual pairs (all permutations)
			x <- pairwise_mannwhitney(matrix,grp_sam,all_comp[1,i],all_comp[2,i])
			MW_stat <- cbind(MW_stat, x)
		}
		colnames(grp_medians) <- names(grp_sam)
	} else {
		MW_stat <- rep(NA,nrow(matrix))
	}
	out <- list(matrix,grp_medians,round(MW_stat,5))
	names(out) <- c("mat","grp_medians","stats")
	return(out)
}

####################################### MAIN #######################################
############ Part 1: Use all items in the table ############
 ### Prefilter empty samples ###
min_nonNAs <- ifelse(min_nonNAs>nrow(df),nrow(df),min_nonNAs) # If the value for min_nonNAs inputted by the user is larger than the total rows present, cap it to match it
df <- df[,(colSums(!apply(df,2, is.na)))>0] # By default, remove colums with only NAs
if(ncol(df)<1){stop(paste0("Aborting: No samples have at least ",min_nonNAs,"% NA values"), call.=FALSE)} #  abort if all are missing
df <- as.data.frame(apply(df,2,function(x) noNAs(x,nrow(df)))) # Reorder the matrix so that NAs are always in the last rows
presets <- prepare_data(df)

if(length(groups)<=length(group_bord)){ # Reassign colors to existing groups (to make them the same between graphs). This is only carried out if the number of defined colors (solid and border) is the same than the number of groups.
	cols <- set_cols(presets$grp_sam,presets$grp_all)
	s_col <- cols$cols
	s_bor <- cols$bord
	s_col2 <- cols$cols_ord
	s_bor2 <- cols$bord_ord
	g_col <- cols$group_cols
	g_bor <- cols$group_bord
} else {
	s_col <- s_col2 <- g_col <- "gray"
	s_bor <- s_bor2 <- g_bor <- "black"
}
# # # # # pdf(paste0(prefix,"_alpha-",sub,"-1_All_obs.pdf"))
# # # # # boxplot(df,col=s_col,border=s_bor,outline=F,main=paste0(basename(prefix), " - ",sub, " - Full Samples"),xaxt='n')
# # # # # mtext("Total obs per sample in x axis")
# # # # # axis(1,las=2,at=seq(1,ncol(df),1),labels=paste(colSums(!apply(df,2, is.na)),names(df)),cex.axis=0.8)
# # # # # boxplot(df[,unlist(presets$grp_sam)],col=s_col2,border=s_bor2,outline=F,main=paste0(basename(prefix), " - ",sub, " - Full Samples"),xaxt='n')
# # # # # mtext("Total obs per sample in x axis")
# # # # # axis(1,las=2,at=seq(1,ncol(df),1),labels=paste(colSums(!apply(df,2, is.na)),names(df[,unlist(presets$grp_sam)])),cex.axis=0.8)
# # # # # dev.off()
# med_user <- apply(plotdata,2,median)
# med_min <- 
# # # # # # Raw_samples <- t(round(stat_s(apply(df,2,function(x) median(x[!is.na(x)])), presets$grp_sam),5)) # UNCOMMENT THIS TO CREATE RAW STATS TABLES
# # # # # write.table(Raw_samples,paste0(prefix,"_alpha-",sub,"-1_All_obs-stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE) # and export it as a tsv file 

############ Part 2: Use min sample size for groups ############
m_obs <- min(colSums(!apply(df,2, is.na)));plotdata <- df[1:m_obs,] # IMPORTANT: Define a new total number of obs and subset the data
# First plot samples
pdf(paste0(prefix,"_alpha-",sub,"-2_",m_obs,"_obs.pdf"))
par(las=1)
# boxplot(plotdata,col=s_col,border=s_bor,outline=F,main=paste0(basename(prefix), " - ",sub, " - Samples at smallest # obs: ",m_obs),xaxt='n')
# axis(1,las=2,at=seq(1,ncol(df),1),labels=,names(df),cex.axis=0.8)
# boxplot(plotdata[,unlist(presets$grp_sam)],col=s_col2,border=s_bor2,outline=F,main=paste0(basename(prefix), " - ",sub, " - Full Samples"),xaxt='n')
# axis(1,las=2,at=c(1:ncol(plotdata[,unlist(presets$grp_sam)])),labels=names(plotdata[,unlist(presets$grp_sam)]),cex.axis=0.8)
sub_mat <- stat_g(plotdata,presets$grp_sam) # Get group medians and statistics
# Next, plot groups
# boxplot(sub_mat$grp_medians,col=g_col,border=g_bor,outline=F,main=paste0(basename(prefix), " - ",sub, " - Group medians at smallest # obs: ",m_obs),xaxt='n') # This gets 10k or whatever medians per group
# axis(1,las=2,at=c(1:length(colnames(sub_mat$grp_medians))),labels=colnames(sub_mat$grp_medians),cex.axis=1)
sd=1.2
p <- seq(100,1100,100)
values <- sapply(presets$grp_sam,function(x) apply(plotdata,2,median)[x])
boxplot(values,col=g_col,border=g_bor,outline=F,main="Phage Contigs - Observed Richness",boxwex=0.7,frame=F,xaxt='n',yaxt='n',lwd=2,cex.axis=1.5,cex.lab=1.5,ylim=range(p),ylab="Unique observed contigs") # Use the list to extract the corresponding medians (one median per sample)
library(beeswarm)
# stripchart(sapply(presets$grp_sam,function(x) apply(plotdata,2,median)[x]), vertical = TRUE, method = "overplot", add = TRUE, pch = 1, col = g_bor)
# stripchart(sapply(presets$grp_sam,function(x) apply(plotdata,2,median)[x]), vertical = TRUE, method = "overplot", add = TRUE, pch = 20, col = g_bor,cex=0.5)
axis(1,las=1,at=c(1:length(colnames(sub_mat$grp_medians))),labels=gsub("-","",colnames(sub_mat$grp_medians)),cex.axis=1)
axis(2, las=1, at=p,cex.axis=1.3)
beeswarm(values,add=T,col="white",method="swarm",cex=2,pch=16,spacing=sd)
beeswarm(values,add=T,col=g_bor,method="swarm",cex=2,pch=1,spacing=sd)
beeswarm(values,add=T,col=g_bor,method="swarm",cex=2,pch=20,spacing=sd)
# axis(1,las=1,at=c(1:3),labels=c("NW","O","OMS"),cex.axis=1.5)
dev.off()
# And output the statistics generated so far
# # # # # # UNCOMMENT THIS TO CREATE WHOLE STATS TABLES
# # # # # write.table(as.data.frame(sub_mat),paste0(prefix,"_alpha-",sub,"-2_",m_obs,"_obs-stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE) # and export it as a tsv file 
# Statistics are also calculated with a single value per sample, with the median. This calculus is carried out once, for there is only one median per sample.
med_min <- t(round(stat_s(apply(plotdata,2,median), presets$grp_sam),5))
write.table(t(med_min[,c(29:34)]),paste0(prefix,"_alpha-",sub,"-2_",m_obs,"_obs-stats_single.tsv"), sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE) # and export it as a tsv file 

############ Part 3: Use min user non-NAs items for groups ############
# # # # # df <- df[,(colSums(!apply(df,2, is.na)))>=min_nonNAs] # Remove columns with less than this total non-NAs items
# # # # # if(ncol(df)<1){stop(paste0("Aborting: No samples have at least ",min_nonNAs,"% NA values"), call.=FALSE)} # abort if none pass
# # # # # df <- df[1:min_nonNAs,] # Now keep only that number of obs (the table had been sorted before so there shouldn't be any NAs afterwards
# # # # # presets <- prepare_data(df)
# # # # # 
# # # # # if(length(groups)<=length(group_cols)){ # Reassign colors to existing groups (to make them the same between graphs). This is only carried out if the number of defined colors (solid and border) is the same than the number of groups.
# # # # # 	cols <- set_cols(presets$grp_sam,presets$grp_all)
# # # # # 	s_col <- cols$cols
# # # # # 	s_bor <- cols$bord
# # # # # 	g_col <- cols$group_cols
# # # # # 	g_bor <- cols$group_bord
# # # # # } else {
# # # # # 	s_col <- g_col <- "gray"
# # # # # 	s_bor <- g_bor <- "black"
# # # # # }
# # # # # pdf(paste0(prefix,"_alpha-",sub,"-3_",min_nonNAs,"_obs.pdf"))
# # # # # boxplot(df,col=s_col,border=s_bor,outline=F,main=paste0(basename(prefix), " - ",sub, " - Samples at smallest # obs: ",min_nonNAs),xaxt='n')
# # # # # mtext("Total obs per sample in x axis")
# # # # # axis(1,las=2,at=seq(1,ncol(df),1),labels=names(df),cex.axis=0.8)
# # # # # boxplot(df[,unlist(presets$grp_sam)],col=s_col2,border=s_bor2,outline=F,main=paste0(basename(prefix), " - ",sub, " - Full Samples"),xaxt='n')
# # # # # axis(1,las=2,at=c(1:ncol(df[,unlist(presets$grp_sam)])),labels=names(df[,unlist(presets$grp_sam)]),cex.axis=0.8)
# # # # # # Also for group medians
# # # # # sub_mat <- stat_g(df,presets$grp_sam) # Get group medians and statistics
# # # # # # Next, plot groups
# # # # # boxplot(sub_mat$grp_medians,col=g_col,border=g_bor,outline=F,main=paste0(basename(prefix), " - ",sub, " - Group medians at smallest # obs: ",min_nonNAs),xaxt='n') # This gets 10k or whatever medians per group
# # # # # axis(1,las=2,at=c(1:length(colnames(sub_mat$grp_medians))),labels=colnames(sub_mat$grp_medians),cex.axis=1)
# # # # # # Statistics are also calculated with a single value per sample, with the median. This calculus is carried out once, for there is only one median per sample.
# # # # # boxplot(sapply(presets$grp_sam,function(x) apply(df,2,median)[x]),col=g_col,border=g_bor,outline=F,main=paste0(basename(prefix), " - ",sub, " - Sample medians at smallest # obs: ",min_nonNAs),xaxt='n') #Use the list to extract the corresponding medians (one median per sample)
# # # # # axis(1,las=2,at=c(1:length(colnames(sub_mat$grp_medians))),labels=colnames(sub_mat$grp_medians),cex.axis=1)
# # # # # dev.off()
# # # # # # And output the statistics generated for the table
# # # # # out <- as.data.frame(sub_mat)
# # # # # # # # # # # UNCOMMENT THIS TO CREATE WHOLE STATS TABLES
# # # # # # # # # # write.table(out,paste0(prefix,"_alpha-",sub,"-3_",min_nonNAs,"_obs-stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE) # and export it as a tsv file 
# # # # # med_user <- t(round(stat_s(apply(df,2,median), presets$grp_sam),5))
# # # # # write.table(med_user,paste0(prefix,"_alpha-",sub,"-3_",min_nonNAs,"_obs-stats_single.tsv"), sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE) # and export it as a tsv file 
