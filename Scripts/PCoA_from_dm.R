# UPDATE 2020-02-03: Changed fix on casting problem to add it as an exception based on class object "matrix"
# UPDATE 2020-01-30: Fixed casting problem when loading the input matrix and preparing the data (now it forces it as type: dataframe). Also removed a minimum sample cap, which shouldn't be here
# Started on 2020-07-18
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# The script is intended to plot PCoA ordinations using a distance/dissimilarity matrix and paint by groups where 0=max simmilarity
# Groups are defined by the user (colors are defined in-script)
# The input is a square simmetrical dissimilarity matrix such as Jaccard, Bray-Curtis, Manhattan or UniFrac (qiime-calculated)

# Run as follows:
# cat dm.tsv|Rscript PCoA_from_dm.R <prefix_output>  <name_of_beta_metric> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>

# Tested with command:
# cat 12_beta_diversity/depth_7000/Ordination/Filt_table_lvl5-bray_dm.tsv|Rscript PCoA_from_dm.R test Bray-Curtis HE IE HL IL HM IM
# Test in R:
# df <- read.table("12_beta_diversity/depth_7000/Ordination/Filt_table_lvl5-bray_dm.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
# prefix="test"
# metric = "Bray-Curtis" # The name of the metric in use, for the graph titles
# # groups <- c("Z1","Z2","Z3","Z4","Z5")
# # groups <- c("2015","2016","2017","2018")
# # groups <- c("H","I","S")
# groups <- c("HE","IE","HL","IL","HM","IM")
# # groups <- c("E","L","M")
# # groups <- c("2015_Z1_H","2015_Z1_I","2015_Z2_H","2015_Z2_I","2016_Z2_H","2016_Z2_I","2015_Z3_I","2015_Z3_S","2016_Z3_S","2017_Z3_H","2017_Z3_I","2017_Z3_S","2018_Z4_I","2018_Z5_H")
# # groups <-("2015_Z1_I", "2015_Z2_H", "2016_Z2_H", "2016_Z2_I", "2017_Z3_H", "2017_Z3_I", "2018_Z4_I", "2018_Z5_H")
# prefix="test"
# metric = "unifrac"
# groups2 <- c("Aga","Inu", "Comer", "Cont")
# groups <- c("int", "hep")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, two arguments should be included: <min_nonNAs> <prefix_output>  <name_of_alpha_metric>
  stop("A minimum of 4 arguments are mandatory: cat table.tsv|Rscript PCoA_from_dm.R <prefix_output>  <name_of_beta_metric> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
metric <- as.character(args[2])  # and the name of the metric (free strings)
groups <- as.character(args[3:length(args)]) # Create a vector with the name of the groups (this should be included in the sample names and be exclusive to each group for this to work; min 2 groups)
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # 
print("Loaded parameters:")
print(paste("Prefix:",prefix))
print(paste("Metric name:",metric))
print("Groups that should be present:")
print(groups)

####################################### Pre-load #######################################
 ### Define colours ### 
# Colors will only be used if the total number of groups is the same or less than the number of colors
# pair or color sets used for different comparisons
# group_cols <- c("coral1","red","cornflowerblue","blue") # I've used this for random 4 group items
 group_cols <- c("black","red","green","blue","magenta") # I've used this for runs (5 grps)
# group_cols <- c("coral1","cornflowerblue","turquoise1") # I've used these for ponds (3 grps)
# group_cols <- c("springgreen4","darkorchid3","tan1") # I've used these for organs (3 grps)
# group_cols <- c("navyblue","hotpink","forestgreen","orangered") # I've used these for years (4 grps)
# group_cols <- c("dodgerblue4","plum","greenyellow","tan1") # Alternative for years (DEPRECATED)
# group_cols <- c("chartreuse2","darkorchid3","darkgreen","magenta","seagreen3","violetred") # I've used these for 2organs/3ponds
# group_cols <- c("chartreuse2","orchid","limegreen","orchid","darkseagreen1","pink") #alternative for organ/ponds (NOW DEPRECATED)
# group_cols <- c("navyblue","orangered4","red1","cornflowerblue","tan4","royalblue4","darkorchid","gold2") # I've used this for years/runs/organs
# group_cols <- c("red4","navyblue","orangered4","darkorchid4","red1","cornflowerblue","darkslateblue","black","gray30","tan4","royalblue4","darkslategray","darkorchid","gold2") # alternative version for years/runs/organs (DEPRECATED)
# group_cols <- rep("gray",length(group_cols))
 ### Define Functions ###
# Predefine some important global objects
prepare_data <- function(df) { # create presets (sample totals and groups # MOD: ignored empty groups
	samples <- names(df) # extract names
	grp_sam <- sapply(groups,function(x) grep(x,samples)) # create index of samples in each group 
	if(class(grp_sam)=="matrix"){ # (UPDATE 2020-02-03: Added an exception to fix the casting problem when input groups have the same number of items (as this produces a matrix instead of a list, we force it again into a list with the corresponding names)
		temp <- colnames(grp_sam)
		grp_sam <- lapply(seq_len(ncol(grp_sam)), function(i) grp_sam[,i])
		names(grp_sam) <- temp
	}
	grp_less <- sapply(grp_sam,length) # and count the totals # This will now consider empty groups
	grp_all <- grp_less # Since this was modified, we still require the old version
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
	len <- length(unlist(grp_sam))
	cols <- rep("gray",len) # set default colors
	for(i in 1:length(grp_sam)){ # This dual cycle gets the actual color for that sample
		for(j in 1:length(grp_sam[[i]])){
		cols[grp_sam[[i]][j]] <- group_cols[i];
		}
	}
	rep <- sapply(grp_sam,length) # ALTERNATIVE VERSION: this uses a block color (for ordered samples)
	cols_ord <- rep(group_cols,times=rep)
	out <-list(group_cols,cols,cols_ord)
	names(out) <- c("group_cols","cols","cols_ord")
	return(out)
}
group_statistics <- function(dm, group_presets){ # This is a general launcher to create ad hoc and post hoc permutations
	# First is the ad hoc test (consider all categories in this group comparison):
	found_grp <- group_presets$grp_less # get which of the requested groups are present
	whole_mat <- dm[unlist(group_presets$grp_sam),] # and subset the matrix columnwise
	whole_mat <- whole_mat[,unlist(group_presets$grp_sam)] # as well as rowwise
	all_factors <- as.factor(rep(names(found_grp),times=found_grp)) # create sample map
	stats <- single_comparison(whole_mat,all_factors) # run adonis and anosim
	if(length(group_presets$grp_less)>2){ # Carry out post hoc test if more than 2 groups are present
		perm <- combn(names(found_grp),2) # Determine all pemutations to be carryed out
		for(it in 1:ncol(perm)){ # For each permutation
			sub_factors <- as.factor(rep(names(group_presets$grp_less[perm[,it]]),times=group_presets$grp_less[perm[,it]]))
			sub_mat <- dm[unlist(group_presets$grp_sam[perm[,it]]),] # and subset the matrix columnwise
			sub_mat <- sub_mat[,unlist(group_presets$grp_sam[perm[,it]])] # as well as rowwise
			stats <-rbind(stats,single_comparison(sub_mat,sub_factors)) # run adonis and anosim
		}
	}
	return(stats)
}
single_comparison <- function(inmat, grp_factors){
# 	1.- Adonis (Get R2 and pvalue):
	adonis_grp <- adonis(inmat ~ grp_factors,method=inmat,permutations=1000)
	adonis_inmat_R2 <- round(adonis_grp$aov.tab[1,5],4)
	adonis_inmat_pval <- round(adonis_grp$aov.tab[1,6],6)
# 	2.- Anosim (Get R and pvalue):
# 	a) Jaccard (composition, presence-absense)
	anosim_grp <- anosim(inmat, grp_factors, permutations=1000) # anosim's iterations are longer than adonis' so adjust accordingly (remember we'll be carrying out this step 
	anosim_inmat_R <- round(anosim_grp$statistic,4)
	anosim_inmat_pval <- round(anosim_grp$signif,6)
	out_table <- cbind("Comp"=toString(levels(grp_factors)), "Adonis-R2"=round(adonis_inmat_R2,3),"Anosim-R"=round(anosim_inmat_R,3),"Adonis-pval"=adonis_inmat_pval, "Anosim-pval"=anosim_inmat_pval)
	return(out_table)
}

create_PCoA <- function(dm){ # This function plots the ordination plots for the first 3 dimensions' permutations
	inmat_PCoA <- cmdscale(dm, eig = T, k=ncol(df)-1)
	inmat_PCoA.axes <- round(inmat_PCoA$eig*100/sum(inmat_PCoA$eig),1)
	pdf(paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"PCoA.pdf", sep="_"))
	par(las=1)
	# Set ranges to expand a little to include complete names (adjust factor "fix" as required)
	r1 <- range(inmat_PCoA$points[,1]) #original ranges
	r2 <- range(inmat_PCoA$points[,2])
	r3 <- range(inmat_PCoA$points[,3])
	fix=0.05
	dim1 <- c(r1[1]-diff(r1)*fix,r1[2]+diff(r1)*fix) # Fixed ranges (adjust fix as desired
	dim2 <- c(r2[1]-diff(r2)*fix,r2[2]+diff(r2)*fix)
	dim3 <- c(r3[1]-diff(r3)*fix,r3[2]+diff(r3)*fix)
	# First plot dimensions 1 and 2
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"),xlim=dim1,ylim=r2) #using original range for y axis limit, fixed for x axis
	mtext(paste("Dimensions 1 and 2 - Anosim R:", compare[1,3], "p =", compare[1,5])) # Print only anosim
	orditorp(inmat_PCoA,display="sites",cex=0.8,air=0.01,col=s_col, choices=c(1,2))
	ordihull(inmat_PCoA, groups=as.factor(s_col), draw="polygon", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,2),alpha=50)
	# Now dimensions 1 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"),xlim=dim1,ylim=r3)
	mtext(paste("Dimensions 1 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	orditorp(inmat_PCoA,display="sites",cex=0.8,air=0.01,col=s_col, choices=c(1,3))
	ordihull(inmat_PCoA, groups=as.factor(s_col), draw="polygon", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,3),alpha=50)
	# Now dimensions 2 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"),xlim=dim2,ylim=r3)
	mtext(paste("Dimensions 2 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	orditorp(inmat_PCoA,display="sites",cex=0.8,air=0.01,col=s_col, choices=c(2,3))
	ordihull(inmat_PCoA, groups=as.factor(s_col), draw="polygon", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(2,3),alpha=50,,xlim=dim2,ylim=dim3)
	# Now with bullets only
	# First plot dimensions 1 and 2
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"))
	mtext(paste("Dimensions 1 and 2 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,1], inmat_PCoA$points[,2], col=s_col, pch=1,cex=1.1)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	# Now dimensions 2 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"))
	mtext(paste("Dimensions 1 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,1], inmat_PCoA$points[,3], col=s_col, pch=1,cex=1.1)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	# Now dimensions 2 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"))
	mtext(paste("Dimensions 2 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,2], inmat_PCoA$points[,3], col=s_col, pch=1,cex=1.1)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	# Now with ellipses (sd with conf 0.95 over chi2 distr)
	# First plot dimensions 1 and 2
	fix=0.35
	dim1 <- c(r1[1]-diff(r1)*fix,r1[2]+diff(r1)*fix)
	dim2 <- c(r2[1]-diff(r2)*fix,r2[2]+diff(r2)*fix)
	dim3 <- c(r3[1]-diff(r3)*fix,r3[2]+diff(r3)*fix)
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"),xlim=dim1,ylim=r2)
	mtext(paste("Dimensions 1 and 2 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,1], inmat_PCoA$points[,2], col=s_col, pch=1,cex=1.1)
	ordiellipse(inmat_PCoA, groups=as.factor(s_col), draw="lines", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,2),kind="sd",conf=0.95)
	ordibar(inmat_PCoA, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,2),kind="sd",conf=0.05)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	# Now dimensions 1 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"),xlim=dim1,ylim=r3)
	mtext(paste("Dimensions 1 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,1], inmat_PCoA$points[,3], col=s_col, pch=1,cex=1.1)
	ordiellipse(inmat_PCoA, groups=as.factor(s_col), draw="lines", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,3),kind="sd",conf=0.95)
	ordibar(inmat_PCoA, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,3),kind="sd",conf=0.05) # get the group centroid
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	# Now dimensions 2 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"),xlim=dim2,ylim=r3)
	mtext(paste("Dimensions 2 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,2], inmat_PCoA$points[,3], col=s_col, pch=1,cex=1.1)
	ordiellipse(inmat_PCoA, groups=as.factor(s_col), draw="lines", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(2,3),kind="sd",conf=0.95)
	ordibar(inmat_PCoA, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(2,3),kind="sd",conf=0.05)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	# Now with spiders
	# First plot dimensions 1 and 2
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"))
	mtext(paste("Dimensions 1 and 2 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,1], inmat_PCoA$points[,2], col=s_col, pch=20,cex=1.1)
	ordispider(inmat_PCoA, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,2), lwd=0.5)
	legend("topleft", legend=names(presets$grp_less), pch=20, col=g_col,cex=1)
	# Now dimensions 1 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("Dim1 - Variation explained:",inmat_PCoA.axes[1],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"))
	mtext(paste("Dimensions 1 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,1], inmat_PCoA$points[,3], col=s_col, pch=20,cex=1.1)
	ordispider(inmat_PCoA, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,3), lwd=0.5)
	legend("topleft", legend=names(presets$grp_less), pch=20, col=g_col,cex=1)
	# Now dimensions 2 and 3
	ordiplot(inmat_PCoA,type="n", main=paste(paste0("PCoA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("Dim2 - Variation explained:",inmat_PCoA.axes[2],"%"), ylab=paste("Dim3 - Variation explained:",inmat_PCoA.axes[3],"%"))
	mtext(paste("Dimensions 2 and 3 - Anosim R:", compare[1,3], "p =", compare[1,5]))
	points(inmat_PCoA$points[,2], inmat_PCoA$points[,3], col=s_col, pch=20,cex=1.1)
	ordispider(inmat_PCoA, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(2,3), lwd=0.5)
	legend("topleft", legend=names(presets$grp_less), pch=20, col=g_col,cex=1)
	dev.off()
}

cluster_consistency <- function(in_matrix,presets) {
	all_labels <- names(presets$grp_less) # This is imported from outside
	beta_dist <- data.frame(matrix(NA,nrow=1,ncol=length(all_labels)*4)) # Create an empty matrix
	group_list <- vector(mode = "list", length = length(all_labels)*2) # as well as an empty list
	names(beta_dist) <- c(paste0(all_labels,"_w"),paste0(all_labels,"_b"),paste0(all_labels,"_pval"),paste0(all_labels,"_qval"))
	for(i in 1:length(all_labels)){
	# 	print(all_labels[i])
		list <- unlist(presets$grp_sam[i]) # get the list of the current group
		sub <- in_matrix[list,list] # subset the matrix to retain only those in the group
		within <- sub[lower.tri(sub)] # keep only the lower triangle (of a symmetric matrix)
		between <- unlist(in_matrix[list,setdiff((1:nrow(in_matrix)),list)], use.names = FALSE) # get the rows of the current group but the complement in the columns (all distances between the group and other samples)
		min_set <- min(length(between),length(within)) # group and non-group comparisons are normally uneven in numbers, keep track of the smallest
		group_list[[seq(1,length(group_list),2)[i]]] <- within # save the corresponding collections
		group_list[[seq(2,length(group_list),2)[i]]] <- between
		beta_dist[1,i] <- median(within) # and also the medians, which may vary a lot
		beta_dist[1,i+length(all_labels)] <- median(between)
		MW_test <- function() {wilcox.test(sample(within, min_set, replace=FALSE),sample(between, min_set, replace=FALSE))$p.value} # create a function to carry out a mann whitney test for pair, with a sample for n repetitions
		MW_all <- replicate(1000, MW_test()) # to save some time, calculate with 1000 items
		beta_dist[1,i+(length(all_labels)*2)] <- median(MW_all) # Save the medians for the pvalue
		beta_dist[1,i+(length(all_labels)*3)] <- median(p.adjust(MW_all, method="fdr")) # and its adjusted pvalue
		beta_dist
	}
	y <- trunc(seq(min(in_matrix),1,(1-min(in_matrix))/20)*1000)/1000
	pdf(paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"between_within_groups.pdf", sep="_"))
	boxplot(group_list, border=c("cornflowerblue","coral1"),col=c("darkslategray2","bisque"),main=paste("Between (b) and within (w) dissimilarities per group.",metric),frame=FALSE,outline=FALSE,yaxt='n',xaxt='n')
	mtext(prefix)
	axis(1,las=2,at=seq(1:(length(all_labels)*2)),labels=as.character(rbind(paste0(all_labels,"_w"), paste0(all_labels,"_b"))))
	axis(2,las=1,at=y)
	dev.off()
	write.table(cbind(prefix,beta_dist),paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"between_within_groups.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

####################################### MAIN #######################################
############ Part 1: Assess the group distribution, calculate stats ############
 ### Prefilter samples ###
library("vegan")
presets <- prepare_data(df) # Define how samples from each group are distributed
# Remove samples not in groups requested by the user:
df <- df[unlist(presets$grp_sam),] # and subset the matrix rowwise
df <- df[,unlist(presets$grp_sam)] # as well as columnwise
presets <- prepare_data(df) # Now recalculate presets
# Create color scheme (this will work if the total colors, defined above, is the correct one)
if(length(groups)<=length(group_cols)){ # Reassign colors to existing groups (to make them the same between graphs). This is only carried out if the number of defined colors (solid and border) is the same than the number of groups.
	cols <- set_cols(presets$grp_sam,presets$grp_all)
	s_col <- cols$cols
	g_col <- cols$group_cols
} else {
	s_col <- g_col <- "gray"
}
# Create a table with ad-hoc and post-hoc statistics
compare <- group_statistics(df, presets)
write.table(cbind(prefix,metric,compare),paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"stats.tsv", sep="_"), sep="\t", quote=FALSE, row.names=F) 
############ Part 2: Eigenvalue based MDS - Principal Coordinate Analysis ############
# Calculate Multidimensional scaling by calculating eigenvectors. Calculate each dimension's contribution and graph both with polygons and bullets
create_PCoA(df)
############ Part 3: Calculate between and within group dissimilarities ############
# Calculate the within and between median for each group
cluster_consistency(df,presets)
print("Done")
