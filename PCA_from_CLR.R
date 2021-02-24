# Started on 2021-02-16
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# The script is intended to plot PCA ordinations using a contingency (incidence) matrix. It was based on a similar algorithm for distance/correlation matrices called PCA_from_dm.R and adapted to use the rectangular matrix directly
# Groups are defined by the user (colors are defined in-script)
# The input is a rectangular (non-squared) matrix, with columns as sites, rows as "species" according to vegan's nomenclature (transposed with respect to the expected input)

# Run as follows:
# cat dm.tsv|Rscript PCA_from_CLR.R <prefix_output>  <name_of_input_metric> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>
# Tested with command:
# cat /home/rod/Documents/02_Collaborations/Phageome_obesity/Round2/CLR/phage-rar_CLR-10000-depth-20000-perm.tsv|Rscript PCA_from_CLR.R test CLR H O- OMC
# In R:
# df <- read.table("CLR/phage-rar_CLR-10000-depth-20000-perm.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
# metric="CLR"
# prefix="CLR/phage-rar_CLR-10000-depth-20000-perm"
# groups <- c("H","O-","OMC")

# Input table
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, two arguments should be included: <min_nonNAs> <prefix_output>  <name_of_alpha_metric>
  stop("A minimum of 4 arguments are mandatory: cat table.tsv|Rscript PCA_from_dm.R <prefix_output>  <name_of_beta_metric> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>", call.=FALSE)
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
 group_cols <- c("blue","red","green4","black","magenta") # I've used this for runs (5 grps)
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
single_comparison <- function(inmat, grp_factors){ # for carrying out single comparisons with a given submatrix
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
create_PCA <- function(mat){ # This function plots the ordination plots for the first 3 dimensions
	vpca_20k = rda(t(mat),scale=TRUE)
	vpca_20k.axes <- round(100*vpca_20k$CA$eig/sum(vpca_20k$CA$eig),2)
	cols=c(rep(4,8),rep(2,10),rep(3,8))
	r1 <- range(scores(vpca_20k,choices=c(1))$sites) #original ranges
	r2 <- range(scores(vpca_20k,choices=c(2))$sites)
	r3 <- range(scores(vpca_20k,choices=c(3))$sites)
	fix=0.05
	dim1 <- c(r1[1]-diff(r1)*fix,r1[2]+diff(r1)*fix) # Fixed ranges (adjust fix as desired
	dim2 <- c(r2[1]-diff(r2)*fix,r2[2]+diff(r2)*fix)
	dim3 <- c(r3[1]-diff(r3)*fix,r3[2]+diff(r3)*fix)
	# First plot dimensions 1 and 2
	pdf(paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"PCA.pdf", sep="_"))
	par(las=1)
	plot(round(cumsum(100*vpca_20k$CA$eig/sum(vpca_20k$CA$eig)),2), ylim=c(0,100), type='o', main="Variation explained by Principal Components", xlab="Principal component", ylab="Cummulative variation explained")
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"),xlim=dim1,ylim=r2) # Create empty plot
	orditorp(vpca_20k,display="sites",cex=0.8,air=0.01,col=s_col, choices=c(1,2)) # add items
	ordihull(vpca_20k, groups=as.factor(s_col), draw="polygon", col=levels(as.factor(s_col)),  border=levels(as.factor(s_col)), label=F, choices=c(1,2),alpha=50) # add polygons
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Next, dimensions 1 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"),xlim=dim1,ylim=r3)
	orditorp(vpca_20k,display="sites",cex=0.8,air=0.01,col=s_col, choices=c(1,3))
	ordihull(vpca_20k, groups=as.factor(s_col), draw="polygon", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,3),alpha=50)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Finally, dimensions 2 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"),xlim=dim2,ylim=r3)
	orditorp(vpca_20k,display="sites",cex=0.8,air=0.01,col=s_col, choices=c(2,3))
	ordihull(vpca_20k, groups=as.factor(s_col), draw="polygon", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(2,3),alpha=50)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now with bullets only
	# First plot dimensions 1 and 2
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"))
	points(scores(vpca_20k,choices=c(1))$sites, scores(vpca_20k,choices=c(2))$sites, col=s_col, pch=1,cex=1.1)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now dimensions 1 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"))
	points(scores(vpca_20k,choices=c(1))$sites, scores(vpca_20k,choices=c(3))$sites, col=s_col, pch=1,cex=1.1)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now dimensions 2 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"))
	points(scores(vpca_20k,choices=c(2))$sites, scores(vpca_20k,choices=c(3))$sites, col=s_col, pch=1,cex=1.1)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now with ellipses (sd with conf 0.95 over chi2 distr)
	# First plot dimensions 1 and 2
	fix=0.35
	dim1 <- c(r1[1]-diff(r1)*fix,r1[2]+diff(r1)*fix)
	dim2 <- c(r2[1]-diff(r2)*fix,r2[2]+diff(r2)*fix)
	dim3 <- c(r3[1]-diff(r3)*fix,r3[2]+diff(r3)*fix)
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"),xlim=dim1,ylim=dim2)
	points(scores(vpca_20k,choices=c(1))$sites, scores(vpca_20k,choices=c(2))$sites, col=s_col, pch=1,cex=1.1)
	ordiellipse(vpca_20k, groups=as.factor(s_col), draw="lines", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,2),kind="sd",conf=0.95)
	ordibar(vpca_20k, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,2),kind="sd",conf=0.05)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now dimensions 1 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"),xlim=dim1,ylim=dim3)
	points(scores(vpca_20k,choices=c(1))$sites, scores(vpca_20k,choices=c(3))$sites, col=s_col, pch=1,cex=1.1)
	ordiellipse(vpca_20k, groups=as.factor(s_col), draw="lines", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(1,3),kind="sd",conf=0.95)
	ordibar(vpca_20k, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,3),kind="sd",conf=0.05) # get the group centroid
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now dimensions 2 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"),xlim=dim2,ylim=r3)
	points(scores(vpca_20k,choices=c(2))$sites, scores(vpca_20k,choices=c(3))$sites, col=s_col, pch=1,cex=1.1)
	ordiellipse(vpca_20k, groups=as.factor(s_col), draw="lines", col=levels(as.factor(s_col)), border=levels(as.factor(s_col)), label=F, choices=c(2,3),kind="sd",conf=0.95)
	ordibar(vpca_20k, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(2,3),kind="sd",conf=0.05)
	legend("topleft", legend=names(presets$grp_less), pch=1, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now with spiders
	# First plot dimensions 1 and 2
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,2),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"))
	points(scores(vpca_20k,choices=c(1))$sites, scores(vpca_20k,choices=c(2))$sites, col=s_col, pch=20,cex=1.1)
	ordispider(vpca_20k, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,2), lwd=0.5)
	legend("topleft", legend=names(presets$grp_less), pch=20, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now dimensions 1 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(1,3),cex.main=.7,xlab=paste("PC1 - Variation explained:",vpca_20k.axes[1],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"))
	points(scores(vpca_20k,choices=c(1))$sites, scores(vpca_20k,choices=c(3))$sites, col=s_col, pch=20,cex=1.1)
	ordispider(vpca_20k, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(1,3), lwd=0.5)
	legend("topleft", legend=names(presets$grp_less), pch=20, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	# Now dimensions 2 and 3
	ordiplot(vpca_20k,type="n", main=paste(paste0("PCA -",metric),prefix,sep="\n"), choices=c(2,3),cex.main=.7,xlab=paste("PC2 - Variation explained:",vpca_20k.axes[2],"%"), ylab=paste("PC3 - Variation explained:",vpca_20k.axes[3],"%"))
	points(scores(vpca_20k,choices=c(2))$sites, scores(vpca_20k,choices=c(3))$sites, col=s_col, pch=20,cex=1.1)
	ordispider(vpca_20k, groups=as.factor(s_col), col=levels(as.factor(s_col)), label=F, choices=c(2,3), lwd=0.5)
	legend("topleft", legend=names(presets$grp_less), pch=20, col=g_col,cex=1)
	mtext(paste("Anosim R:", compare[1,3], "p =", compare[1,5]," / Adonis R2:", compare[1,2], "p =", compare[1,4]))
	dev.off()
}

####################################### MAIN #######################################
############ Part 1: Assess the group distribution, calculate stats ############
 ### Prefilter samples ###
library("vegan")
presets <- prepare_data(df) # Define how samples from each group are distributed
# Remove samples not in groups requested by the user:
df <- df[,unlist(presets$grp_sam)] # and subset the matrix columnwise
presets <- prepare_data(df) # Now recalculate presets
# Create color scheme (this will work if the total colors, defined above, is the correct one)
if(length(groups)<=length(group_cols)){ # Reassign colors to existing groups (to make them the same between graphs). This is only carried out if the number of defined colors (solid and border) is the same than the number of groups.
	cols <- set_cols(presets$grp_sam,presets$grp_all)
	s_col <- cols$cols
	g_col <- cols$group_cols
} else {
	s_col <- g_col <- "gray"
}
eucl <- dist(t(df), method="euclidean", upper=T, diag=T) # Produce an euclidean distance matrix from the input data
eucl.sim <- as.matrix(eucl/max(eucl)) # Change it to 0-1 similarity (where 1 is max sim)
write.table(as.matrix(eucl.sim),paste(prefix,metric,"eucl_mat.tsv", sep="_"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # print the output euclidean matrix
compare <- group_statistics(eucl.sim, presets)
write.table(cbind(prefix,metric,compare),paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"stats.tsv", sep="_"), sep="\t", quote=FALSE, row.names=F) 
############ Part 2: Principal Components Analysis ############
# Carry out multidimensional scaling by calculating eigenvectors. Calculate each dimension's contribution and graph both with polygons and bullets
create_PCA(df)
print("Done")
