***Started on Saturday, October 24th, 2020 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

This script was used to obtain multiple Spearman's correlations with a contingency table as input.

#### Creation of Rscript

The name of this file was: spearman_correlations.R

```R
#!/usr/bin/env Rscript
## Purpose of script: Spearman's correlations
## Last modification: 24.10.2020
## Author: luigui gallardo-becerra (bfllg77@gmail.com)

# Args input
args = commandArgs(trailingOnly=TRUE)

# Package installation
#install.packages("Hmisc")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("reshape")
#install.packages("rlist")
#install.packages("igraph")
#install.packages("ggpubr")
library(Hmisc)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rlist)
library(igraph)
library(ggpubr)

# Data input
table <- read.delim(file = args[1],
    header = TRUE,
    row.names = 1)
table_transposed <- t(table)

# Creation of correlation matrix
# Parameters of Spearman's correlations
cor.cutoff = 0.3 # This is the R value, could be changed
p.cutoff = 0.05 # This is the p-value, could be changed

# Correlation matrix
res.cor <- rcorr(t(table_transposed),
    type="spearman")
matrix.cor <- res.cor$r # This is the R value matrix
matrix.cor.p <- res.cor$P # This is the p-value matrix

# Filter of R and p values
matrix.cor[which(matrix.cor >= (-cor.cutoff) & matrix.cor <= cor.cutoff)] = 0 # Filter of R (positive and negative)
matrix.cor[which(matrix.cor.p >= p.cutoff)] = 0 # Filter of p-value

# Reshape the correlation matrix
melted_cormat <- melt(matrix.cor)
melted_cormat_pvalue <- melt(matrix.cor.p)

# Final filter: remove zeros, self-correlations and order taxa and function
final_cormat <- melted_cormat %>% # Correlation output
#    filter(value != 0) %>% 0 # Filter 0
    filter(Var1 != Var2) %>% # Remove auto-correlations 
    filter(!grepl("scaffold", Var1)) %>% # Remove scaffolds in Var1
    filter(grepl("k__", Var1)) %>% # Mantain taxonomy in Var1
    filter(grepl("scaffold", Var2)) # Mantain scaffolds in Var2

final_cormat_pvalue <- melted_cormat_pvalue %>% # p-value output
#    filter(value != 0) %>% 0 # Filter 0
    filter(Var1 != Var2) %>% # Remove auto-correlations 
    filter(!grepl("scaffold", Var1)) %>% # Remove scaffolds in Var1
    filter(grepl("k__", Var1)) %>% # Mantain taxonomy in Var1
    filter(grepl("scaffold", Var2)) # Mantain scaffolds in Var2

# Write table with Rho
write.table(final_cormat,
    file = args[2],
    sep = "\t"
)

# Write table with p-value
write.table(final_cormat_pvalue,
    file = args[3],
    sep = "\t"
)
```

#### Example of usage

```bash
Rscript spearman_correlations.R input_table.txt spearman_rho.txt spearman_p-value.txt
```