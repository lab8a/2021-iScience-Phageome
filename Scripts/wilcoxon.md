***Started on Saturday, October 24th, 2020 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

This script was used to obtain multiple Wilcoxon test with a contingency table as input.

#### Creation of Rscript

The name of this file was: wilcoxon.R

```R
#!/usr/bin/env Rscript
## Purpose of script: Multiple Wilcoxon test
## Last modification: 24.10.2020
## Author: luigui gallardo-becerra (bfllg77@gmail.com)

# Args input
args = commandArgs(trailingOnly=TRUE)

# Package installation
#install.packages("matrixTests")
library(matrixTests)

# Data input
table <- read.delim(file = args[1],
    header = TRUE,
    row.names = 1)

h <- table[table$group=="H",2:ncol(table)]
o <- table[table$group=="O",2:ncol(table)]
omc <- table[table$group=="OMC",2:ncol(table)]

# Multiple Wilcoxon test
wilcoxon_h_o <- col_wilcoxon_twosample(h, o)
wilcoxon_h_omc <- col_wilcoxon_twosample(h, omc)
wilcoxon_o_omc <- col_wilcoxon_twosample(o, omc)

# Write table with results
write.table(wilcoxon_h_o,
    file = paste0(args[2],"_H-vs-O.txt"),
    sep = "\t"
)

write.table(wilcoxon_h_omc,
    file = paste0(args[2],"_H-vs-OMC.txt"),
    sep = "\t"
)

write.table(wilcoxon_o_omc,
    file = paste0(args[2],"_O-vs-OMC.txt"),
    sep = "\t"
)
```

#### Example of usage

```bash
Rscript wilcoxon.R input_table.txt output_prefix
```