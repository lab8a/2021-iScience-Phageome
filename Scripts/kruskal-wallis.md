***Started on Saturday, October 24th, 2020 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

This script was used to obtain multiple Kruskal-Wallis rank sum test with a contingency table as input.

#### Creation of Rscript

The name of this file was: kruskal-wallis.R

```R
#!/usr/bin/env Rscript
## Purpose of script: Multiple Kruskal-Wallis rank sum test
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

# Multiple Kruskal-Wallis rank sum test
kruska_wallis_result <- col_kruskalwallis(table[2:ncol(table)], table$group)

# Write table with results
write.table(kruska_wallis_result,
    file = args[2],
    sep = "\t"
)
```

#### Example of usage

```bash
Rscript kruskal-wallis.R input_table.txt kruskal-wallis_result.txt
```

