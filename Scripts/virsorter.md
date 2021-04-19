***Started on Thursday, march 4th, 2021 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

For this analysis we clone the latest version available of virsorter (2.2.1) from GitHub (https://github.com/jiarong/VirSorter2). After, both the quality reads and the assembly were classificated with virsorter following the online tutorial.


#### Instalation virsorter (version 2.2.1)

```bash
conda create -n vs2 -c conda-forge -c bioconda "python>=3.6" scikit-learn=0.22.1 imbalanced-learn pandas seaborn hmmer==3.3 prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click mamba

conda activate vs2

git clone https://github.com/jiarong/VirSorter2.git

cd VirSorter2

pip install -e .
```
#### Classification of reads

```bash
for s in Final-All_ViralReads_*id.fa;
	do virsorter run -w virsorter_$s -i $s -j 64;
done
```

#### Classification of assembly

```bash
virsorter run -w virsorter_FinalViralScaffolds_larger4Kb -i FinalViralScaffolds_larger4Kb.fasta -j 64
```