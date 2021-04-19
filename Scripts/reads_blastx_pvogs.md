***Started on Thursday, march 18th, 2021 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

The quality-reads were mapped against the Prokaryotic Virus Orthologous Groups (pVOGs) database using BLASTX.

#### Creation of blast protein database with makeblastdb (blast version 2.6.0)

```bash
makeblastdb -in pVOGs.faa -dbtype prot -out pVOGs -parse_seqids
```

#### Blastx of all samples against the NR RefSeq viral database (blast version 2.6.0)

```bash
for s in Final-All_ViralReads_*id.fa;
	do blastx -db pVOGs -query $s -out $s\pvogs.blastx -num_threads 64 -evalue 1e-3 -max_target_seqs 50 -outfmt 7;
done
```
