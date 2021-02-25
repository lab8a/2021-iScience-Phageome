***Started on Thursday, February 28th, 2019 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

The data after quality filtering was blasted against the complete NR RefSeq viral database (downloaded in January 23 of 2019). After, we searched the lowest-common ancestor (LCA) with Meta Genome ANalyzer (MEGAN6). 

#### Creation of blast protein database with makeblastdb (blast version 2.6.0)

```bash
makeblastdb -in viral.all.protein.single-line.faa -dbtype prot -out viral_refseq -parse_seqids
```

#### Blastx of all samples against the NR RefSeq viral database (blast version 2.6.0)

```bash
for s in Final-All_ViralReads_*id.fa;
	do blastx -db viral_refseq -query $s -out $s.blastx -num_threads 64  -evalue 1e-3;
done
```

#### Parameters used in MEGAN6 (Community Edition, version 6_21_1)

- Min Score: 40.0
- Max Expected: 0.01
- Min Percent Identity: 0.0
- Top Percent: 10.0
- Min Support Percent: 0.05
- Min Support: 1
- Min Read Length: 0
- LCA Algorithm: naive
- Percent to cover: 100
- Read Assignment Mode: readCount







