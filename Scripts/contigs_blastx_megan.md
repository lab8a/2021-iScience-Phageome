***Started on Monday, March 4th, 2019 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

The filtered contigs were blasted against the complete NR RefSeq viral database (downloaded in January 23 of 2019). After, we searched the lowest-common ancestor (LCA) with Meta Genome ANalyzer (MEGAN6). 

#### Creation of blast nucleotide and protein database with makeblastdb (blast version 2.6.0)

```bash
makeblastdb -in viral_refseq_nucl.fasta -dbtype db_nucl_refseq.faa -out db_nucl_refseq -parse_seqids

makeblastdb -in viral.all.protein.single-line.faa -dbtype prot -out viral_refseq -parse_seqids
```

#### Dc-Megablast of all samples against the NR RefSeq viral database (blast version 2.6.0)

```bash
blastn -task dc-megablast -db db_nucl_refseq -query FinalViralScaffolds_larger4Kb.fasta -out FinalViralScaffolds_larger4Kb.megablast.blastn -template_type coding_and_optimal -template_length 16 -evalue 1e-3 -num_threads 64
```

#### Blastx of all samples against the NR RefSeq viral database (blast version 2.6.0)

```bash
blastx -db viral_refseq -query inalViralScaffolds_larger4Kb.fasta -out FinalViralScaffolds_larger4Kb.fasta.blastx -evalue 1e-3 -num_threads 64
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







