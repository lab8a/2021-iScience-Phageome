[TOC]

# Gut Virome Shows Diversity and Richness Alterations Associated to Childhood Obesity and Metabolic Syndrome



***Protocol written by Gamaliel López-Leal***

***gamlopez@ccg.unam.mx***



### Raw reads filtered by quality###

```
fastx_trimmer -v -f 20  -i Reads.fastq  -o Reads_trimmed_First_20nt.fastq #this line remove the adpaters 

trim_galore -q 30  --paired    R1-Reads_trimmed_First20nt.fastq    R2-Reads_trimmed_First20nt.fastq  # this line filter the reads by quality using a quality threshold of ≥ 30 (phred quality scores)
```



### Revoming human and bacterial reads ###

First we removed all reads that contained kmers associated to any bacterial taxonomic node using Kraken

```
kraken --db krakenDB/bacteria_RefSeq --fastq-input --paired quality-reads-R1.fq quality-reads-R2.fq --threads 10  --unclassified-out Unclassified-Reads  --classified-out Classified-Reads --output out.wdir
```

Then, the Unclassified-Reads were filtered by read mapping against *Homo sapines* genome v38

```
bwa mem Hs_Genome-v38.fasta -t 7 -p interleaved-Unclassified-Reads.fq  > Hs_aln_unclass.sam
```

We extract the unmapped reads by following the next command line

```
samtools view -b -S -f 4 Hs_aln_unclass.sam > ViralReads.bam # extract unmapped reads

bamToFastq -i ViralReads.bam -fq ViralReads.bam.fq # convert the the information in the .bam to a .fastq

grep -c "@N" ViralReads.fq > Ids_ViralReads_firts.txt # get all the unmapped ids-reads

sort -u Ids_ViralReads_redundant.txt > Ids_ViralReads_final.txt 

seqtk subseq Unclassified-Reads.fq Ids_ViralReads_final.txt  > Final-ViralReads.fq #using the unmapped ids-reads we collected the reads in pair mode

```

### *De Novo* viral assembly ###

Total viral reads from all samples were pooled for *de novo* assembly using IDBA-UD using the following command line

```
idba_ud -r Interleave_ViralReads.fasta --num_threads 8 --mink 20 --maxk 125 --pre_correction -o out.dir
```

Then we evaluated the coverage of each contig following the command lines

```
bowtie2-build ViralContigs.fa ViralContigs # index the viral contigs

bowtie2 --no-unal  --sensitive  -x ViralContigs -1 ViralReads-for-eacch-sample1.fastq -2 ViralReads-for-eacch-sample2.fastq -S bowtie.sam -p 8 # read mapping

samtools view -bS bowtie.sam -o  bowtie.bam 
samtools sort bowtie.bam bowtie_sorted.bam

samtools depth bowtie_sorted.bam > aln-ViralContigs.coverage # get the depth coverage for each viral contig

samtools idxstats bowtie_sorted.bam > ReadCount_ViralContigs.txt # get the contig abundance

```










































