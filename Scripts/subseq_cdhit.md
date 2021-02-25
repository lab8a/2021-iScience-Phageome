***Started on Thursday, April 4th, 2019 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

The input for this analysis was the data after quality filtering. Random selection was adjusted to the maximum of the smallest sample (149 000 sequences). First we used *seqtk sample* to get the random subsequences per sample and after we clustered each of the 1 000 files per sample with *cd-hit*. Finally, the mean per sample was calculated with *awk*.

#### Random selection of 149 000 sequences per sample with setqk (setqk version 1.3-r106)

```bash
for s in {1..1000};
	do seqtk sample -s$s Final-All_ViralReads_H-10_R1.fq 149000 > sub_$s\_H-10_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-118_R1.fq 149000 > sub_$s\_H-118_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-119_R1.fq 149000 > sub_$s\_H-119_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-120_R1.fq 149000 > sub_$s\_H-120_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-124_R1.fq 149000 > sub_$s\_H-124_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-147_R1.fq 149000 > sub_$s\_H-147_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-161_R1.fq 149000 > sub_$s\_H-161_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-169_R1.fq 149000 > sub_$s\_H-169_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-193_R1.fq 149000 > sub_$s\_H-193_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_H-314_R1.fq 149000 > sub_$s\_H-314_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-121_R1.fq 149000 > sub_$s\_O-121_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-122_R1.fq 149000 > sub_$s\_O-122_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-123_R1.fq 149000 > sub_$s\_O-123_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-152_R1.fq 149000 > sub_$s\_O-152_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-39_R1.fq 149000 > sub_$s\_O-39_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-418_R1.fq 149000 > sub_$s\_O-418_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-420_R1.fq 149000 > sub_$s\_O-420_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-434_R1.fq 149000 > sub_$s\_O-434_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-445_R1.fq 149000 > sub_$s\_O-445_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_O-90_R1.fq 149000 > sub_$s\_O-90_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-124_R1.fq 149000 > sub_$s\_OMC-124_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-125_R1.fq 149000 > sub_$s\_OMC-125_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-126_R1.fq 149000 > sub_$s\_OMC-126_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-288_R1.fq 149000 > sub_$s\_OMC-288_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-446_R1.fq 149000 > sub_$s\_OMC-446_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-55_R1.fq 149000 > sub_$s\_OMC-55_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-64_R1.fq 149000 > sub_$s\_OMC-64_149k_R1.fq;
	seqtk sample -s$s Final-All_ViralReads_OMC-87_R1.fq 149000 > sub_$s\_OMC-87_149k_R1.fq;
done
```

#### Cluster generation with CD-HIT (indentity=95%, CD-HIT version 4.6)

```bash
for s in sub_*fq ;
	do cd-hit -i $s -o $s.cd-hit_95.id.fq -c 0.95 -M 0 -T 64 > $s.out;
done
```

#### Get the mean of reads per sample (awk version 4.1.4)

```bash
printf 'sample\tcount\n' > readcount_seqtk_cdhit.txt

printf 'H-10\t' >> readcount_seqtk_cdhit.txt
for s in *H-10*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$1}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-118\t' >> readcount_seqtk_cdhit.txt
for s in *H-118*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-119\t' >> readcount_seqtk_cdhit.txt
for s in *H-119*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-120\t' >> readcount_seqtk_cdhit.txt
for s in *H-120*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-124\t' >> readcount_seqtk_cdhit.txt
for s in *H-124*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-147\t' >> readcount_seqtk_cdhit.txt
for s in *H-147*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-161\t' >> readcount_seqtk_cdhit.txt
for s in *H-161*cd-hit_95.id.fq ;	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-169\t' >> readcount_seqtk_cdhit.txt
for s in *H-169*cd-hit_95.id.fq ;	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-193\t' >> readcount_seqtk_cdhit.txt
for s in *H-193*cd-hit_95.id.fq ;	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'H-314\t' >> readcount_seqtk_cdhit.txt
for s in *H-314*cd-hit_95.id.fq ;	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-121\t' >> readcount_seqtk_cdhit.txt
for s in *O-121*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-122\t' >> readcount_seqtk_cdhit.txt
for s in *O-122*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-123\t' >> readcount_seqtk_cdhit.txt
for s in *O-123*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-152\t' >> readcount_seqtk_cdhit.txt
for s in *O-152*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-39\t' >> readcount_seqtk_cdhit.txt
for s in *O-39*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-418\t' >> readcount_seqtk_cdhit.txt
for s in *O-418*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-420\t' >> readcount_seqtk_cdhit.txt
for s in *O-420*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-434\t' >> readcount_seqtk_cdhit.txt
for s in *O-434*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-445\t' >> readcount_seqtk_cdhit.txt
for s in *O-445*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'O-90\t' >> readcount_seqtk_cdhit.txt
for s in *O-90*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-124\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-124*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-125\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-125*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-126\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-126*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-288\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-288*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-446\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-446*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-55\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-55*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-64\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-64*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$s}END{print sum/1000}' >> readcount_seqtk_cdhit.txt

printf 'OMC-87\t' >> readcount_seqtk_cdhit.txt
for s in *OMC-87*cd-hit_95.id.fq ;
	do grep -c '^@NS' $s ; 
done | awk '{sum =+$1}END{print sum/1000}' >> readcount_seqtk_cdhit.txt
```


