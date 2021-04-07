***Started on Monday, February 18th, 2019 by Luigui Gallardo-Becerra* for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.**

***Disclaimer: These commands were used for specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.***

For this analysis, we use amplicons of the 16S rRNA (V4 hypervariable region). The QIIME (version 1.9.1) input was the data after quality filtering, that included remotion of barcodes, scanning of the read with a 4-base wide sliding window with a minimum quality of Q20 and remotion of N bases. This pretreatment was performed with Trimmomatic (version 0.36).

#### Pick closed reference, Greengenes database (version gg_13_8_otus, 97_otus) and filtering low abundance OTUs

```bash
# Pick OTUs
parallel_pick_otus_uclust_ref.py -i amplicons_16s_v4.fasta -o pick_closed_gg13_8/uclust_ref_picked_otus -r /a1/luigui/.local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta -T --jobs_to_start 10 --max_rejects 500 --stepwords 20 --enable_rev_strand_match --word_length 12 --max_accepts 20

# Make OTU table
make_otu_table.py -i pick_closed_gg13_8/uclust_ref_picked_otus/amplicons_16s_v4.fasta_otus.txt -t /a1/luigui/.local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -o pick_closed_gg13_8/otu_table.biom 

cd pick_closed_gg13_8/

# Filtering low abundance OTUs
filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_0.005.biom --min_count_fraction 0.00005
```

#### Summarize taxonomy

```bash
# Collapse samples in OTU table by categories
collapse_samples.py -m metadata.txt -b otu_table_0.005.biom --output_biom_fp sum_taxa_0.005_group/Group_otu_table.biom --output_mapping_fp sum_taxa_0.005_group/Group_map.txt --collapse_fields 'Group'

# Summarize Taxonomy
summarize_taxa.py -i sum_taxa_0.005_group/Group_otu_table.biom -o sum_taxa_0.005_group/ 

# Plot Taxonomy Summary
plot_taxa_summary.py -i sum_taxa_0.005_group/Group_otu_table_L2.txt,sum_taxa_0.005_group/Group_otu_table_L3.txt,sum_taxa_0.005_group/Group_otu_table_L4.txt,sum_taxa_0.005_group/Group_otu_table_L5.txt,sum_taxa_0.005_group/Group_otu_table_L6.txt -o sum_taxa_0.005_group//taxa_summary_plots/ 
```

#### Alpha diversity (Shannon, Chao1, Goods Coverage, Observed OTUs and Phylogenetic Diversity whole tree)

```bash
# Alpha rarefaction
parallel_multiple_rarefactions.py -T -i otu_table_0.005.biom -o arare_m31891_r10000//rarefaction/ --max 31891 --step 1 --num_reps 10000 --min 31891 --jobs_to_start 64

# Alpha diversity on rarefied OTU tables
parallel_alpha_diversity.py -T -i arare_m31891_r10000//rarefaction/ -o arare_m31891_r10000//alpha_div/ --metrics shannon,chao1,goods_coverage,observed_otus,PD_whole_tree -t 97_otus.tree --jobs_to_start 64

# Collate alpha
collate_alpha.py -i arare_m31891_r10000//alpha_div/ -o arare_m31891_r10000//alpha_div_collated/

# Rarefaction plot: All metrics
make_rarefaction_plots.py -i arare_m31891_r10000//alpha_div_collated/ -m mapping_file_all_21082018.txt -o arare_m31891_r10000//alpha_rarefaction_plots/
```

#### Beta diversity (Weighted-UniFrac, Unweighted-UniFrac and Bray-Curtis)

```bash
# Weighted Unifrac
# Sample OTU table at 31891 seqs/sample
single_rarefaction.py -i otu_table_0.005.biom -o bdiv_0.005_m31891/otu_table_0.005_even31891.biom -d 31891

# Beta Diversity (weighted_unifrac)
parallel_beta_diversity.py -i bdiv_0.005_m31891/otu_table_0.005_even31891.biom -o bdiv_0.005_m31891 --metrics weighted_unifrac -T  -t 97_otus.tree --jobs_to_start 64

# Rename distance matrix (weighted_unifrac)
mv bdiv_0.005_m31891/weighted_unifrac_otu_table_0.005_even31891.txt bdiv_0.005_m31891/weighted_unifrac_dm.txt

# Principal coordinates (weighted_unifrac)
principal_coordinates.py -i bdiv_0.005_m31891/weighted_unifrac_dm.txt -o bdiv_0.005_m31891/weighted_unifrac_pc.txt 

# Weighted Unifrac
# Beta Diversity (unweighted_unifrac)
parallel_beta_diversity.py -i bdiv_0.005_m31891/otu_table_0.005_even31891.biom -o bdiv_0.005_m31891 --metrics unweighted_unifrac -T  -t 97_otus.tree --jobs_to_start 64

# Rename distance matrix (unweighted_unifrac)
mv bdiv_0.005_m31891/unweighted_unifrac_otu_table_0.005_even31891.txt bdiv_0.005_m31891/unweighted_unifrac_dm.txt

# Principal coordinates (unweighted_unifrac)
principal_coordinates.py -i bdiv_0.005_m31891/unweighted_unifrac_dm.txt -o bdiv_0.005_m31891/unweighted_unifrac_pc.txt

# Bray-Curtis
# Sample OTU table at 31891 seqs/sample
single_rarefaction.py -i otu_table_0.005.biom -o bdiv_0.005_m31891_braycurtis/otu_table_0.005_even31891.biom -d 31891

# Beta Diversity (bray_curtis) 
parallel_beta_diversity.py -i otu_table_0.005.biom -o bdiv_0.005_m31891_braycurtis --metrics bray_curtis -T --jobs_to_start 64

# Rename distance matrix (bray_curtis)
mv bdiv_0.005_m31891_braycurtis/bray_curtis_otu_table_0.005.txt bdiv_0.005_m31891_braycurtis/bray_curtis_dm.txt

# Principal coordinates (bray_curtis)
principal_coordinates.py -i bdiv_0.005_m31891_braycurtis/bray_curtis_dm.txt -o bdiv_0.005_m31891_braycurtis/bray_curtis_pc.txt
```





