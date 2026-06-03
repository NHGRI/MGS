
# Description
This is the core pipeline to perform an allelic imbalance analsis of ATAC-seq reads at heterozygous SNP sites in one more subjects.  

A sample set of data is included in the data/example_data. Note that this data differs different from actual data used in the manuscrtip. Two types of files are needed in minimum to run this pipeline. They are described in the next section.

# Prerequisites
There are a couple of prerequisite input files: 1) subject level bam file of paired-end raw ATAC-seq sequence reads to GRCh38 reference genome; and 2) vcf file of imputed genotypes of the subjects using a tool like TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov). 

# Components
This workflow performs following tasks: 1) identification of ATAC-seq reads overlapping a SNP in a subject; 2) correction of reference allele mapping bias for each subject where a read overlaps a SNP alternate allele; 3) generation of read counts mapping to reference and alternate alleles across all SNPs for each subject; 4) beta-binomial test of allelic imbalance of ATAC-reads at heterozygous sites for each subject; and 5) meta-analysis of beta-binomial test results at each SNP site across all subjects; 6) multiple testing correction on the beta-binomial test results across all tested SNP sites.

# Workflow
### 1. Software needed
1. WASP (https://github.com/bmvdgeijn/WASP)
2. Samtools (https://www.htslib.org)
3. bamUtil (https://github.com/statgen/bamUtil)
   
### 2. Quick start
ATAC-seq allelic imbalance analysis:
```
snakemake  --snakefile pipeline/Snakemake_wasp_celltype_bb --cluster-config configs/config_cluster_sge.json --configfile configs/config.yaml  --cluster configs/config_slurm.py --printshellcmds --latency-wait 600 --jobs 250 --rerun-incomplete --scheduler greedy -k
```

### 3. Output file
This pipeline will generate a file called "ATAC_peak_SNPs_stouffer_allelic_imbalance_betabinom.txt". This is a bed style file with following columns:   
   1. #chr   
   2. start   
   3. stop   
   4. var_id   
   5. ref_base   
   6. alt_base   
   7. meta_betabinomial_p: meta-analysis beta binomial p value across samples with heterozygous genotypes for the SNP.   
   8. meta_betabinomial_q: BH adjusted meta-analysis beta binomial q value. We consider a SNP with this value <0.05 to be imbalanced.   
   9. pseudo.ref.count: total number of ATAC-seq reads with reference allele across samples.   
   10. pseudo.total.count: total number of ATAC-seq reads with reference and alternate alleles across samples.   
   11. pseudo_betabinomial_p: beta binomial test for read counts of reference and alternate alleles across all samples with hets.   
   12. pseudo_betabinomial_q: BH adjusted beta binomial q value for pseudo bulk read counts.  
   13. pseudo_betabinomial_p_flag: The same as pseudo_betabinomial_p. In extremely rare cases they are different due to rounding errors.
   14. ref.effect: reads_covering_reference_allele/reads_covering_both_alleles - 0.5   

### 3. Reference


