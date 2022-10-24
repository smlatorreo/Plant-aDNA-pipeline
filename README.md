[![DOI](https://zenodo.org/badge/333746176.svg)](https://zenodo.org/badge/latestdoi/333746176)

This is a demonstration of the pipeline described in:  
**Latorre S.M., Lang, P.L.M., Burbano, H.A., Gutaker, R.M. 2020.  
Isolation, Library Preparation, and Bioinformatic Analysis of Historical and Ancient Plant DNA**  
d.o.i.: https://doi.org/10.1002/cppb.20121

To run the pipeline in your own machine, just download the file: plant\_aDNA\_pipeline.sh  
Then, simply type the following command and follow the instructions:
```
bash plant_aDNA_pipeline.sh --help
```

Alternatively, you can follow by copy-pasting the following demonstration steps

## SOFTWARE REQUIREMENTS:

Program           | Download Link
----------------- | --------------------------
AdapterRemoval v2 | (https://github.com/mikkelschubert/adapterremoval)
FastQC            | (https://github.com/s-andrews/FastQC)
BWA               | (https://github.com/lh3/bwa)
samtools          | (http://www.htslib.org/download/)
Dedup             | (https://github.com/apeltzer/DeDup)
MapDamage         | (https://github.com/ginolhac/mapDamage)

## ENVIRONMENT AND DATA ACQUISITION

1. Create the following folders:
* 1\_initial\_data/
* 1\_initial\_data/reference\_genome/
* 2\_trimmed\_merged/
* 3\_quality\_control/
* 4\_mapping/
* 5\_aDNA\_characteristics/

```
mkdir -p 1_initial_data/reference_genome 2_trimmed_merged 3_quality_control 4_mapping 5_aDNA_characteristics
```
	
2. We will use MiSeq pair-end Illumina sequences from a herbarium *Arabidopsis thaliana* sample (Weiss et al., 2016 ; https://doi.org/10.1098/rsos.160239) Download the PAIR 1 of reads to the folder 1_initial_data/

```
wget -P 1_initial_data ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR964/ERR964430/s_1_s8_R1.fastq.gz
```

3. Download the PAIR 2 of reads to the folder 1\_initial\_data/

```
wget -P 1_initial_data ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR964/ERR964430/s_1_s8_R2.fastq.gz
```
	
## TRIMMING AND MERGING
	
4. AdapterRemoval will be used to both remove Illumina adapters and merge paired reads when possible.
	
OPTIONS:  
--file1 1\_initial\_data/s\_1\_s8\_R1.fastq.gz : Location of PAIR 1 file  
--file2 1\_initial\_data/s\_1\_s8\_R2.fastq.gz : Location of PAIR 2 file  
--qualitybase 64 : Input Quality Score encoding is 64  
--qualitybase-output 33 : Output Quality Score encoding is 33 (standard)  
--basename 3\_trimmed\_merged/historicAthaliana : Location and prefix for the output files  
--collapse : Merge trimmed reads when possible  
--gzip : Output is compressed  
--threads $nthreads : Number of threads to be used

```	
AdapterRemoval --file1 1_initial_data/s_1_s8_R1.fastq.gz --file2 1_initial_data/s_1_s8_R2.fastq.gz --qualitybase 64 --qualitybase-output 33 --basename 2_trimmed_merged/historicAthaliana --collapse --gzip --threads <n_threads>
```
	
## QUALITY CONTROL

5. To assess the read quality of the sample, we will use the program FastQC.
	
OPTIONS:  
-o 3\_quality\_control/ : The results will be stored in 2\_QC/  
-t <n_threads> : Number of threads to be used  
2\_trimmed\_merged/historicAthaliana.collapsed.gz : Location of the trimmed and merged reads  
2\_trimmed\_merged/historicAthaliana.pair1.truncated.gz : Location of the trimmed PAIR 1 reads  
2\_trimmed\_merged/historicAthaliana.pair2.truncated.gz : Location of the trimmed PAIR 2 reads

```	
fastqc -o 3_quality_control/ -t <n_threads> 2_trimmed_merged/historicAthaliana.collapsed.gz 2_trimmed_merged/historicAthaliana.pair1.truncated.gz 2_trimmed_merged/historicAthaliana.pair2.truncated.gz
```

6. Inspect the quality control files

## MAPPING
	
7. Since the reads originate from an A. thaliana herbarium sample, the most suitable reference genome to use is the TAIR-10 *A. thaliana* assembly. The following command will download it and place it inside 1\_initial\_data/reference\_genome/

```	
curl ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz | zcat > 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```
	
8. The reference genome must be indexed both by BWA and samtools.

```
bwa index 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```
```
samtools faidx 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```

9. Use BWA aln to map the preprocessed (trimmed and merged) reads to the A. thaliana reference genome.
	
OPTIONS:
-l 1024 : Set the seed to a very high value (1024) to inactivate it  
-f mapping/historicAthaliana.collapsed.sai : Path and name for the output file  
-t $nthreads : Number of threads to be used  
1\_initial\_data/reference\_genome/Arabidopsis\_thaliana.TAIR10.dna.toplevel.fa : Location of the indexed reference genome  
2\_trimmed\_merged/historicAthaliana.collapsed.gz : Location of the trimmed and merged reads  
-t <n_threads> : Number of threads to be used

```	
bwa aln -l 1024 -f 4_mapping/historicAthaliana.collapsed.sai -t <n_threads> 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 2_trimmed_merged/historicAthaliana.collapsed.gz
```

10. Use BWA samse to convert the mapped reads (.sai file) into an alignment in the standard SAM format.
	
OPTIONS:  
-r \@\RG\\tID:sample1\\tSM:sample1 : Read Group tag as "sample1". This will be important for downstream analyses  
-f 4\_mapping/historicAthaliana.collapsed.sam : Location for the aligned sam file (output)  
1\_initial\_data/reference\_genome/Arabidopsis\_thaliana.TAIR10.dna.toplevel.fa : Location of the indexed reference genome  
4\_mapping/historicAthaliana.collapsed.sai : Location of the .sai file  
2\_trimmed\_merged/historicAthaliana.collapsed.gz : Location of the trimmed and merged reads

```	
bwa samse -r @RG\\tID:sample1\\tSM:sample1 -f 4_mapping/historicAthaliana.collapsed.sam 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 4_mapping/historicAthaliana.collapsed.sai 2_trimmed_merged/historicAthaliana.collapsed.gz'
```

11. Use samtools flagstat to estimate the endogenous DNA as the proportion of mapped reads.
	
OPTIONS:  
-@ <n_threads> : Number of threads to be used  
4\_mapping/historicAthaliana.collapsed.sam : Location of the alignment (input)  
\> 4\_mapping/historicAthaliana.collapsed.stats : Location for the output

```	
samtools flagstat -@ <n_threads> 4_mapping/historicAthaliana.collapsed.sam > 4_mapping/historicAthaliana.collapsed.stats
```
	
12. Since the SAM (Sequence Alignment Map) format is a plain text file, it can be efficiently compressed into a BAM (Binary Alignment Map) format using samtools. This allows to save disk space. Moreover, unmapped reads can be discarded.
	
OPTIONS:  
-@ <n_threads> : Number of threads to be used  
-Sb : The input is (S)AM and the output is (b)am  
-h : Include the header in the output  
-F 4 : Discard unmapped reads  
-o 4\_mapping/historicAthaliana.mapped.bam : Location for the output file  
4\_mapping/historicAthaliana.collapsed.sam : Location of the input file

```	
samtools view -@ <n_threads> -Sb -h -F 4 -o 4_mapping/historicAthaliana.mapped.bam 4_mapping/historicAthaliana.collapsed.sam
```

13. At this stage, it is possible to remove the uncompressed .sam file since the mapped reads are already compressed in the .bam file

```	
rm 4_mapping/historicAthaliana.collapsed.sam
```
	
14. Most downstream analyses require coordinate-based sorted files. Use samtools to sort the mapped bam file.
	
OPTIONS:  
-@ <n_threads> :  Number of threads to be used  
-o 4\_mapping/historicAthaliana.mapped.sorted.bam : Location for the sorted output  
4\_mapping/historicAthaliana.mapped.bam : Location of the unsorted mapped reads (input)

```	
samtools sort -@ <n_threads> -o 4_mapping/historicAthaliana.mapped.sorted.bam 4_mapping/historicAthaliana.mapped.bam
```
	
15. The unsorted .bam file can be removed

```	
rm 4_mapping/historicAthaliana.mapped.bam
```

# DUPLICATES REMOVAL

16. For downstream analyses it is important to remove the optical PCR duplicates.
	
OPTIONS:  
-i 4\_mapping/historicAthaliana.mapped.sorted.bam : Location of the mapped and sorted file (input)  
-m : The input contains merged reads  
-o 4\_mapping/ : Location of the output folder

```	
java -jar deDup.jar -i 4_mapping/historicAthaliana.mapped.sorted.bam -m -o 4_mapping/
```

# ANCIENT DNA PROPERTIES ASSESSMENT

17. Use mapDamage2 to assess the aDNA damage on the BAM file	
	
OPTIONS:  
-i 4\_mapping/historicAthaliana.mapped.sorted\_rmdup.bam : Location of the input file  
-r 1\_initial\_data/reference\_genome/Arabidopsis\_thaliana.TAIR10.dna.toplevel.fa : Location of the reference genome  
--no-stats : Disable the statistical calculation  
-y 0.05 : Plot Y-axis for nucleotide misincorporation frequencies plot up to 0.05  
-d 5\_aDNA\_characteristics/ : Location of the output folder

```	
mapDamage -i 4_mapping/historicAthaliana.mapped.sorted_rmdup.bam -r 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --no-stats -y 0.05 -d 5_aDNA_characteristics
```

18. Inspect the generated files
