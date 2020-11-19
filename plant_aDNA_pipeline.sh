#! /bin/bash

header() {
	clear
	echo "
################################################################################
########   This is a demonstration of the pipeline described in:        ########
########                                                                ########
########  Latorre S.M., Lang, P.L.M., Burbano, H.A., Gutaker, R.M. 2020 ########
########   Isolation, Library Preparation, and Bioinformatic Analysis   ########
########               of Historical and Ancient Plant DNA              ########
########           d.o.i.: https://doi.org/10.1002/cppb.20121           ########
################################################################################
"
}

help_message() {
	echo -e "
	USAGE:
	-h,--help               Print this message
        -e,--environment        Create a conda environment and install all the required software [Requires Conda]
	-r,--run                Run the pipeline [Make sure all the required software is installed]
	-c,--clear              Clear the space from residual files and folders from past runs


	SOFTWARE REQUIREMENTS*:
        - AdapterRemoval v2	(https://github.com/mikkelschubert/adapterremoval)
        - FastQC		(https://github.com/s-andrews/FastQC)
        - BWA			(https://github.com/lh3/bwa)
        - samtools		(http://www.htslib.org/download/)
        - Dedup			(https://github.com/apeltzer/DeDup)
        - MapDamage		(https://github.com/ginolhac/mapDamage)

        * By running bash $0 --environment you can automatically install all the required software via Conda
	"
}

set_env() {
	conda create -n plant_aDNA_pipeline -c bioconda AdapterRemoval FastQC bwa samtools dedup mapdamage2
	conda install -n plant_aDNA_pipeline -c conda-forge r-rcppgsl r-rcpp r-inline r-ggplot2 r-gam # MapDamage R requirements
}

requirements() {
	# CHECKING SOFTWARE REQUIREMENTS
	echo "Checking software requirements.
	
	Press ENTER if the following paths are correct. Otherwise provide the full path for the program:
	"
	read -p "AdapterRemoval v2 binary ($(which AdapterRemoval)): " AdapterRemoval_bin
	if [ -z $AdapterRemoval_bin ]; then AdapterRemoval_bin=AdapterRemoval; fi
	rc=$($AdapterRemoval_bin 2>/dev/null || echo -e "$?")
	if [ $rc -eq 127 ]; then echo "[ERROR] AdapterRemoval v2 is not available (https://github.com/mikkelschubert/adapterremoval)"; exit ;fi
	
	read -p "FastQC binary ($(which fastqc)): " fastqc_bin
	if [ -z $fastqc_bin ]; then fastqc_bin=fastqc; fi
	$fastqc_bin -h >/dev/null 2>/dev/null
	rc=$(echo $?)
	if [ $rc -eq 127 ]; then echo "[ERROR] FastQC is not available (https://github.com/s-andrews/FastQC)"; exit ; fi
	
	read -p "bwa binary ($(which bwa)): " bwa_bin
	if [ -z $bwa_bin ]; then bwa_bin=bwa; fi
	rc=$($bwa_bin 2>/dev/null || echo -e "$?")
	if [ $rc -eq 127 ]; then echo "[ERROR] BWA is not available (https://github.com/lh3/bwa)"; exit; fi
	
	read -p "samtools binary ($(which samtools)): " samtools_bin
	if [ -z $samtools_bin ]; then samtools_bin=samtools; fi
	rc=$($samtools_bin 2>/dev/null || echo -e "$?")
	if [ $rc -eq 127 ]; then echo "[ERROR] samtools is not available (http://www.htslib.org/download/)"; exit; fi
	
	read -p "dedup conda_wrapper/.jar ($(which dedup)): " dedup_bin
	if [ -z $dedup_bin ]; then dedup_bin=dedup; fi
	$dedup_bin -v >/dev/null 2>/dev/null
	rc=$(echo $?)
	if [ $rc -eq 126 ]; then dedup_bin="java -jar $dedup_bin" ; fi
	$dedup_bin -v >/dev/null 2>/dev/null
	rc=$(echo $?)
	if [ $rc -eq 127 ]; then echo "[ERROR] DeDup is not available (https://github.com/apeltzer/DeDup)"; exit ; fi
	
	read -p "mapDamage binary ($(which mapDamage)): " mapDamage_bin
	if [ -z $mapDamage_bin ]; then mapDamage_bin=mapDamage; fi
	rc=$($mapDamage_bin 2>/dev/null || echo -e "$?")
	if [ $rc -eq 127 ]; then echo "[ERROR] mapDamage is not available (https://github.com/ginolhac/mapDamage)"; exit ; fi
	
	header
	echo ""
	read -p "How many cores do you want to use? (Default 4): " nthreads
	if [ -z $nthreads ]; then nthreads=4; fi
	echo ""
}

openplot() { # Tries open (Mac) or xdg-open (Linux)
	open $1 2> /dev/null
	rc=$(echo $?)
	if [ $rc -eq 127 ]; then
		xdg-open $1
	fi
}

pipeline() {
	RED="\033[0;31m"
	NC="\033[0m"
	# CREATE FOLDERS
	header
	echo -e "
	################################################################################
	################################ ${RED}CREATE FOLDERS ${NC}################################
	################################################################################
	
	The following folders will be created:
	- 1_initial_data/
	- 1_initial_data/reference_genome/
	- 2_trimmed_merged/
	- 3_quality_control/
	- 4_mapping/
	- 5_aDNA_characteristics/
	
	COMMAND:
	mkdir -p 1_initial_data/reference_genome 2_trimmed_merged 3_quality_control 4_mapping 5_aDNA_characteristics
	
	(...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	mkdir -p 1_initial_data/reference_genome 2_trimmed_merged 3_quality_control 4_mapping 5_aDNA_characteristics &&
	
	# DOWNLOAD RAW DATA
	header
	echo -e "
	################################################################################
	############################## ${RED}DOWNLOAD RAW DATA  ${NC}##############################
	################################################################################
	
	We will use MiSeq pair-end Illumina sequences from a herbarium Arabidopsis thaliana sample (Weiss et al., 2016 ; https://doi.org/10.1098/rsos.160239)
	
	Download the PAIR 1 of reads to the folder 1_initial_data/
	
	COMMAND:
	wget -P 1_initial_data ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR964/ERR964430/s_1_s8_R1.fastq.gz
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	wget -P 1_initial_data ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR964/ERR964430/s_1_s8_R1.fastq.gz &&
	header
	echo -e "
	
	Download the PAIR 2 of reads to the folder 1_initial_data/
	
	COMMAND:
	wget -P 1_initial_data ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR964/ERR964430/s_1_s8_R2.fastq.gz
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	wget -P 1_initial_data ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR964/ERR964430/s_1_s8_R2.fastq.gz &&
	
	header
	echo -e "
	The raw reads have been downloaded and stored in 1_initial_data/
	"
	ls -lh 1_initial_data/
	echo -e "
	(...Press ${RED}ENTER ${NC}to continue...)"
	read
	
	# TRIMMING AND MERGING
	header
	echo -e "
	################################################################################
	############################# ${RED}TRIMMING AND MERGING ${NC}#############################
	################################################################################
	
	AdapterRemoval will be used to both remove Illumina adapters and merge paired reads when possible.
	
	OPTIONS:
	--file1 1_initial_data/s_1_s8_R1.fastq.gz : Location of PAIR 1 file
	--file2 1_initial_data/s_1_s8_R2.fastq.gz : Location of PAIR 2 file
	--qualitybase 64 : Input Quality Score encoding is 64
	--qualitybase-output 33 : Output Quality Score encoding is 33 (standard)
	--basename 3_trimmed_merged/historicAthaliana : Location and prefix for the output files
	--collapse : Merge trimmed reads when possible
	--gzip : Output is compressed
	--threads $nthreads : Number of threads to be used
	
	COMMAND:
	$AdapterRemoval_bin --file1 1_initial_data/s_1_s8_R1.fastq.gz --file2 1_initial_data/s_1_s8_R2.fastq.gz --qualitybase 64 --qualitybase-output 33 --basename 2_trimmed_merged/historicAthaliana --collapse --gzip --threads $nthreads
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$AdapterRemoval_bin --file1 1_initial_data/s_1_s8_R1.fastq.gz --file2 1_initial_data/s_1_s8_R2.fastq.gz --qualitybase 64 --qualitybase-output 33 --basename 2_trimmed_merged/historicAthaliana --collapse --gzip --threads $nthreads &&
	
	header
	echo "
	The reads have been trimmed and merged.
	"
	ls -lh 2_trimmed_merged/
	echo -e "
	(...Press ${RED}ENTER ${NC}to continue...)"
	read
	
	# QUALITY CONTROL
	header
	echo -e "
	################################################################################
	############################### ${RED}QUALITY CONTROL  ${NC}###############################
	################################################################################
	
	To assess the read quality of the sample, we will use the program FastQC.
	
	OPTIONS:
	-o 3_quality_control/ : The results will be stored in 2_QC/
	-t $nthreads : Number of threads to be used
	2_trimmed_merged/historicAthaliana.collapsed.gz : Location of the trimmed and merged reads
	2_trimmed_merged/historicAthaliana.pair1.truncated.gz : Location of the trimmed PAIR 1 reads
	2_trimmed_merged/historicAthaliana.pair2.truncated.gz : Location of the trimmed PAIR 2 reads
	
	COMMAND:
	$fastqc_bin -o 3_quality_control/ -t $nthreads 2_trimmed_merged/historicAthaliana.collapsed.gz 2_trimmed_merged/historicAthaliana.pair1.truncated.gz 2_trimmed_merged/historicAthaliana.pair2.truncated.gz

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$fastqc_bin -t $nthreads -o 3_quality_control/ 2_trimmed_merged/historicAthaliana.collapsed.gz 2_trimmed_merged/historicAthaliana.pair1.truncated.gz 2_trimmed_merged/historicAthaliana.pair2.truncated.gz &&
	
	echo -e "
	
	Quality control finished.
	Inspect the quality control files
	(...Press ${RED}ENTER ${NC}to continue...)"
	read
	openplot 3_quality_control/historicAthaliana.collapsed_fastqc.html 2>/dev/null &
	openplot 3_quality_control/historicAthaliana.pair1.truncated_fastqc.html 2>/dev/null &
	openplot 3_quality_control/historicAthaliana.pair2.truncated_fastqc.html 2>/dev/null &

	# Mapping
	header
	echo -e "
	################################################################################
	################################### ${RED}MAPPING  ${NC}###################################
	################################################################################
	
	Since the reads originate from an A. thaliana herbarium sample, the most suitable reference genome to use is the TAIR-10 A. thaliana assembly. The following command will download it and place it inside 1_initial_data/reference_genome/
	
	COMMAND:
	curl ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz | zcat > 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	curl ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz | zcat > 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa &&
	
	header
	echo -e "
	The reference genome must be indexed both by BWA and samtools.
	Command to generate index with bwa.
	
	COMMAND:
	$bwa_bin index 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$bwa_bin index 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa &&
	
	header
	echo -e "
	Command to generate index with samtools
	
	COMMAND:
	$samtools_bin faidx 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$samtools_bin faidx 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa &&
	
	echo "
	
	Indexing finished"
	ls -lh 1_initial_data/reference_genome/
	echo -e "
	(...Press ${RED}ENTER ${NC}to continue...)
	"
	read
	
	header
	echo -e "
	BWA aln will be used to map the preprocessed (trimmed and merged) reads to the A. thaliana reference genome.
	
	OPTIONS:
	-l 1024 : Set the seed to a very high value (1024) to inactivate it
	-f mapping/historicAthaliana.collapsed.sai : Path and name for the output file
	-t $nthreads : Number of threads to be used
	1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa : Location of the indexed reference genome
	2_trimmed_merged/historicAthaliana.collapsed.gz : Location of the trimmed and merged reads
	-t $nthreads : Number of threads to be used
	
	COMMAND:
	$bwa_bin aln -l 1024 -f 4_mapping/historicAthaliana.collapsed.sai -t $nthreads 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 2_trimmed_merged/historicAthaliana.collapsed.gz

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$bwa_bin aln -l 1024 -f 4_mapping/historicAthaliana.collapsed.sai -t $nthreads 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 2_trimmed_merged/historicAthaliana.collapsed.gz &&
	
	echo -e "
	Mapping finished
	(...Press ${RED}ENTER ${NC}to continue...)"
	read
	
	header
	echo '
	BWA samse will be used to convert the mapped reads (.sai file) into an alignment in the standard SAM format.
	
	OPTIONS:
	-r @RG\\tID:sample1\\tSM:sample1 : Read Group tag as "sample1". This will be important for downstream analyses
	-f 4_mapping/historicAthaliana.collapsed.sam : Location for the aligned sam file (output)
	1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa : Location of the indexed reference genome
	4_mapping/historicAthaliana.collapsed.sai : Location of the .sai file
	2_trimmed_merged/historicAthaliana.collapsed.gz : Location of the trimmed and merged reads
	
	COMMAND:
	$bwa_bin samse -r @RG\\tID:sample1\\tSM:sample1 -f 4_mapping/historicAthaliana.collapsed.sam 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 4_mapping/historicAthaliana.collapsed.sai 2_trimmed_merged/historicAthaliana.collapsed.gz'
	echo -e "
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$bwa_bin samse -r @RG\\tID:sample1\\tSM:sample1 -f 4_mapping/historicAthaliana.collapsed.sam 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 4_mapping/historicAthaliana.collapsed.sai 2_trimmed_merged/historicAthaliana.collapsed.gz &&
	echo "
	SAM-formatted alignment generated.
	"
	ls -lh 4_mapping/
	echo -e "
	(...Press ${RED}ENTER ${NC}to continue...)
	"
	read
	
	header
	echo -e "
	samtools flagstat will be used to estimate the endogenous DNA as the proportion of mapped reads.
	
	OPTIONS:
	-@ $nthreads : Number of threads to be used
	4_mapping/historicAthaliana.collapsed.sam : Location of the alignment (input) 
	> 4_mapping/historicAthaliana.collapsed.stats : Location for the output
	
	COMMAND:
	$samtools_bin flagstat -@ $nthreads 4_mapping/historicAthaliana.collapsed.sam > 4_mapping/historicAthaliana.collapsed.stats
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$samtools_bin flagstat -@ $nthreads 4_mapping/historicAthaliana.collapsed.sam > 4_mapping/historicAthaliana.collapsed.stats
	
	echo -e "
	Inspect the content of the produced file and assess the endogenous DNA as the proportion of mapped reads.
	
	COMMAND:
	cat 4_mapping/historicAthaliana.collapsed.stats
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	cat 4_mapping/historicAthaliana.collapsed.stats
	
	echo -e "
	(...Press ${RED}ENTER ${NC}to continue...)"
	read
	
	header
	echo -e "
	Since the SAM (Sequence Alignment Map) format is a plain text file, it can be efficiently compressed into a BAM (Binary Alignment Map) format using samtools. This allows to save disk space.
	Moreover, unmapped reads can be discarded.
	
	OPTIONS:
	-@ $nthreads : Number of threads to be used
	-Sb : The input is (S)AM and the output is (b)am
	-h : Include the header in the output
	-F 4 : Discard unmapped reads
	-o 4_mapping/historicAthaliana.mapped.bam : Location for the output file
	4_mapping/historicAthaliana.collapsed.sam : Location of the input file
	
	COMMAND:
	$samtools_bin view -@ $nthreads -Sb -h -F 4 -o 4_mapping/historicAthaliana.mapped.bam 4_mapping/historicAthaliana.collapsed.sam
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$samtools_bin view -@ $nthreads -Sb -h -F 4 -o 4_mapping/historicAthaliana.mapped.bam 4_mapping/historicAthaliana.collapsed.sam
	
	header
	echo -e "
	At this stage, it is possible to remove the uncompressed .sam file since the mapped reads are already compressed in the .bam file
	
	COMMAND:
	rm 4_mapping/historicAthaliana.collapsed.sam
	
	(Press ${RED}ENTER ${NC}to execute the command...)"
	read
	rm 4_mapping/historicAthaliana.collapsed.sam
	
	header
	echo -e "
	Most downstream analyses require coordinate-based sorted files.
	samtools can be used to sort the mapped bam file.
	
	OPTIONS:
	-@ $nthreads :  Number of threads to be used
	-o 4_mapping/historicAthaliana.mapped.sorted.bam : Location for the sorted output
	4_mapping/historicAthaliana.mapped.bam : Location of the unsorted mapped reads (input)
	
	COMMAND:
	$samtools_bin sort -@ $nthreads -o 4_mapping/historicAthaliana.mapped.sorted.bam 4_mapping/historicAthaliana.mapped.bam
	
        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$samtools_bin sort -@ $nthreads -o 4_mapping/historicAthaliana.mapped.sorted.bam 4_mapping/historicAthaliana.mapped.bam
	
	echo -e "
	
	The unsorted .bam file can be removed
	
	COMMAND:
	rm 4_mapping/historicAthaliana.mapped.bam

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	rm 4_mapping/historicAthaliana.mapped.bam
	
	header
	echo -e "
	################################################################################
	############################## ${RED}REMOVE DUPLICATES  ${NC}##############################
	################################################################################
	
	For downstream analyses it is important to remove the optical PCR duplicates.
	
	OPTIONS:
	-i 4_mapping/historicAthaliana.mapped.sorted.bam : Location of the mapped and sorted file (input)
	-m : The input contains merged reads
	-o 4_mapping/ : Location of the output folder
	
	COMMAND:
	$dedup_bin -i 4_mapping/historicAthaliana.mapped.sorted.bam -m -o 4_mapping/

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$dedup_bin -i 4_mapping/historicAthaliana.mapped.sorted.bam -m -o 4_mapping/
	
	echo -e "
	
	Duplicates removed.
	(...Press ${RED}ENTER ${NC}to continue)"
	read
	
	
	header
	echo -e "
	################################################################################
	####################### ${RED}ASSESSING ANCIENT DNA PROPERTIES ${NC}#######################
	################################################################################
	
	
	OPTIONS:
	-i 4_mapping/historicAthaliana.mapped.sorted_rmdup.bam : Location of the input file
	-r 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa : Location of the reference genome
	--no-stats : Disable the statistical calculation
	-y 0.05 : Plot Y-axis for nucleotide misincorporation frequencies plot up to 0.05
	-d 5_aDNA_characteristics/ : Location of the output folder
	
	COMMAND:
	$mapDamage_bin -i 4_mapping/historicAthaliana.mapped.sorted_rmdup.bam -r 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --no-stats -y 0.05 -d 5_aDNA_characteristics

        (...Press ${RED}ENTER ${NC}to execute the command...)"
	read
	$mapDamage_bin -i 4_mapping/historicAthaliana.mapped.sorted_rmdup.bam -r 1_initial_data/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --no-stats -y 0.05 -d 5_aDNA_characteristics
	
	echo -e "
	
	Inspect the generated files
	
	(...Press ${RED}ENTER ${NC}to continue)"
	read
	openplot 5_aDNA_characteristics/Fragmisincorporation_plot.pdf 2>/dev/null &
	
	echo "
	Finished...

	I you used the -e/--environment option, you can remove the environment by typing:
	conda deactivate
	conda env remove -n plant_aDNA_pipeline
	"
}


main() {
        RED="\033[0;31m"
	NC="\033[0m"

	if [ -z $1 ] || [ $1 == "-h" ] || [ $1 == "--help" ]; then
		help_message
		exit
	elif [ $1 == "-c" ] || [ $1 == "--clear" ]; then
		rm -r 1_initial_data 2_trimmed_merged 3_quality_control 4_mapping 5_aDNA_characteristics 2>/dev/null
		echo -e "
		The space has been cleared
		You can proceed by typing bash $0 --run
		"
		exit
	elif [ $1 == "-r" ] || [ $1 == "--run" ]; then
		requirements
		pipeline
	elif [ $1 == "-e" ] || [ $1 == "--environment" ]; then
		set_env &&
		echo -e "
		##########################################################################################
		####										      ####
		#### ${RED}Activate the conda environment and run the pipeline with the following commands:${NC} ####

		     ${RED}conda activate plant_aDNA_pipeline
		     bash $0 --run${NC}								      ####
		####										      ####
		##########################################################################################
		"
		exit
	else
		help_message
		exit
	fi
}

header
main $1
