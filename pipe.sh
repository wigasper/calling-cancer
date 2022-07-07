#!/bin/bash

set -Eeuo pipefail
trap cleanup SIGINT SIGTERM ERR

cleanup() {
	trap - SIGINT SIGTERM ERR
	echo "Received signal at $step"
	echo "UID $uid"
}

samples_tsv=samples.tsv

# Path to an unzipped GTF file
genome_gtf=

# Number of thraeds for STAR to use during indexing
# and alignment
star_n_threads=

# Number of threads for trim galore to use, diminishing
# returns after 4
trim_galore_n_threads=4

# Number of threads to use for variant calling, 
# suggested number_available / 3, since everything is piped
# and multiple process will be created and SnpEff will
# create a ton of threads. i.e. if you have 16 threads maybe
# set this value to 6
varcall_n_threads=

# Number of threads to use for coverage checking, 
# suggested number_available / 2, since there will be
# 2 processes created for each file.
covmachine_n_threads=

# Path to snpsubtract binary, if in $PATH just put 'snpsubtract'
snpsubtract_path= 

# Path to the common variants vcf, obtainable from NCBI dbSNP
common_variants_vcf=

# Path to an unzipped genome FASTA
export genome_fasta=

# String designating the SnpEff target genome (e.g., hg38), SnpEff will
# pull the DB for this if you don't have it. It needs to
# match your $genome_fasta and $genome_gtf
export snpeff_genome_target=

# Path to SnpEff jar
export snpeff_path=

# Make indices if they don't exist
if [[ ! -d indices ]]; then
	mkdir indices
	STAR --runThreadN $star_n_threads \
		--runMode genomeGenerate \
		--genomeDir ./indices \
		--genomeFastaFiles $genome_fasta \
		--sjdbGTFfile $genome_gtf
fi

# Make required dirs if they don't exist
mkdir -p fastqs
mkdir -p temp
mkdir -p final_variants
mkdir -p trimmed_fastqs
mkdir -p temp_snpeff_out
mkdir -p coverage_data

# Variant calling function, I do it like this so that it
# can easily be parallelized with parallel. Takes a uid
# as the only argument
varcall() {
	local uid=$1

	if [[ ! -f final_variants/$uid.final.vcf ]]; then
		step="varcall:$uid"
		echo "Calling variants on $uid"
		
		bcftools mpileup -f $genome_fasta temp/$uid.sorted.bam | \
			bcftools call --ploidy 2 -mv -Ob | \
			bcftools view -i '%QUAL>=30' -Ov -o temp/$uid.vcf

		vcfutils.pl varFilter temp/$uid.vcf > final_variants/$uid.final.vcf
		
		rm temp/$uid.vcf
	fi
	
	if [[ ! -f temp_snpeff_out/$uid.ann.vcf ]]; then
		step="snpeff"
		# do snpeff
		java -Xmx8g -jar $snpeff_path ann \
			-noStats \
			$snpeff_genome_target \
			final_variants/$uid.final.vcf > temp_snpeff_out/$uid.ann.vcf
		
	fi 
	step="Pipeline complete"
}

export -f varcall

# Function to assess coverage at driver residue bases
covmachine() {
	local uid=$1

	if [[ ! -f coverage_data/$uid.cov ]]; then
		bcftools mpileup -O v -f $genome_fasta temp/$uid.sorted.bam | \
			python3 covmachine.py analyze -o coverage_data/$uid.cov
	fi
}

export -f covmachine

# Main pipeline
#
# Reads SRA identifiers and UIDs from the $samples_tsv
# Retrieve and QC FASTQs if they don't exist
# Align and sort FASTQs
while read line; do
	step="pick_accession"

	srr=`echo "$line" | cut -f 1`
	uid=`echo "$line" | cut -f 2`

	r1=$uid\_R1.fastq
	r2=$uid\_R2.fastq

	if [[ ! -f trimmed_fastqs/$uid\_R1_val_1.fq ]]; then
		step="prefetch_fasterq:$srr:$uid"
		prefetch $srr
		fasterq-dump --split-files $srr/$srr.sra
		
		mv $srr\_1.fastq fastqs/$r1
		mv $srr\_2.fastq fastqs/$r2

		rm -rf $srr
		
		step="trim:$srr:$uid"
		echo "Trimming $r1 $r2"
		trim_galore --paired \
			--no_report_file \
			-o trimmed_fastqs \
			-j $trim_galore_n_threads \
			fastqs/$r1 fastqs/$r2
		
		rm fastqs/$r1
		rm fastqs/$r2
	fi

	if [[ ! -f temp/$uid.sorted.bam ]]; then		
		step="map:$srr:$uid"
		echo "Mapping"
		STAR --runThreadN $star_n_threads \
			--genomeDir ./indices \
			--outSAMtype BAM SortedByCoordinate \
			--readFilesIn trimmed_fastqs/$uid\_R1_val_1.fq trimmed_fastqs/$uid\_R2_val_2.fq
		mv Aligned.sortedByCoord.out.bam temp/$uid.sorted.bam
	fi

done < $samples_tsv

# Perform variant calling
cat $samples_tsv | cut -f 2 | parallel -j $varcall_n_threads varcall {} 

mkdir -p final_variants_subbed

$snpsubtract_path temp_snpeff_out $common_variants_vcf final_variants_subbed

# Annotate variants using OncoKB API
if [[ ! -f onco_annotations.json ]]; then
	python3 annotate_prot_variants.py final_variants_subbed onco_annotations.json
fi

# Identify the genomic positions of all identified driver alterations
if [[ ! -f driver_positions.json ]]; then
	python3 covmachine.py init -a onco_annotations.json -g $genome_gtf
fi

# Store coverage values for bases corresponding to all driver alterations
cat $samples_tsv | cut -f 2 | parallel -j $covmachine_n_threads covmachine {}
