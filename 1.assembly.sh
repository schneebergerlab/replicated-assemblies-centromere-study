#!/bin/bash

################################################################################
# Genome Assembly and Evaluation Pipeline
# ------------------------------------------------------------------------------
# This script runs de novo genome assembly using multiple tools (HiFiasm, IPA,
# Peregrine, Canu, Flye), performs reference-based scaffolding, and evaluates 
# assemblies using NG50, Merqury, and BUSCO. It also generates a consensus 
# reference genome for coordinate-based mutation comparison.
#
# Author: xdong@mpipz.mpg.de
################################################################################


################################################################################
# Step 1: De Novo Assembly with Multiple Tools
################################################################################

# Assemble with HiFiasm (best contiguity; used for final assemblies)
# -l0 disables duplication (for homozygous/inbred genomes)
hifiasm -l0 -t 20 -o A1 A1.hifi.fastq.gz

# Convert HiFiasm GFA output to FASTA
awk '/^S/{print ">"$2;print $3}' A1.hifiasm_l0.bp.p_ctg.gfa > A1.hifiasm_l0.bp.p_ctg.fasta


# Assemble with IPA (default parameters)
ipa local --nthreads 24 --njobs 1 -i A1.hifi.fastq.gz


# Assemble with Peregrine (default parameters)
pg_run.py asm \
    ./peregrine/hifidata.lst 24 24 24 24 24 24 24 24 24 \
    --with-consensus --shimmer-r 3 --best_n_ovlp 8 \
    --output ./


# Assemble with Canu
canu -p A1_canu -d A1_canu_asm genomeSize=135m -pacbio-hifi A1.hifi.fastq.gz \
    executiveThreads=24 &> A1_canu.log &


# Assemble with Flye
flye --pacbio-hifi A1.hifi.fastq.gz --threads 20 -o ./


# Assembly statistics (NG50 etc.) on primary contigs
bin/calc_CN50.pl A1.hifiasm_l0.bp.p_ctg.fasta 135000000 5 > \
    A1.hifiasm_l0.bp.p_ctg.fasta.N50.calc.result.txt



################################################################################
# Step 2: Reference-Based Scaffolding (Using Col-CEN)
################################################################################

# Detect organellar contamination
minimap2 -x asm5 TAIR10.organel.fa A1_hifiasm_l0.bp.p_ctg.fasta > A1.organel.out

# Merge overlapping regions
cut -f 1,3,4 A1.organel.out | sort -k1,1 -k2,2n | bedtools merge -d 1 > A1.organel.bed

# Create whole-contig BED file and calculate coverage
cut -f1,2 A1.organel.out | uniq | awk '{print $1"\t0\t"$2}' > A1.bed
bedtools coverage -a A1.bed -b A1.organel.bed > A1.organel.ratio.bed

# List contigs with >50% organellar content
awk '$7>=0.5{print $1}' A1.organel.ratio.bed > A1.organel.contig.list

# Reference-based scaffolding (non-organelle contigs only)
ragtag.py scaffold -o ragtag_output Col-CEN.fa A1_l0.non-organel.contig.fa -t 20



################################################################################
# Step 3: Assembly Evaluation with Merqury
################################################################################

conda activate merqury
cd $path2meryldb

# 1. Build meryl DBs for Illumina reads
for read in $(ls $path2shortdata/ | grep gz | awk -F "." '{print $1}'); do
    meryl k=18 count output ${read}.illumina.meryl ${read}.illumina.fq.gz
done

# 2. Merge meryl DBs
meryl union-sum output ${sample}.meryl A*.meryl

# 3. Run Merqury
merqury.sh ./meryl_dbs/${sample}.meryl ../A1.hifiasm_l0.bp.p_ctg.fasta output &



################################################################################
# Step 4: Assembly Evaluation with BUSCO
################################################################################

conda activate busco5
busco -i A1.hifiasm_l0.bp.p_ctg.fasta \
      -l ~/busco_downloads/embryophyta_odb10 \
      -o busco_embryo_odb10 -m genome -c 10 --offline



################################################################################
# Step 5: Generate Consensus Genome for Coordinate Comparison
################################################################################

# Convert variant table to VCF format
awk '{print $1"\t"$2"\t.\t"$3"\t"$4"\t.\tPASS"}' change.txt > change.vcf.txt
cat ./data/bcftools_consensus.head change.vcf.txt > change.vcf

# Compress and index VCF
bgzip change.vcf
bcftools index change.vcf.gz

# Apply VCF changes to original genome to create consensus
conda activate bcftools-1.16
cat B1.ragtag.chr.fa | bcftools consensus -c chain-1.txt change.vcf.gz > F0.fa

