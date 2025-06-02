#!/bin/bash

################################################################################
# Assembly Comparison and Error Validation Pipeline
# ------------------------------------------------------------------------------
# This script compares genome assemblies using SYRI, categorizes structural 
# variants, validates candidate assembly errors using Illumina and HiFi reads,
# and evaluates correction tools such as Pilon, DeepVariant, pbsv, and Sniffles.
#
# Author: xdong@mpip.mpg.de
################################################################################


################################################################################
# Step 1: Identify Assembly Differences Using SYRI
################################################################################

ref=A2
qry=A1

# Align assemblies
minimap2 -ax asm5 -t 20 --eqx ./${ref}.chr.fasta ./${qry}.chr.fasta | samtools sort -o ${qry}-${ref}.sorted.bam

# Run SYRI
conda activate syri
syri -c ${qry}-${ref}.sorted.bam -r ./${ref}.chr.fasta -q ./${qry}.chr.fasta -F B --prefix ${qry}-${ref}.

# Filter structural variants from SYRI output
cat ${qry}-${ref}.syri.out | awk '
$11 !~ /SYN|DUPAL|INVAL|INVDPAL|TRANSAL/ {
    if ($11 == "SNP" && ($4 == "N" || $5 == "N"))
        print $1, $2, $3, $4, $6, $7, $8, $5, $9, $11, "N";
    else if ($11 == "INS" && $4 == "N" && length($5) > 50)
        print $1, $2, $3, $4, $6, $7, $8, $5, $9, $11, "N;big";
    else if ($11 == "INS" && $4 != "N" && length($5) > 50)
        print $1, $2, $3, $4, $6, $7, $8, $5, $9, $11, "big";
    else if ($11 == "DEL" && $5 == "N" && length($4) > 50)
        print $1, $2, $3, $4, $6, $7, $8, $5, $9, $11, "N;big";
    else if ($11 == "DEL" && $5 != "N" && length($4) > 50)
        print $1, $2, $3, $4, $6, $7, $8, $5, $9, $11, "big";
    else
        print $1, $2, $3, $4, $6, $7, $8, $5, $9, $11, "-";
}' OFS='\t' > ${qry}-${ref}.error.list

# Generate summary stats
cut -f10,11 ${qry}-${ref}.error.list | sort | uniq -c > ${qry}-${ref}.error.stat.txt


################################################################################
# Step 2: Quantify Types of Structural Differences
################################################################################

cat ref_A1/A2-A1.syri.out | awk '$11=="INS" && length($5)<=50 {sum += length($5)-1} END {print "DEL_small:", sum}'
cat ref_A1/A2-A1.syri.out | awk '$11=="INS" && length($5)>50 {sum += length($5)-1} END {print "DEL_big:", sum}'
cat ref_A1/A2-A1.syri.out | awk '$11=="DEL" && length($4)<=50 {sum += length($4)-1} END {print "INS_small:", sum}'
cat ref_A1/A2-A1.syri.out | awk '$11=="DEL" && length($4)>50 {sum += length($4)-1} END {print "INS_big:", sum}'
cat ref_A1/A2-A1.syri.out | awk '$11=="DUP" {sum += length($4)-1} END {print "DUP:", sum}'


################################################################################
# Step 3: Validate Assembly Errors via Read Alignments
################################################################################

# Align HiFi reads to assembly
minimap2 -ax map-hifi -t 20 -R '@RG\tID:B1\tSM:B1' A1.all.fasta \
  A1.hifi.reads.fa.gz | \
  samtools sort -@ 20 -o A1.hifi.sorted.bam
samtools index A1.hifi.sorted.bam

# Align Illumina reads to assembly
bwa index A1.chr.fasta
bwa mem -t 20 -R '@RG\tID:A\tSM:A' A1.all.fasta \
  A.illumina_1.fq.gz \
  A.illumina_2.fq.gz | \
  samtools sort -@ 20 -o A.illu.A1.sorted.bam


################################################################################
# Step 4: Evaluate Assembly Correction Tools
################################################################################

# Step 4.1 Run Pilon 
java -Xmx200G -jar /path/to/pilon.jar \
   --genome A1.all.fasta \
   --frags A.illu.A1.sorted.bamm \
   --fix snps,indels \
   --vcf --changes --tracks --output pilon


# Step 4.2 Run DeepVariant (multiple modes: PACBIO, WGS and HYBRID)
for i in $(seq 1 5); do
  singularity run --bind ${PWD} /path/to/deepvariant.img \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref ragtag.scaffold.chr.fasta \
    --reads A1.hifi.A1.all.sorted.bam \
    --intermediate_results_dir ./tmp \
    --output_vcf Mock-A.hifi-illu.Chr${i}_RagTag.vcf.gz \
    --num_shards 4 \
    --regions Chr${i}_RagTag
done

# Filter low-quality variant calls
  bcftools view -f "PASS" -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=30' -Oz \
    ${i}-A.hifi-illu.vcf.gz > ${i}-A.hifi-illu.PASS.vaf0.5.gq30.vcf.gz


################################################################################
# Step 5: Structural Variant Calling with pbsv and Sniffles
################################################################################
# Step 4.3 Run pbsv
conda activate pbsv
for i in A1 A2 B1 B2; do
  pbmm2 align -j 10 ./${i}.all.fasta \
    ./${i}.hifi.fq.gz ${i}/${i}.pbmm2.sorted.bam \
    --sort --preset CCS --sample ${i} --rg "@RG\tID:${i}"

  pbsv discover --hifi ${i}/${i}.pbmm2.sorted.bam ${i}/${i}.pbmm2.svsig.gz
  tabix -c '#' -s 3 -b 4 -e 4 ${i}/${i}.pbmm2.svsig.gz
  pbsv call -j 10 ../../../3.compare/check_hifi_alignment/${i}_ragtag/${i}.ragtag.scaffold.all.fasta \
    ${i}/${i}.pbmm2.svsig.gz ${i}/${i}.vcf

# Sniffles SV calling
for i in A1 A2 B1 B2; do
  sniffles --input ../../../3.compare/check_hifi_alignment/${i}_ragtag/all/${i}.hifi.${i}.all.sorted.bam \
    --vcf ${i}/${i}.vcf \
    --reference ../../2.ragtag_scaf_Col-CEN/${i}/ragtag.scaffold.chr.fasta -t 10

  bcftools view -f "PASS" -e 'FORMAT/GQ<=30' -Oz ${i}/${i}.vcf > ${i}/${i}.PASS.gq30.vcf.gz
done

