#!/bin/bash

################################################################################
# Repeat Annotation Pipeline
# ------------------------------------------------------------------------------
# This script performs repeat annotations on the reference genome (F0.fa):
# - CEN178 annotation via TRASH
# - rDNA and telomeric repeat annotation via RepeatMasker with a custom library
# - Simple sequence repeat (SSR) annotation using Arabidopsis-specific library
#
# Dependencies: TRASH, RepeatMasker
################################################################################

# Step 1: Centromere Annotation Using TRASH
# Annotate centromeric tandem repeats (e.g., CEN178) using TRASH with custom HOR class
~/software/TRASH/TRASH_run.sh F0.fa \
    --seqt ./data/CEN178.csv \
    --horclass CEN178 \
    --par 5 \
    --horonly \
    --o


# Step 2: rDNA and Telomeric Repeat Annotation with RepeatMasker
# Use a custom repeat library to annotate 5S rDNA, 45S rDNA, and telomere arrays
RepeatMasker -lib ./data/rDNA_NaishCEN_telomeres.fa \
    -nolow -gff -xsmall \
    -cutoff 200 \
    F0.fa


# Step 3: Simple Sequence Repeat (SSR) Annotation with RepeatMasker
# Perform SSR annotation using the built-in Arabidopsis repeat library
RepeatMasker -species "arabidopsis thaliana" \
    -s \
    -no_is \
    -cutoff 255 \
    -frag 20000 \
    F0.fa

# Step 4: Centromere Annotation Using CentroAnno
~/centroAnno/centroAnno ./F0.fa -o ./ -x anno-asm -m ./cen178.fa -t 20

# Step 5: TE annotation Using EDTA
#Annotation of TEs of TAIR10 was downloaded from https://www.arabidopsis.org/download/list?dir=Genes%2FAraport11_genome_release; while annotation of genes were lifted off from TAIR12 annotation under the BioProject accession PRJEB100887.
perl ~/EDTA/bin/EDTA.pl --genome F0.fa --cds ~/liftoff/F0.TAIR12_liftoff.cds.fa --curatedlib ./Araport11_GFF3.transposable_element_gene.fa --force 1

# Step 6: Annotation of ATHILA elements Using ATHILAfiner
bash ~/ATHILAfinder/ATHILAfinder.sh ./F0.fa ATHILA_PBSjunction.txt ATHILA_PPTjunction.txt 11000 2000 11000 2000 20 5 1000 2000 Atha Athila 0.80 500 2500 5 5 19 0.90 orfis.Ty3.updated.hmmdb 0.01
