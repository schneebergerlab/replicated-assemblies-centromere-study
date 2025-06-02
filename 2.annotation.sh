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

