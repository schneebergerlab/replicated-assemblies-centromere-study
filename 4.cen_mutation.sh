#!/bin/bash

################################################################################
# Centromere Mutation Analysis Pipeline
# ------------------------------------------------------------------------------
# This script is used to study the pattern of mutations in centromeric regions.
# Due to their highly repetitive nature, centromeres often yield misalignments
# that may lead to false positives/negatives in mutation calling.
# 
# To resolve this:
# - Word-based alignment is used to find the most parsimonious alignment
# - Left-alignment is performed on mutations
#
# To study the mutation pattern:
# - HOR scores and alignment to the consensus (CEN178) are used to study
#   mutational patterns across samples.
#
# Author: xdong@mpipz.mpg.de
################################################################################


################################################################################
# Step 1: Identify centromeric mutations based on TRASH annotations
# ------------------------------------------------------------------------------
# TRASH output offsets CEN178 units by 1 bp forward; correct for this.
################################################################################

cat CEN178.repeats.from.F0.fa.csv | grep Chr | \
awk -F "," '{print $8"\t"$1-2"\t"$2-1"\t"$5"\t"$4"\t"$3"\t"$9"\t"$10}' \
> CEN178.repeats.from.F0.fa.bed

bedtools intersect -wa -wb \
-a CEN178.repeats.from.F0.fa.bed \
-b mutation.bed \
> mutation.CEN178.repeats.from.F0.fa.bed


################################################################################
# Step 2: Refine alignments with word-based alignment
# ------------------------------------------------------------------------------
# Left-align indels, reduce ambiguity in clustered mutation regions.
################################################################################

python ./bin/word_based_alignment.rev_comp.py -l 150 -i 0 -o word_based.flank3K.output.txt ref.flank3k.fa qry.flank3k.fa

python ./bin/merge_segments_plot.py \
  -i word_based.flank3K.output.txt \
  -o word_based.flank3K.merged_segmetns.out \
  -pdf --highlight_reverse
# Optional: --label_internal_ends


################################################################################
# Step 3: Calculate mutation positions relative to the CEN178 consensus
# ------------------------------------------------------------------------------
# Relative positions help identify recurring mutation hotspots in repeats.
################################################################################

python ./bin/calculate.CEN178_relative_pos.py \
  ./F0.fa \
  mutation.pos \
  CEN178.repeats.from.F0.fa.bed \
  mutation.CEN178_relative_pos.txt


################################################################################
# Step 4: Compute HOR (Higher Order Repeat) score for centromeric mutations
# ------------------------------------------------------------------------------
# HOR score estimates how repetitive each region is.
################################################################################

# Extract repeat info and clean up
cat all.repeats.from.F0.fa.csv | grep CEN178 | \
awk -F "," '{print $8"\t"$1"\t"$2"\t"$10}' \
> all.repeats.from.F0.repetitiveness

sed -i 's/"//g' all.repeats.from.F0.repetitiveness

# Get maximum repetitiveness per chromosome
awk '{
  chr = $1;
  if (!max[chr] || $4 > max[chr]) {
    max[chr] = $4;
  }
} END {
  for (chr in max) {
    print chr, max[chr];
  }
}' all.repeats.from.F0.repetitiveness > max_per_chr.txt

# Normalize each line by its chromosome's max repetitiveness
awk 'BEGIN {
  while ((getline < "max_per_chr.txt") > 0) {
    max[$1] = $2;
  }
}
{
  norm = ($4 == 0 || max[$1] == 0) ? 0 : $4 / max[$1];
  print $0"\t"norm;
}' all.repeats.from.F0.repetitiveness > all.repeats.from.F0.repetitiveness.HORscore

# Map HOR score to mutation positions
bedtools map -a snp.bed -b CEN178.bed.repetitiveness.HORscore -c 5 -o mean > snp.HORscore.out

for i in dup del; do
  bedtools map -a ${i}.start.bed -b CEN178.bed.repetitiveness.HORscore -c 5 -o mean > ${i}.start.HORscore.out
done

# Map HOR score to 500 manually validated SNPs
cat 500SNP.list | grep -v Pos | \
awk -F ":" '{print $1"\t"$2}' | \
awk -F "-" '{print $1"\t"$2}' | \
awk '{print $1"\t"$2+$4-2"\t"$2+$4-1"\t"$5"\t"$6}' > 500SNP.bed

rm 500SNP.sorted.bed
for i in `cut -f1 500SNP.bed | sort | uniq`; do
  cat 500SNP.bed | grep $i | sort -nk2 >> 500SNP.sorted.bed
done

bedtools map -a 500SNP.sorted.bed -b CEN178.bed.repetitiveness.HORscore -c 5 -o mean > 500SNP.HORscore.out


################################################################################
# Step 5: Analyze in-frame mutation patterns via NUCmer alignment to CEN178
# ------------------------------------------------------------------------------
# Helps to visualize if insertions/deletions align with repeat structure.
################################################################################

for sample in A B
do
  nucmer --maxmatch -l 10 -c 20 $sample.cen.Indel.fa CEN178.fa -p $sample.cen.Indel.CEN178.nucmer
  show-coords -c $sample.cen.Indel.CEN178.nucmer.delta > $sample.cen.Indel.CEN178.nucmer.coords

  ./bin/mummerCoordsDotPlotly.R \
    -i $sample.cen.Indel.CEN178.nucmer.coords \
    -o $sample.cen.Indel.CEN178.nucmer \
    -m 10 -q 10 -k 5 -l -x
done
