#!/bin/bash

################################################################################
# Main Pipeline: Centromeric Mutation Simulation (CEN178)
# ------------------------------------------------------------------------------
# This script simulates mutations in a tandem repeat array of CEN178 sequences
# over 6 million generations, with downstream analyses:
# - Mutation accumulation
# - Unit-level extraction and consensus generation
# - Statistical summaries and visualization
# - Heatmap generation and video synthesis
#
# Output Summary:
# - 1.fasta/: Full-length centromere arrays per generation
# - 2.records/: Mutation records per generation step
# - 3.generation_unit_pos/: Unit boundary positions per generation
# - 4.generation_unit_fasta/: Extracted unit FASTA files and BED files
# - 5.generation_consensus_fasta/: Aligned unique units and consensus sequences
# - 6.moddotplot/: ModDotPlot BED files, heatmap PNGs and video
# - all.generations.*: Aggregated stats (AT content, edit distance, counts)
# - 5tests.all.generations.*: Merged stats across 5 simulations
#
# External Dependencies: Python, seqtk, seqkit, mafft, usearch, ModDotPlot, ffmpeg
# Author: xdong@mpipz.mpg.de
################################################################################

##############################
# Step 1: Initial Simulation #
##############################
python ./bin/simulate.py \
    ./data/15000copy_cen178.seq \
    1000 \
    ./data/15000copy_cen178.178bp.bed.pos \
    1.fasta/1000generation.out.fa \
    2.records/1000generation.record.txt \
    3.generation_unit_pos/1000generation.unit.pos

#############################################################
# Step 2: Continue Simulation (from 2000 to 6,000,000 gens) #
#############################################################
for i in $(seq 2000 1000 6000000); do
    prev_gen=$((i - 1000))
    python ./bin/simulate.py \
        1.fasta/${prev_gen}generation.out.fa \
        1000 \
	3.generation_unit_pos/${prev_gen}generation.unit.pos
        1.fasta/${i}generation.out.fa \
        2.records/${i}generation.record.txt \
        3.generation_unit_pos/${i}generation.unit.pos
done

########################################
# Step 3: Add Headers to FASTA Files   #
########################################
for test in $(seq 1 5); do
    for gen in $(seq 1000 1000 6000000); do
        echo ">centro_${gen}gen" | cat - test${test}/1.fasta/${gen}generation.out.fa \
            > test${test}/1.fasta/${gen}generation.out.fa.tmp
        mv test${test}/1.fasta/${gen}generation.out.fa.tmp test${test}/1.fasta/${gen}generation.out.fa
    done
done

########################################
# Step 4: Extract Individual Units     #
########################################
for gen in $(seq 1000 1000 6000000); do
    sed '$d' 3.generation_unit_pos/${gen}generation.unit.pos > 4.generation_unit_fasta/${gen}generation.unit.pos1
    sed '1d' 3.generation_unit_pos/${gen}generation.unit.pos > 4.generation_unit_fasta/${gen}generation.unit.pos2

    paste 4.generation_unit_fasta/${gen}generation.unit.pos1 4.generation_unit_fasta/${gen}generation.unit.pos2 \
        | awk -v gen=$gen '{print "centro_"gen"gen\t"$2"\t"$4}' \
        > 4.generation_unit_fasta/${gen}generation.unit.bed

    seqtk subseq 1.fasta/${gen}generation.out.fa 4.generation_unit_fasta/${gen}generation.unit.bed \
        > 4.generation_unit_fasta/${gen}generation.unit.fa

    rm 4.generation_unit_fasta/${gen}generation.unit.pos1 4.generation_unit_fasta/${gen}generation.unit.pos2
done

############################################
# Step 5: Generate Consensus Sequences     #
#         + Collect Stats (GC, edits, etc) #
############################################
for gen in $(seq 1000 1000 6000000); do
    usearch -fastx_uniques 4.generation_unit_fasta/${gen}generation.unit.fa \
        -fastaout 4.generation_unit_fasta/${gen}generation.unit.uniq.fa

    mafft --thread 10 --retree 2 \
        4.generation_unit_fasta/${gen}generation.unit.uniq.fa \
        > 5.generation_consensus_fasta/${gen}generation.unit.uniq.fa.mafft.aln

    python ./bin/repeat_uniq.consensus.stat.py \
        5.generation_consensus_fasta/${gen}generation.unit.uniq.fa.mafft.aln \
        5.generation_consensus_fasta/${gen}generation.unit.uniq.fa.mafft.stat

    python ./bin/get_repeat_uniq_consensus_seq.py \
        5.generation_consensus_fasta/${gen}generation.unit.uniq.fa.mafft.stat \
        5.generation_consensus_fasta/${gen}generation.unit.uniq.consensus.fa

    awk -v gen=$gen '{print gen"\t"$0}' \
        5.generation_consensus_fasta/${gen}generation.unit.uniq.consensus.fa

    python ./bin/fasta_stat.py 4.generation_unit_fasta/${gen}generation.unit.fa \
        >> all.generations.fasta.stat

    python ./bin/count_cen178.py 4.generation_unit_fasta/${gen}generation.unit.fa ./data/cen178.fa \
        >> all.generations.cen178.count

    python ./bin/edit_distance_Levenshtein.py \
        4.generation_unit_fasta/${gen}generation.unit.fa \
        ./data/cen178.fa \
        4.generation_unit_fasta/${gen}generation.unit.fa.edit_distance \
        >> all.generations.avg.cen178.edit_distance

    seqkit fx2tab --gc 1.fasta/${gen}generation.out.fa \
        | awk '{print $1"\t"100-$3}' \
        >> all.generations.fasta.at_content
done

##############################################
# Step 6: Cleanup File Paths in Stats Tables #
##############################################
sed -i 's/^4.generation_unit_fasta\///g' all.generations.*
sed -i 's/generation.unit.fa//g' all.generations.*

##########################################
# Step 7: Merge All 5 Test Statistics    #
##########################################
for test in $(seq 1 5); do
    for file in fasta.stat cen178.count avg.cen178.edit_distance unit.uniq.consensus.fa fasta.at_content; do
        awk -v test=$test '$1 <= 6000000 {print $0"\tTest"test}' \
            test${test}/all.generations.${file} \
            >> 5tests.all.generations.${file}
    done
done

########################################
# Step 8: ModDotPlot Heatmaps         #
########################################
source ~/software/ModDotPlot/venv/bin/activate
maxseq=4600000

# Initial (15K) reference heatmap
moddotplot static -f data/15000copy_cen178.seq -w 10000 -id 0 --no-plot -o 6.moddotplot/1.bed_iden0/
python ./bin/heatmap_moddotplot.py \
    6.moddotplot/1.bed_iden0/15000copy_cen178.bed \
    6.moddotplot/3.png/centro_15K_CEN178.heatmap.png \
    --palette-colors '#313695' '#4575B4' '#74ADD1' '#ABD9E9' '#E0F3F8' '#FEE090' '#FDAE61' '#F46D43' '#D73027' '#A50026' \
    --breakpoints 70 73 76 79 82 85 88 91 94 97 100 \
    --width 7 --height 7 -xlab "Array Size (bp)" -ylab "Array Size (bp)" \
    --title 0gen --xlim 0 $maxseq --ylim 0 $maxseq

# Heatmaps for generations
for gen in $(seq 10000 10000 6000000); do
    moddotplot static -f 1.fasta/${gen}generation.out.fa -w 10000 -id 0 --no-plot -o 6.moddotplot/1.bed_iden0/

    awk '$7 >= 70' 6.moddotplot/1.bed_iden0/centro_${gen}gen.bed \
        > 6.moddotplot/2.bed_iden70/centro_${gen}gen.iden70.bed

    python ./bin/heatmap_moddotplot.py \
        6.moddotplot/2.bed_iden70/centro_${gen}gen.iden70.bed \
        6.moddotplot/3.png/centro_${gen}gen.iden70.heatmap.png \
        --palette-colors '#313695' '#4575B4' '#74ADD1' '#ABD9E9' '#E0F3F8' '#FEE090' '#FDAE61' '#F46D43' '#D73027' '#A50026' \
        --breakpoints 70 73 76 79 82 85 88 91 94 97 100 \
        --width 7 --height 7 -xlab "Array Size (bp)" -ylab "Array Size (bp)" \
        --title ${gen}gen --xlim 0 $maxseq --ylim 0 $maxseq
done

########################################
# Step 9: Generate Heatmap Video       #
########################################
ln -s 6.moddotplot/3.png/centro_15K_CEN178.heatmap.png 6.moddotplot/4.ffmpeg/0000.png

for i in $(seq 1 1200); do
    gen=$((i * 5000))
    printf -v padded "%04d" $i
    ln -s 6.moddotplot/3.png/centro_${gen}gen.iden70.heatmap.png 6.moddotplot/4.ffmpeg/${padded}.png
done

ffmpeg -framerate 30 -i 6.moddotplot/4.ffmpeg/%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p 6.moddotplot/4.ffmpeg/simulation.mp4

