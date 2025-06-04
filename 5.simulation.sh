#!/bin/bash

###############################################################################
# Non-allelic Gene Conversion and Mutaiton Simulations in CEN178 Arrays
# -----------------------------------------------------------------------------
# This pipeline performs two major simulations:
#
# 1. Non-Allelic Gene Conversion Simulation:
#    Simulates 500 random point mutations in centromeric repeat arrays and 
#    identifies potential gene conversion donors using k-mer-based analysis.
#
# 2. Centromeric Tandem Repeat Mutation Simulation:
#    Simulates evolution of 15,000-copy CEN178 arrays over 6 million generations.
#    Includes downstream analyses such as unit extraction, consensus calling,
#    identity heatmaps, and animation of repeat homogenization and diversification.
#
# Author: xdong@mpipz.mpg.de
###############################################################################

###############################################
# Simulation 1: Non-Allelic Gene Conversion
###############################################

# Step 1: Simulate 500 random point mutations
python ./bin/simulate_SNP_in_bed.py CEN178.repeats.from.F0.fa.bed.fa 500 \
    CEN178.repeats.from.F0.fa.bed.fa.500SNP.fa 500SNP.list

# Step 2: Prepare mutation positions
cat 500SNP.list | grep -v Pos | awk -F ":" '{print $1"\t"$2}' \
    | awk '{print $1"_"$2"_"$3"\t20"}' > list.kmer.count

# Step 3: Generate k-mer search folders
for i in $(cut -f1 list.kmer.count); do
    mkdir -p $i
    echo $i | awk -F "_" '{print $1"_"$2":"$3"\t"$4}' > $i/pos.txt
done

# Step 4: Find longest identical k-mers
for i in $(cut -f1 list.kmer.count); do
    python ./bin/find_longest_all_kmer.py \
        F0.fa.centromere.repeat.500SNP.fa \
        ${i}/pos.txt ${i}/pos.kmer.out 20 20
done

# Step 5: Extract k-mer coordinates
for i in $(cut -f1 list.kmer.count); do
    cd $i
    python ./bin/get_kmer_position.py pos.kmer.out \
        F0.fa.centromere.repeat.500SNP.fa kmer
    cd ..
done

# Step 6: Compute flanking identity length
for i in $(cut -f1 list.kmer.count); do
    cd $i
    for bed in kmer*bed; do
        python ./bin/bed_flank_identical_len.py \
            F0.fa.centromere.repeat.500SNP.fa \
            $bed $bed.iden.len --bed
    done
    cd ..
done

# Step 7: Annotate matches and identify potential donors
for i in $(cut -f1 list.kmer.count); do
    cd $i
    cat *len | awk '{print $0"\t"$6-$5"\t"$5-$2}' \
        | sort | uniq | sort -nrk7,7 -nk8,8 \
        | awk '{if($1==$4){print $0"\tYES"}else{print $0"\tNO"}}' \
        > all.kmer.iden.len
    cd ..
done

# Step 8: Find nearest in-frame match (potential gene conversion donor)
for i in $(cut -f1 list.kmer.count); do
    cd $i
    python ./bin/find_nearest_inframe_kmer.py \
        CEN178.repeats.from.F0.fa.bed \
        all.kmer.iden.len nearest.inframe.20mer.len
    cd ..
done

# Step 9: Collate all results
for i in $(cut -f1 list.kmer.count); do
    cat $i/nearest.inframe.20mer.len | awk -v i=$i '{print i"\t"$0}' \
        >> all.nearest.inframe.20mer.len
done


###############################################
# Simulation 2: Centromeric Mutation Simulation
###############################################

# Step 1: Initial Simulation (1000 generations)
python ./bin/simulate.py \
    ./data/15000copy_cen178.seq 1000 \
    ./data/15000copy_cen178.178bp.bed.pos \
    1.fasta/1000generation.out.fa \
    2.records/1000generation.record.txt \
    3.generation_unit_pos/1000generation.unit.pos

# Step 2: Continue Simulation (to 6 million generations)
for i in $(seq 2000 1000 6000000); do
    prev=$((i - 1000))
    python ./bin/simulate.py \
        1.fasta/${prev}generation.out.fa 1000 \
        3.generation_unit_pos/${prev}generation.unit.pos \
        1.fasta/${i}generation.out.fa \
        2.records/${i}generation.record.txt \
        3.generation_unit_pos/${i}generation.unit.pos
done

# Step 3: Add FASTA headers
for test in $(seq 1 5); do
    for gen in $(seq 1000 1000 6000000); do
        echo ">centro_${gen}gen" | cat - test${test}/1.fasta/${gen}generation.out.fa \
            > tmp && mv tmp test${test}/1.fasta/${gen}generation.out.fa
    done
done

# Step 4: Extract Units
for gen in $(seq 1000 1000 6000000); do
    sed '$d' 3.generation_unit_pos/${gen}generation.unit.pos > tmp1
    sed '1d' 3.generation_unit_pos/${gen}generation.unit.pos > tmp2
    paste tmp1 tmp2 | awk -v gen=$gen '{print "centro_"gen"gen\t"$2"\t"$4}' \
        > 4.generation_unit_fasta/${gen}generation.unit.bed
    seqtk subseq 1.fasta/${gen}generation.out.fa 4.generation_unit_fasta/${gen}generation.unit.bed \
        > 4.generation_unit_fasta/${gen}generation.unit.fa
    rm tmp1 tmp2
done

# Step 5: Generate Consensus and Stats
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

    python ./bin/fasta_stat.py 4.generation_unit_fasta/${gen}generation.unit.fa \
        >> all.generations.fasta.stat

    python ./bin/count_cen178.py 4.generation_unit_fasta/${gen}generation.unit.fa ./data/cen178.fa \
        >> all.generations.cen178.count

    python ./bin/edit_distance_Levenshtein.py \
        4.generation_unit_fasta/${gen}generation.unit.fa ./data/cen178.fa \
        4.generation_unit_fasta/${gen}generation.unit.fa.edit_distance \
        >> all.generations.avg.cen178.edit_distance

    seqkit fx2tab --gc 1.fasta/${gen}generation.out.fa \
        | awk '{print $1"\t"100-$3}' >> all.generations.fasta.at_content
done

# Step 6: Clean paths
sed -i 's/^4.generation_unit_fasta\///g' all.generations.*
sed -i 's/generation.unit.fa//g' all.generations.*

# Step 7: Merge Test Statistics
for test in $(seq 1 5); do
    for file in fasta.stat cen178.count avg.cen178.edit_distance unit.uniq.consensus.fa fasta.at_content; do
        awk -v test=$test '$1 <= 6000000 {print $0"\tTest"test}' \
            test${test}/all.generations.${file} >> 5tests.all.generations.${file}
    done
done

# Step 8: ModDotPlot Heatmaps
source ~/software/ModDotPlot/venv/bin/activate
maxseq=4600000

moddotplot static -f data/15000copy_cen178.seq -w 10000 -id 0 --no-plot -o 6.moddotplot/1.bed_iden0/
python ./bin/heatmap_moddotplot.py \
    6.moddotplot/1.bed_iden0/15000copy_cen178.bed \
    6.moddotplot/3.png/centro_15K_CEN178.heatmap.png \
    --palette-colors '#313695' '#4575B4' '#74ADD1' '#ABD9E9' '#E0F3F8' '#FEE090' '#FDAE61' '#F46D43' '#D73027' '#A50026' \
    --breakpoints 70 73 76 79 82 85 88 91 94 97 100 \
    --width 7 --height 7 -xlab "Array Size (bp)" -ylab "Array Size (bp)" \
    --title 0gen --xlim 0 $maxseq --ylim 0 $maxseq

for gen in $(seq 10000 10000 6000000); do
    moddotplot static -f 1.fasta/${gen}generation.out.fa -w 10000 -id 0 --no-plot -o 6.moddotplot/1.bed_iden0/
    awk '$7 >= 70' 6.moddotplot/1.bed_iden0/centro_${gen}gen.bed > 6.moddotplot/2.bed_iden70/centro_${gen}gen.iden70.bed
    python ./bin/heatmap_moddotplot.py \
        6.moddotplot/2.bed_iden70/centro_${gen}gen.iden70.bed \
        6.moddotplot/3.png/centro_${gen}gen.iden70.heatmap.png \
        --palette-colors '#313695' '#4575B4' '#74ADD1' '#ABD9E9' '#E0F3F8' '#FEE090' '#FDAE61' '#F46D43' '#D73027' '#A50026' \
        --breakpoints 70 73 76 79 82 85 88 91 94 97 100 \
        --width 7 --height 7 -xlab "Array Size (bp)" -ylab "Array Size (bp)" \
        --title ${gen}gen --xlim 0 $maxseq --ylim 0 $maxseq
done

# Step 9: Compile Heatmap Video
ln -s 6.moddotplot/3.png/centro_15K_CEN178.heatmap.png 6.moddotplot/4.ffmpeg/0000.png
for i in $(seq 1 1200); do
    gen=$((i * 5000))
    printf -v padded "%04d" $i
    ln -s 6.moddotplot/3.png/centro_${gen}gen.iden70.heatmap.png 6.moddotplot/4.ffmpeg/${padded}.png
done

ffmpeg -framerate 30 -i 6.moddotplot/4.ffmpeg/%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p \
    6.moddotplot/4.ffmpeg/simulation.mp4

