# Genome Assembly and Centromeric Mutation Analysis Pipelines

This repository contains five modular, reproducible Bash pipelines designed for a comprehensive study of genome assembly, annotation, mutation identification, and mutation simulation in *Arabidopsis thaliana*. The scripts automate key steps from raw assembly to simulated evolution of tandem repeat arrays.

## Repository Structure

```
project-root/
├── 1.assembly.sh               # Genome assembly and evaluation
├── 2.annotation.sh            # Repeat annotation (CEN178, rDNA, SSRs)
├── 3.assembly_difference.sh   # Assembly comparison, validation, and correction
├── 4.cen_mutation.sh          # Centromeric mutation analysis and alignment refinement
├── 5.simulation.sh            # Simulation of gene conversion and tandem repeat evolution
├── bin/                       # Contains Python and R scripts used in the pipeline
└── data/                      # Input reference genomes, annotations, and simulation templates
```


## Pipelines Overview

### 1. `1.assembly.sh`: Genome Assembly and Evaluation
- De novo assembly with HiFiasm, IPA, Peregrine, Canu, Flye
- Reference-based scaffolding
- Assembly quality assessment using NG50, Merqury, and BUSCO
- Consensus genome generation for downstream mutation comparison

### 2. `2.annotation.sh`: Repeat Annotation
- CEN178 annotation using TRASH
- rDNA and telomeric repeat detection with RepeatMasker and a custom library
- Simple sequence repeat annotation using Arabidopsis-specific repeat libraries

### 3. `3.assembly_difference.sh`: Assembly Comparison and Validation
- Structural variant calling using SYRI
- Misassembly validation using Illumina and HiFi reads based alignment
- Error correction with Pilon, DeepVariant, pbsv, and Sniffles

### 4. `4.cen_mutation.sh`: Centromeric Mutation Analysis
- Word-based alignment for optimal matching in repeat-dense regions
- Mutation left-alignment and false-positive filtering
- HOR score analysis and in-frame mutation pattern checking

### 5. `5.simulation.sh`: Simulations of Centromeric Evolution
- Gene conversion simulation: Introduces random point mutations and detects non-allelic conversion events via k-mer overlap
- Tandem repeat mutation simulation: Evolves a 15,000-copy CEN178 array over 6 million generations, followed by analysis of homogenization, consensus generation, heatmap creation, and video animation

## Dependencies

Each script has its own software dependencies, which are listed at the top of the respective file. Common tools and packages include:

- Genome assemblers: HiFiasm, IPA, Peregrine, Canu, Flye
- Annotation tools: TRASH, RepeatMasker
- Variant analysis: SYRI, Pilon, DeepVariant, Sniffles, pbsv
- Supporting scripts: Python, R, Perl (located in `bin/`)
- Custom input files: (located in `data/`)

## Usage

Make sure all dependencies are installed and paths are properly configured. Then, execute each script in order or independently as needed:

<!--
## Citation

Please cite the following work when using this repository or any of its components in your research:

> Dong, X. et al. "Title of preprint." *bioRxiv*, https://doi.org/xxxxxxx
-->

## Contact

For questions, please contact: xdong@mpipz.mpg.de

