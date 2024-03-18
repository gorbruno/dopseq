# dopseq
Isolated chromosome sequencing analysis

## Introduction

dopseq is a pipeline for reference-based prediction of regions present on chromosomes based on high-throughput sequencing data generated from isolated (flow sorted or microdissected) chromosomes. This project is a reimplementation of https://github.com/ilyakichigin/DOPseq_analyzer.  

Current version includes chromosomal region prediction pipeline: 
- read trimming and qc with cutadapt and fastqc,
- alignment to reference genome with bwa mem,
- PCR duplicate and quality filtering with picard and samtools,
- genome segmentation with DNAcopy R package, plotting and statistics accumulation.

This software relies on [Snakemake](https://snakemake.readthedocs.io/en/v5.2.0/) for workflow and on [conda](https://conda.io/docs/) for dependencies management. The pipeline implementation is based on [Snakemake workflow for dna-seq](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling).

## Quick start

Clone dopseq and go to its directory

```
git clone git@github.com:lca-imcb/dopseq.git
cd dopseq
```

Create environment for the pipeline - we suggest using mamba instead of conda to speed up dependencies resolution
```
mamba env create -f env.yaml
```

If environment resolution fails or some dependencies don't work correctly, you can try using explicit specification files available for Linux and OS X 64-bit

```
conda create --name dopseq --file spec-linux-64.txt
```

or

```
conda create --name dopseq --file spec-osx-64.txt
```


Activate the environment

```
conda activate dopseq
```

Set up your analysis by editing the following files:

- `samples.tsv` for per-sample data, 
- `config.yaml` for shared analysis parameters. 

After that, you can test the pipeline with a dry run

```
snakemake -n
```

And then run the analysis 

```
snakemake -j 4
```

## Parameter setting

### config.yaml

Configuration applied to the entire analysis in YAML format. Sections:

- `samples` - path to tab-separated file with sample data.

- `workflow` - amend workflow logic:
  - `do_plot_reg` - produce ditance-region plots for seqids specified in `chrom_list` column of samples tsv. 

- `params` section provides software-related parameters:
  - `threads` - number of threads for each of `cutadapt`, `bwa mem` and `samtools` processes. 
  You should specify number of cores used by the pipeline with snakemake `-j` parameter.
  - `cutadapt` - cutadapt general filtering options for paired-end (PE) and single-end (SE) reads. Default - trim terminal Ns and remove reads shorter than 20 bp.
  - `picard_MarkDuplicates` - by default set to remove duplicates instead of marking them.
  - `region` - ditance-region plots size parameters

### samples.tsv

Tab-separated file with sample data. Columns:

- `sample` - sample name used as prefix for the output files.

- `unit` - multiple input units (lanes, libraries, biological replicates etc) can be specified for each sample. 
Each unit is trimmed, aligned and filtered separately, so file prefixes are `{sample}-{unit}` for steps prior to bam merging.
After bam merging, file prefixes only include `{sample}`.

- `platform` - sequencing platform, most likely `ILLUMINA`. used to set `@RG PL` field of the per-unit BAMs, merged BAM contains all `@RG` lines. 

- `adapters` - sequencing or library adapters to be removed from reads:
  - `dop` - DOP-PCR MW6 primer (Telenius et al., 1992), reads without primer match at 5\` end are discarded;
  - `dop_relaxed` - same as `dop`, but reads without primers are not discarded;
  - `wga` - WGA1 GP primer, reads without primer match at 5\` end are discarded;
  - `wga_relaxed` - same as `wga`, but reads without primers are not discarded;
  - `illumina` - only remove Illumina standard adapter from 3\` end (readthrough);
  - `none` - do not trim adapters, just filter reads using parameters from the configuration file

For paired-end reads trimming with `dop` and `wga`, only pairs with primer matches in both reads are retained. 
To increase the amount of retained reads, you can specify forward and reverse reads as separate units of one sample.

- `reference` - path to unpacked reference genome in fasta format.  
Note that index files will be saved in the same directory as the reference itself.
For non-model species, selection of the reference balances between evolutionary proximity to the sample species and assembly quality. 

- `min_len` - minimum alignment length to be retained for region inference. 

- `min_q` - minimum mapping quality of the alignment to be retained for region inference.

- `chrom_list` - list of sequids and their lengths to be visualised in ditance-region plots - fai files are accepted

- `fq1` - forward or single-end reads fastq file (can be plain or gzipped).

- `fq2` - reverse reads fastq file. Keep blank for single-end reads.


## Pipeline steps

Output files are stored in the `dopseq/results` directory, 
with subdirectories numbered according to the analysis order and named by the output type. 

- `0_fastqc_init` - FastQC assessment of the input reads;
- `1_trimmed` - reads trimmed and filtered with cutadapt;
- `2_fastqc_trim` - FastQC assessment of the trimmed reads, helps to identify the remaining problems with reads;
- `3_mapped` - reads mapped to reference genome with bwa mem;
- `4_dedup` - deduplicated read alignments with Picard MarkDuplicates;
- `5_filtered` - alignments filtered for mapping quality and aligned fragment length using samtools and awk;
- `6_merged` - merging all unit BAMs per sample with samtools;
- `7_positions` - merging overlapping reads into read positions with pybedtools;
- `8_regions` - genome segmentation based on distances between read positions with DNAcopy;
- `stats.xlsx` - statistics for sequencing, mapping, filtering, and positions.

## Output interpretation 

The key outputs of the pipeline are per-sample genome segmentations located in 
`results/8_regions/{sample}.tsv` files. 
The segmentation is based on the end-to-start distances between consecutive genome positions covered by reads, or PD (Pairwise Distances). 
The mean values of PD per segment are summarized in `pd_mean` column. 
Regions with lower `pd_mean` are more likely to be present on sampled chromosome and will be subsequently called _target regions_. 
Regions with higher `pd_mean` tend to represent a noise that can result from erroneous mappings or contamination with whole-genome or external DNA. 
Thus, problem of target regions identification can formulated as the problem of drawing the borderline between low `pd_mean` target regions and high `pd_mean` contamination. 
Currently, we address this problem by testing for significance of position density enrichment relative to genome-averaged position density.


Still, the step of target regions identification is not automated as there are multiple confounding factors:

- additional filters based on mean position size or mean position coverage for mapped positions can be useful in some cases;
- segmentation for small contigs/scaffolds is unreliable, we often used the minimum scaffold size threshold of 50 kbp, but it cannot be recommended for all cases;
- for highly fragmented reference genome assemblies there can be no obvious distinction between high-PD and low-PD regions;
- same problem arises with increase of evolutionary distance between sampled and reference species;
- some target regions are over-segmented due to mapping efficiency varying along the chromosome and should be joined upon visual inspection;
- sometimes, region margins should be manually corrected.

General recommendation is to sort region prediction files by `pd_mean`, 
set up filters that fit your data 
(max `p_value`, min `reg_pos`, min `chrom_len`, 
optionally - min `pos_len_mean` and min `prop_reg_covered`), 
and extract top most significant target regions by looking for significant gap in `pd_mean` 
(e.g., several target regions have 1-5 kbp `pd_mean`, 
then first non-target regions has 20 kbp `pd_mean`, 
then values increase gradually). 
For manual correction of region margins, visualization of BAMs (step 6) and BEDs (step 7) in genome browser can be useful. 

## Notes

This section can be useful if you want to change parameters not listed above, as well as to add or remove steps.

`Snakemake` file sets the desired output files and links to the other smk files: 
- `rules/dopseq.smk` includes rules for the pipeline itself;
- `rules/common.smk` includes functions for filename and parameter setting based on sample data, which is located in `samples.tsv`;
For further reading, please refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/v5.2.0/).

`script` directory: 
- `script/regions.py` converts filtered BAM to BED with overlapping reads merged into positions (`results/7_positions`) and performs genome segmentation (`results/8_regions`).
- `script/stats.py` collects statistics per unit and per sample in a single xlsx file.

`env.yaml` lists all conda packages and their versions. 
A single environment is generated for all tasks except fastqc 
(which relies on smk wrapper).

`schemas` directory includes files for testing `config.yaml` and `samples.tsv` format.

`test` directory contains tests of dopseq functionality, including whole pipeline run on fox B chromosome data and sanity checks of the output files. Test can be run with 
`cd dopseq && bash test/test.sh` - if you already created the `dopseq` environment, do `touch test/env.installed` first.

## Citation

- Makunin AI, Kichigin IG, Larkin DM, O’Brien PCM, Ferguson-Smith MA, Yang F, et al. Contrasting origin of B chromosomes in two cervids (Siberian roe deer and grey brocket deer) unravelled by chromosome-specific DNA sequencing. BMC Genomics. 2016;17: 618. doi:10.1186/s12864-016-2933-6
