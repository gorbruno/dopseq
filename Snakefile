include: "rules/common.smk"

##### Target rules #####
rule all:
    input:
        expand("results/8_regions/{sample}.tsv", sample=units["sample"].unique()),
        expand("results/0_fastqc_init/{prefix}.html", prefix=units.prefix),
        expand("results/2_fastqc_trim/{prefix}.html", prefix=units.prefix),
        "results/stats.xlsx"
        # "results/2_fastqc_trim/{sample}-{unit}.html".format(units.index)
        # "results/qc/multiqc.html"
        # 'trimmed/A-lane2.fastq.gz',
        # 'trimmed/A-lane1.1.fastq.gz',
        # 'trimmed/A-lane1.2.fastq.gz'
        # "annotated/all.vcf.gz",
        # "tables/calls.tsv.gz",
        # "plots/depths.svg",
        # "plots/allele-freqs.svg"


##### Modules #####
include: "rules/dopseq.smk"
# include: "rules/qc.smk"
