include: "rules/common.smk"

##### Target rules #####
rule all:
    input:
        expand("{genome}.bwt", genome=units["reference"].unique()),
        expand("{genome}.fai", genome=units["reference"].unique()),
        expand("results/8_regions/{sample}.reg.tsv", sample=units["sample"].unique()),
        expand("results/0_fastqc_init/{prefix}.qc_init.html", prefix=units.prefix),
        expand("results/2_fastqc_trim/{prefix}.qc_trim.html", prefix=units.prefix),
        "results/stats.xlsx"


##### Modules #####
include: "rules/dopseq.smk"