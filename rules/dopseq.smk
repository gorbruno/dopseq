# trimming reads
rule fastqc_init:
    input:
        get_fastq
    output:
        html="results/0_fastqc_init/{sample}-{unit}.html",
        zip="results/0_fastqc_init/{sample}-{unit}.zip"
    # params:
    #     dir="results/0_fastqc_init/"
    # conda:
    #     "../env.yaml"
    wrapper:
        "0.27.1/bio/fastqc"

rule trim_reads_se:
    input:
        get_fastq
    output:
        fastq="results/1_trimmed/{sample}-{unit}.fastq.gz",
        qc="results/1_trimmed/{sample}-{unit}.trim.txt"
    params:
        ampl_to_cutadapt_se, config["params"]["cutadapt"]["se"]
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    conda:
        "../env.yaml"
    shell:
        "cutadapt"
        " {params}"
        " -o {output.fastq}"
        " {input}"
        " > {output.qc} 2> {log}"

rule trim_reads_pe:
    input:
        get_fastq
    output:
        fastq1="results/1_trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/1_trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/1_trimmed/{sample}-{unit}.trim.txt"
    params:
        ampl_to_cutadapt_pe, config["params"]["cutadapt"]["pe"]
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    # conda:
    #     "../env.yaml"
    shell:
        "cutadapt"
        " {params}"
        " -o {output.fastq1}"
        " -p {output.fastq2}"
        " {input}"
        " > {output.qc} 2> {log}"

rule fastqc_trim:
    input:
        get_trimmed_reads
    output:
        html="results/2_fastqc_trim/{sample}-{unit}.html",
        zip="results/2_fastqc_trim/{sample}-{unit}.zip"
    # params:
    #     dir="results/0_fastqc_init/"
    # conda:
    #     "../env.yaml"
    wrapper:
        "0.27.1/bio/fastqc"

# genome preparation and alignment
rule bwa_index:
    input:
        "{genomepath}"
    output:
        "{genomepath}.amb",
        "{genomepath}.ann",
        "{genomepath}.bwt",
        "{genomepath}.pac",
        "{genomepath}.sa"
    conda:
        "../env.yaml"
    shell:
        "bwa index"
        " -p {input}"
        " {input}"

rule samtools_faidx:
    input:
        "{genomepath}"
    output:
        "{genomepath}.fai"
    conda:
        "../env.yaml"
    shell:
        "samtools faidx {input}"

rule map_reads_bwa_mem:
    input:
        reads=get_trimmed_reads
    output:
        "results/3_mapped/{sample}-{unit}.sorted.bam"
    log:
        mem="results/logs/bwa_mem/{sample}-{unit}.log",
        sort="results/logs/samtools_sort/{sample}-{unit}.log",
    params:
        index=get_ref,
        extra=get_read_group
    threads: config["params"]["threads"]
    conda:
        "../env.yaml"
    shell:
        "bwa mem"
        " -t {threads} "
        "{params.extra} "
        "{params.index} "
        "{input.reads} "
        "2> {log.mem} | "
        "samtools sort - "
        "-o {output} &> {log.sort}"

# alignment filtering and merging
rule mark_duplicates:
    input:
        "results/3_mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("results/4_dedup/{sample}-{unit}.bam"),
        metrics="results/4_dedup/{sample}-{unit}.dedup.txt"
    log:
        "results/logs/picard/dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    conda:
        "../env.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} "
        "&> {log}"

rule samptools_filter:
    input:
        # get_dedup_bams
        "results/4_dedup/{sample}-{unit}.bam"
    output:
        bam="results/5_filtered/{sample}-{unit}.bam",
        metrics="results/5_filtered/{sample}-{unit}.filter.txt"
    params:
        minq=get_min_q,
        minl=get_min_len
    conda:
        "../env.yaml"
    shell:
        "samtools view -h -q {params.minq} {input} | "
        "awk 'length($10) > {params.minl} || $1 ~ /^@/' | "
        "samtools view -bS - > {output.bam}; "
        "samtools stats {output.bam} > {output.metrics}"

rule samtools_merge:
    input:
        get_filtered_bams
    output:
        "results/6_merged/{sample}.bam"
    params:
        "" 
    threads: config["params"]["threads"]
    conda:
        "../env.yaml"
    shell:
        "samtools merge --threads {threads} {params} "
        "{output} {input}"

# regions
rule regions:
    input:
        "results/6_merged/{sample}.bam",
        genome_fai=get_ref_fai
    output:
        pos="results/7_positions/{sample}.bed",
        reg="results/8_regions/{sample}.tsv"
    params:
        sample="{sample}"
    conda:
        "../env.yaml"
    script:
        "../scripts/regions.py"

# statistics
rule stats:
    input:
        get_position_beds
    output:
        "results/stats.xlsx"
    params:
        samples=config["samples"],
        trim="results/1_trimmed",
        dedup="results/4_dedup",
        flt="results/5_filtered",
        pos="results/7_positions/"
    conda:
        "../env.yaml"
    script:
        "../scripts/stats.py"