rule fastqc_init:
    input:
        get_fastq
    output:
        flag="results/0_fastqc_init/{sample}-{unit}.done"
    params:
        tempdir="work/fastqc_init_{sample}_{unit}",
        outdir="results/0_fastqc_init/"
    conda:
        "../env.yaml"
    log:
        "results/logs/fastqc/{sample}-{unit}_init.log"
    threads:
        config["params"]["threads"]
    shell:
        "rm -rf {params.tempdir} && "
        "mkdir -p {params.tempdir} && "
        "fastqc"
        " --threads {threads}"
        " --outdir {params.tempdir}"
        " {input}"
        " &> {log} && "
        "mv {params.tempdir}/*html {params.outdir} && "
        "mv {params.tempdir}/*zip {params.outdir} && "
        "touch {output.flag} && rm -r {params.tempdir}"

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
    threads: 
        config["params"]["threads"]
    conda:
        "../env.yaml"
    shell:
        "cutadapt"
        " {params}"
        " -j {threads}"
        " -o {output.fastq}"
        " {input}"
        " > {output.qc} 2> {log}"

rule trim_reads_pe:
    input:
        get_fastq
    output:
        fastq1="results/1_trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/1_trimmed/{sample}-{unit}.2.fastq.gz",
        stats="results/1_trimmed/{sample}-{unit}.trim.txt"
    params:
        ampl_to_cutadapt_pe, config["params"]["cutadapt"]["pe"]
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    threads: 
        config["params"]["threads"]
    conda:
        "../env.yaml"
    shell:
        "cutadapt"
        " {params}"
        " -j {threads}"
        " -o {output.fastq1}"
        " -p {output.fastq2}"
        " {input}"
        " > {output.stats} 2> {log}"

rule fastqc_trimmed:
    input:
        get_trimmed_reads
    output:
        flag="results/2_fastqc_trim/{sample}-{unit}.done"
    params:
        tempdir="work/fastqc_trim_{sample}_{unit}",
        outdir="results/2_fastqc_trim/"
    conda:
        "../env.yaml"
    log:
        "results/logs/fastqc/{sample}-{unit}_trim.log"
    threads: 
        config["params"]["threads"]
    shell:
        "rm -rf {params.tempdir} && "
        "mkdir -p {params.tempdir} && "
        "fastqc"
        " --threads {threads}"
        " --outdir {params.tempdir}"
        " {input}"
        " &> {log} && "
        "mv {params.tempdir}/*html {params.outdir} && "
        "mv {params.tempdir}/*zip {params.outdir} && "
        "touch {output.flag} && rm -r {params.tempdir}"

# genome preparation done in the input dir
rule bwa_index:
    input:
        "{genome}"
    output:
        # "{genome}.amb",
        # "{genome}.ann",
        "{genome}.bwt"
        # "{genome}.pac",
        # "{genome}.sa"
    conda:
        "../env.yaml"
    log:
        "{genome}.bwa_index.log"
    shell:
        "bwa index"
        " -p {input}"
        " {input} &> {log}"

rule samtools_faidx:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    conda:
        "../env.yaml"
    shell:
        "samtools faidx {input}"

rule map_reads_bwa_mem:
    input:
        reads=get_trimmed_reads,
        index=ancient(get_ref_bwt)
    output:
        bam="results/3_mapped/{sample}-{unit}.sorted.bam"
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
        " -t {threads}"
        " {params.extra}"
        " {params.index}"
        " {input.reads}"
        " 2> {log.mem} | "
        "samtools sort - "
        " -o {output} &> {log.sort}"

# alignment filtering and merging
rule mark_duplicates:
    input:
        "results/3_mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("results/4_dedup/{sample}-{unit}.dedup.bam"),
        metrics="results/4_dedup/{sample}-{unit}.dedup.txt"
    log:
        "results/logs/picard_dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard_MarkDuplicates"]
    conda:
        "../env.yaml"
    shell:
        "picard MarkDuplicates {params} "
        " REMOVE_DUPLICATES=true "
        " INPUT={input} OUTPUT={output.bam} "
        " METRICS_FILE={output.metrics} "
        " &> {log}"

rule samtools_filter:
    input:
        "results/4_dedup/{sample}-{unit}.dedup.bam" 
        # if config["workflow"]["do_rmdup"] else 
        # "results/3_mapped/{sample}-{unit}.sorted.bam")    
    output:
        bam="results/5_filtered/{sample}-{unit}.filter.bam",
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
        " {output} {input}"

rule regions:
    input:
        "results/6_merged/{sample}.bam",
        genome_fai=get_ref_fai,
    output:
        pos="results/7_positions/{sample}.bed",
        reg="results/8_regions/{sample}.reg.tsv"
    params:
        sample="{sample}",
        do_plot_reg=config["workflow"]["do_plot_reg"],
        plot="results/8_regions/{sample}.reg.pdf",
        plot_ncols=config["params"]["region"]["plot_ncols"],
        plot_chrom_height=config["params"]["region"]["plot_chrom_height"],
        plot_chrom_width=config["params"]["region"]["plot_chrom_width"],
        chrom_list=get_chrom_list
    conda:
        "../env.yaml"
    script:
        "../scripts/regions.py"

rule stats: 
    input:
        get_position_beds
    output:
        "results/stats.xlsx"
    params:
        samples=config["samples"],
        trim="results/1_trimmed",
        dedup="results/4_dedup", #if config["workflow"]["do_rmdup"] else None),
        flt="results/5_filtered",
        pos="results/7_positions/"
    conda:
        "../env.yaml"
    script:
        "../scripts/stats.py"