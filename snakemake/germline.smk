rule all:
    input:
        "results/eval/CHM1-CHM13.hs38DH.bwa-mem.octopus.pass.vcfeval"

rule download_hs38DH:
    output:
        "data/references/hs38DH.fa",
        "data/references/hs38DH.fa.alt"
    container:
        "docker://biocontainers/bwakit:v0.7.15_cv1"
    shell:
        """
        run-gen-ref hs38DH
        mv hs38DH.fa* $(dirname {output[0]})
        """

rule download_fastq:
    output:
        "data/reads/raw/CHM1-CHM13_{strand}.fastq.gz",
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_{strand}.fastq.gz",
    shell:
        "curl {params.url} | gzip > {output}"

rule download_syndip:
    output:
        vcf = "data/truth/CHM1-CHM13.hs38DH.vcf.gz",
        vcf_index = "data/truth/CHM1-CHM13.hs38DH.vcf.gz.tbi",
        bed = "data/truth/CHM1-CHM13.hs38DH.bed.gz",
        bed_index = "data/truth/CHM1-CHM13.hs38DH.bed.gz.tbi",
        tmp = temp("data/truth/CHM-eval.kit")
    params:
        url = "https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar"
    shell:
        '''
        curl -L {params.url} | tar xf -
        mv {output.tmp}/full.38.vcf.gz {output.vcf}
        mv {output.tmp}/full.38.vcf.gz.tbi {output.vcf_index}
        mv {output.tmp}/full.38.bed.gz {output.bed}
        mv {output.tmp}/full.38.bed.gz.tbi {output.bed_index}
        '''

rule generate_hs38DH_chromosomes_bed:
    input:
        "data/references/hs38DH.fa.fai"
    output:
        "data/references/hs38DH.chromosomes.bed"
    shell:
        "head -25 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

rule samtools_index_fasta:
    input:
        "{fasta}"
    output:
         "{fasta}.fai"
    container:
        "docker://mblanche/bwa-samtools"
    shell:
        "samtools faidx {input}"

rule bwa_index:
    input:
        fa = "data/references/{reference}.fa",
        fai = "data/references/{reference}.fa.fai"
    output:
        "data/references/{reference}.fa.amb",
        "data/references/{reference}.fa.ann",
        "data/references/{reference}.fa.bwt",
        "data/references/{reference}.fa.pac",
        "data/references/{reference}.fa.sa"
    container:
        "docker://mblanche/bwa-samtools"
    shell:
        "bwa index {input.fa}"

rule bwa_map:
    input:
        rules.bwa_index.output,
        fa = "data/references/{reference}.fa",
        fai = "data/references/{reference}.fa.fai",
        fq1 = "data/reads/raw/{library}_1.fastq.gz",
        fq2 = "data/reads/raw/{library}_2.fastq.gz"
    output:
        "data/reads/mapped/{library}.{reference}.bwa-mem.bam"
    params:
        rg = r"@RG\tID:{library}\tSM:{library}\tLB:{library}\tPU:Illumina",
        sort_threads = 4,
        sort_memory_per_thread = "4G"
    log:
        "logs/bwa/{library}.{reference}.log"
    threads: 16
    container:
        "docker://mblanche/bwa-samtools"
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.fa} {input.fq1} {input.fq2} | \
          samtools view -bh | \
          samtools sort -@ {params.sort_threads} -m {params.sort_memory_per_thread} -o {output}) \
         2> {log}"

rule samtools_index_bam:
    input:
        "data/reads/mapped/{prefix}.bam"
    output:
        "data/reads/mapped/{prefix}.bam.bai"
    container:
        "docker://mblanche/bwa-samtools"
    shell:
        "samtools index {input}"

rule octopus:
    input:
        reference = "data/references/{reference}.fa",
        bam = "data/reads/mapped/{library}.{reference}.{mapper}.bam",
        bai = "data/reads/mapped/{library}.{reference}.{mapper}.bam.bai",
        bed = "data/references/{reference}.chromosomes.bed"
    output:
        vcf="results/calls/{library}.{reference}.{mapper}.octopus.vcf.gz",
        vcf_index="results/calls/{library}.{reference}.{mapper}.octopus.vcf.gz.tbi"
    params:
        err_model = "PCRF.X10",
        forest = "/opt/octopus/resources/forests/germline.v0.7.4.forest"
    log:
        "logs/octopus/{library}.{reference}.{mapper}.log"
    benchmark:
        "results/benchmarks/octopus/{library}.{reference}.{mapper}.tsv"
    threads: 16
    container:
        "docker://dancooke/octopus"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bam} \
         -t {input.bed} \
         --sequence-error-model {params.err_model} \
         --forest {params.forest} \
         -o {output} \
         --threads {threads} \
         2> {log}"

rule rtg_format:
    input:
        "{prefix}.fa"
    output:
        directory("{prefix}.sdf")
    container:
        "docker://realtimegenomics/rtg-tools"
    shell:
        "rtg format {input} -o {output}"

def _get_score_field(wildcards):
    if wildcards.caller == "octopus":
        return "RFGQ"
    else:
        return "GQ"

rule vcfeval:
    input:
        reference = "data/references/{reference}.sdf",
        baseline_vcf = "data/truth/{library}.{reference}.vcf.gz",
        baseline_vcf_index = "data/truth/{library}.{reference}.vcf.gz.tbi",
        evaluation_regions = "data/truth/{library}.{reference}.bed.gz",
        calls_vcf = "results/calls/{library}.{reference}.{mapper}.{caller}.vcf.gz",
        calls_vcf_index = "results/calls/{library}.{reference}.{mapper}.{caller}.vcf.gz.tbi"
    output:
        directory("results/eval/{library}.{reference}.{mapper}.{caller}.{filter}.vcfeval")
    log:
        "logs/eval/{library}.{reference}.{mapper}.{caller}.{filter}.vcfeval.log"
    params:
        score_field = _get_score_field,
        all_records = lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
        output_mode = "annotate",
        memory = "40g"
    threads: 20
    resources:
        mem_gb = 40
    container:
        "docker://realtimegenomics/rtg-tools"
    shell:
        "(rtg \
            RTG_MEM={resources.mem_gb}g \
            vcfeval \
            -t {input.reference} \
            -b {input.baseline_vcf} \
            --evaluation-regions {input.evaluation_regions} \
            -c {input.calls_vcf} \
            -o {output} \
            -f {params.score_field} \
            --ref-overlap \
            --output-mode {params.output_mode} \
            --threads {threads} \
            {params.all_records} \
            )2> {log}"