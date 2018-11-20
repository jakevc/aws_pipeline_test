"""Align and call peaks on nine sample fastqs with 10K paired-end reads.""" 

# import pre-configured s3 remote
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(access_key_id="AKIAI2SP5BM5QZN5RQLQ", secret_access_key="KUntFWs0ZA753oQw5zlY1ggW4348TVdLDfdZxLdD")

# samples 1-9
SAMPLES =  ["sample_{}".format(num) for num in range(1,9)]

# access remote genome
GENOME = S3.remote("jv2245-test-data/genome/test-genome.fa.gz") 

rule all: 
    input:
        S3.remote(expand("jv2245-test-data/peaks/{sample}_peaks.broadPeak", sample=SAMPLES))

rule bwa_map:
    input:
        genome = GENOME,
        r1 = S3.remote("jv2245-test-data/fastq_test/{sample}_1.fastq"),
        r2 = S3.remote("jv2245-test-data/fastq_test/{sample}_2.fastq") 
    output:
        temp("data/align/bwa/{sample}.sam" )
    threads: 4
    conda: 
	"envs/bwa.yml"
    params:
        rg = "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina" 
    shell:
        """
        bwa mem -t {threads} -M -R '{params.rg}' {input.genome} {input.r1} {input.r2} > {output}
        """
        
rule sam_to_bam: 
    input:
        rules.bwa_map.output
    output: 
        temp("data/align/bwa/{sample}.bam")
    threads: 4    
    conda: 
	"envs/samtools.yml"
    shell:
        """
        samtools view -@ {threads} -Sb {input} -o {output}
        """

rule sort:
    input:
        rules.bwa_map.output
    output:
        "data/align/sort/{sample}.sorted.bam"
    threads: 4     
    conda: 
	"envs/samtools.yml"
    shell:
        "samtools sort -@ {threads} -T {wildcards.sample} -o {output} {input}"

rule index:
    input:
        rules.sort.output
    output:
        "data/align/sort/{sample}.sorted.bam.bai"
    conda:
	"envs/samtools.yml"
    shell: 
        "samtools index {input} {output}" 
        
###########  Call Peaks  #############
rule macs2:
    input:
        "data/align/sort/{sample}.sorted.bam"
    output:
        S3.remote("jv2245-test-data/peaks/{sample}_peaks.broadPeak")
    conda: 
	"envs/macs2.yml"
    shell:
        """
        macs2 callpeak \
                    -t {input} \
                    -n {wildcards.sample} \
                    --outdir S3.remote("jv2245-test-data/peaks") \
                    -g 1.3e9 \
                    -q 0.01 \
                    -f BAMPE \
                    --broad \
                    --nomodel \
        """
