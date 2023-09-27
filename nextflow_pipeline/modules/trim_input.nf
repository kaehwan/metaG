process FASTP {

    tag "$meta.id"

    conda (params.enable_conda ? 'bioconda::fastp=0.20.1' : null)
    container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("out_fastp/*_{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("out_fastp/*.fastp.json")    , emit: json
    tuple val(meta), path("out_fastp/*.fastp.html")    , emit: html

    script:
    def prefix = "$meta.id"
    """
    mkdir out_fastp
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o out_fastp/${prefix}_1.fastq.gz \\
        -O out_fastp/${prefix}_2.fastq.gz \\
        -q 30 \\
        --json out_fastp/${prefix}.fastp.json \\
        --html out_fastp/${prefix}.fastp.html \\
        --thread 2 \\
        --detect_adapter_for_pe \\
        --overrepresentation_analysis \\
        --correction \\
        --cut_right
    """
}