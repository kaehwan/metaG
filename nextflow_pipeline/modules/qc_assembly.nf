process QUAST {

    tag "$meta.id"
    publishDir "${params.outdir}"

    conda (params.enable_conda ? 'bioconda::quast=5.0.2' : null)
    container 'quay.io/biocontainers/quast:5.0.2--py37pl5321hfecc14a_6'

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("out_quast/${prefix}.quast.tsv"), emit: tsv

    script:
    prefix = "$meta.id"
    """
    mkdir out_quast

    quast.py \\
		-o out_quast/${prefix} \\
		-t 16 \\
        --space-efficient \\
		--pe1 ${reads[0]} \\
        --pe2 ${reads[1]} \\
        ${assembly}

    mv out_quast/${prefix}/report.tsv out_quast/${prefix}.quast.tsv
    """
}

process QUAST_NANOPORE {

    tag "$meta.id"
    publishDir "${params.outdir}"

    conda (params.enable_conda ? 'bioconda::quast=5.0.2' : null)
    container 'quay.io/biocontainers/quast:5.0.2--py37pl5321hfecc14a_6'

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("out_quast/${prefix}.quast.tsv"), emit: tsv

    script:
    prefix = "$meta.id"
    """
    mkdir out_quast

    quast.py \\
		-o out_quast/${prefix} \\
		-t 16 \\
        --space-efficient \\
		--nanopore ${reads} \\
        ${assembly}

    mv out_quast/${prefix}/report.tsv out_quast/${prefix}.quast.tsv
    """
}