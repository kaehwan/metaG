process MLST {

    tag "$meta.id"
    publishDir "${params.outdir}"

    conda (params.enable_conda ? 'bioconda::mlst=2.22.0' : null)
    container 'quay.io/biocontainers/mlst:2.22.0--hdfd78af_0'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("out_mlst/${prefix}.mlst.tsv"), emit: mlst

    script:
    prefix = "$meta.id"
    """
    mkdir out_mlst
    mlst \\
        --legacy \\
        --scheme sagalactiae \\
        ${assembly} \\
        > out_mlst/${prefix}.mlst.tsv
    """
}