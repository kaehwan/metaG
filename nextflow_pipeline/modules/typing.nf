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

process VF_AMR {

    tag "$meta.id"
    publishDir "${params.outdir}"

    conda (params.enable_conda ? 'bioconda::abricate=1.0.1' : null)
    container 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("out_vfamr/${prefix}.vf.tsv"), emit: vf
    tuple val(meta), path("out_vfamr/${prefix}.amr.tsv"), emit: amr

    script:
    prefix = "$meta.id"
    """
    mkdir out_vfamr
    abricate \\
        --db vfdb \\
        --quiet \\
        $assembly \\
        > out_vfamr/${prefix}.vf.tsv
    abricate \\
        --db ncbi \\
        --quiet \\
        $assembly \\
        > out_vfamr/${prefix}.amr.tsv
    """
}