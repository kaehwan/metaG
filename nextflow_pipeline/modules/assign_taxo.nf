process METAPHLAN4 {

    tag "$meta.id"
    publishDir "${params.outdir}/out_metaphlan4"

    conda (params.enable_conda ? 'bioconda::metaphlan=4.0.2' : null)
    container 'quay.io/biocontainers/metaphlan:4.0.2--pyhca03a8a_0'

    input:
    tuple val(meta), path(reads)
    path metaphlan4_db

    output:
    tuple val(meta), path("out_metaphlan4/${prefix}.metaphlan4.txt"), emit: metaphlan4_txt

    script:
    prefix = "$meta.id"
    """
    mkdir out_metaphlan4
    metaphlan \\
        ${reads[0]},${reads[1]} \\
        --input_type fastq \\
        --tax_lev s \\
        --bowtie2db ${metaphlan4_db} \\
        --bowtie2out ${prefix}.bowtie2.bz2 \\
        --nproc 8 \\
        --output_file out_metaphlan4/${prefix}.metaphlan4.txt
    """
}

process KRAKEN2 {

    tag "$meta.id"
    publishDir "${params.outdir}/out_kraken2_bracken"

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.2 bioconda::bracken=2.8' : null)
    container 'kaehwan88/taxoprofiling_env'

    input:
    tuple val(meta), path(read)
    path kraken2_db

    output:
    tuple val(meta), path("out_kraken2_bracken/${prefix}.kraken2.report.txt"), emit: kraken2_txt
    tuple val(meta), path("out_kraken2_bracken/${prefix}.bracken.report.txt"), emit: bracken_txt

    script:
    prefix = "$meta.id"
    """
    mkdir out_kraken2_bracken

    kraken2 \\
		--use-names \\
		--db ${kraken2_db} \\
		--report ${prefix}.kraken2.report.txt \\
		--gzip-compressed \\
        --threads 8 \\
		${read} \\
		> ${prefix}.kraken2

	bracken \\
		-d ${kraken2_db} \\
		-i ${prefix}.kraken2.report.txt \\
		-o ${prefix}.bracken.report.txt \\
		-l S

    mv ${prefix}.*.report.txt out_kraken2_bracken
    """
}

process GTDB_PREP {

    tag "${database}"

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container 'ubuntu:20.04'

    input:
    path database

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*")

    script:
    """
    mkdir database
    tar -xzf ${database} -C database --strip 1
    """

}

process GTDBTK {

    tag "$meta.id"
    publishDir "${params.outdir}/out_gtdbtk"

    conda (params.enable_conda ? 'bioconda::gtdbtk=2.3.2' : null)
    container 'quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0'

    input:
    tuple val(meta), path(assembly)
    tuple val(db_name), path("database/*")

    output:
    tuple val(meta), path("out_gtdbtk/${prefix}.gtdbtk.tsv"), emit: tsv

    script:
    prefix = "$meta.id"
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    mkdir out_gtdbtk
    
    gtdbtk classify_wf \\
        --genome_dir . \\
        --prefix "${prefix}" \\
        --out_dir "\${PWD}" \\
        --cpus ${task.cpus} \\
        --pplacer_cpus 1 \\
        --skip_ani_screen \\
        --extension fasta
    
    mv classify/${prefix}.bac120.summary.tsv out_gtdbtk/${prefix}.gtdbtk.tsv
    """
}