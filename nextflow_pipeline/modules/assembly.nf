process SPADES {

    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::spades=3.15.4' : null)
    container 'quay.io/biocontainers/spades:3.15.4--h95f258a_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("out_spades/${prefix}.scaffolds.fasta"), emit: assembly

    script:
    prefix = "$meta.id"
    """
    mkdir out_spades
    spades.py \\
        --careful \\
        -o out_spades/${prefix}/ \\
        -1 ${reads[0]} \\
        -2 ${reads[1]}
    
    rm_short_contigs.py out_spades/${prefix}/scaffolds.fasta 1000 > out_spades/${prefix}.scaffolds.fasta
    """
}

process METAWRAP {

    tag "$meta.id"
    conda (params.enable_conda ? 'ursky::metawrap-mg' : null)
    container 'kaehwan88/metawrap_env'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("out_mag/*.fasta"), emit: mag

    script:
    prefix = "$meta.id"
    """
    mkdir raw_fastq
    mkdir out_mag
    pigz -dkc ${reads[0]} > raw_fastq/${prefix}_1.fastq
    pigz -dkc ${reads[1]} > raw_fastq/${prefix}_2.fastq

    # Assembly using MetaSPAdes, followed by MegaHIT
	metawrap assembly -1 raw_fastq/${prefix}_1.fastq -2 raw_fastq/${prefix}_2.fastq -t 16 -m 256 --metaspades -o out_assembly_${prefix}

	# Binning assemblies using metabats, maxbin2, concoct
	metawrap binning -o out_binning_${prefix} -a out_assembly_${prefix}/final_assembly.fasta -t 16 -m 256 --metabat2 --maxbin2 raw_fastq/${prefix}_*.fastq

	# Consolidate bin sets using bin_refinement module
	metawrap bin_refinement -o out_refinement_${prefix} -t 16 -m 256 -A out_binning_${prefix}/metabat2_bins/ -B out_binning_${prefix}/maxbin2_bins/ -c 70 -x 10

	# Improve bins with reassembly: Collect reads belonged to each bin, reassemble them separately with a 'permissive' and 'strict' algorithm
	metawrap reassemble_bins -o out_reassembly_${prefix} -1 raw_fastq/${prefix}_1.fastq -2 raw_fastq/${prefix}_2.fastq -t 16 -m 256 -c 70 -x 10 -b out_refinement_${prefix}/metawrap_70_10_bins

    for i in out_reassembly_${prefix}/reassembled_bins/*.fa
    do
        newName=\$(echo \$(basename \${i%.fa}))
        mv \${i} out_mag/${prefix}.\${newName}.fasta
    done
    """
}

process FLYE {

    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::flye=2.9' : null)
    container 'quay.io/biocontainers/flye:2.9--py39h39abbe0_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("out_flye/${prefix}.fasta"), emit: assembly

    script:
    prefix = "$meta.id"
    """
    mkdir out_flye

    flye \\
        --nano-hq ${reads} \\
        --out-dir out_flye \\
        -g 2.1m \\
        -t 16

    mv out_flye/assembly.fasta out_flye/${prefix}.fasta
    """
}

process MEDAKA {

    tag "$meta1.id"
    conda (params.enable_conda ? 'bioconda::medaka=1.4.4' : null)
    container 'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0'

    input:
    tuple val(meta1), path(reads)
    tuple val(meta2), path(assembly)

    output:
    tuple val(meta1), path("out_medaka/${prefix}.fasta"), emit: assembly

    script:
    prefix = "$meta1.id"
    """
    mkdir out_medaka

    medaka_consensus \\
        -i ${reads} \\
        -d ${assembly} \\
        -o out_medaka \\
        -t 16 \\
        -m r941_min_sup_g507

    mv out_medaka/consensus.fasta out_medaka/${prefix}.fasta
    """
}
