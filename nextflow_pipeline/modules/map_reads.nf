process BWA_INDEX {

    tag "$fasta.baseName"

    conda (params.enable_conda ? 'bioconda::bwa=0.7.17' : null)
    container 'quay.io/biocontainers/bwa:0.7.17--h7132678_9'

    input:
    path fasta 

    output:
    path bwa, emit: bwa_index

    script:
    """
    mkdir bwa
    bwa index \\
        -p bwa/${fasta.baseName} \\
        ${fasta}
    """
}

process SAMTOOLS_FAIDX {

    tag "$fasta.baseName"

    conda (params.enable_conda ? 'bioconda::samtools=1.16.1' : null)
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

    input:
    path fasta

    output:
    path samtools, emit: samtools_faidx

    script:
    """
    mkdir samtools
    samtools faidx \\
        ${fasta}

    mv *.fasta* samtools
    """
}

process BWA_MEM {

    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::bwa=0.7.17 bioconda::samtools=1.17 bioconda::samclip=0.4.0' : null)
    container 'kaehwan88/mapping_env'

    input:
    tuple val(meta), path(reads)
    path bwa_index
    path samtools_faidx

    output:
    tuple val(meta), path("out_mapping/*.fastq.gz"), emit: reads

    script:
    prefix = "$meta.id"
    """
    mkdir out_mapping
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    FAIDX=`find -L ./ -name "*.fai" | sed 's/\\.fai\$//'`

    bwa mem \${INDEX} ${reads} > ${prefix}.mapped.sam
	samclip --max 5 --ref \${FAIDX} < ${prefix}.mapped.sam > ${prefix}.clipped.sam
	samtools sort -@ 4 -n -O sam ${prefix}.clipped.sam | samtools fixmate -m -O bam - ${prefix}.clipped.fixmate.bam
	samtools sort -@ 4 -O bam -o ${prefix}.clipped.sorted.bam ${prefix}.clipped.fixmate.bam
	samtools markdup -@ 4 -r -S ${prefix}.clipped.sorted.bam ${prefix}.clipped.sorted.dedup.bam
	samtools view -@ 4 -h -b -f 3 ${prefix}.clipped.sorted.dedup.bam > ${prefix}.clipped.sorted.dedup.concord.bam
	samtools sort -@ 4 -n -O bam -o ${prefix}.clipped.sorted.dedup.concord.sorted.bam ${prefix}.clipped.sorted.dedup.concord.bam
	samtools fastq -1 out_mapping/${prefix}_1.fastq.gz -2 out_mapping/${prefix}_2.fastq.gz ${prefix}.clipped.sorted.dedup.concord.sorted.bam

    rm ${prefix}.*.sam
    rm ${prefix}.*.bam
    """
}

process BWA_SNP_ALIGN {

    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::bwa=0.7.17 bioconda::samtools=1.17 bioconda::samclip=0.4.0' : null)
    container 'kaehwan88/mapping_env'

    input:
    tuple val(meta), path(reads)
    path bwa_index
    path samtools_faidx

    output:
    tuple val(meta), path("out_map2ref/${prefix}.bam"), emit: bam

    script:
    prefix = "$meta.id"
    """
    mkdir out_map2ref
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    FAIDX=`find -L ./ -name "*.fai" | sed 's/\\.fai\$//'`

    bwa mem -M -t 4 \${INDEX} ${reads} > ${prefix}.mapped.sam
	samclip --max 0 --ref \${FAIDX} < ${prefix}.mapped.sam > ${prefix}.clipped.sam
	samtools sort -@ 4 -n -O sam ${prefix}.clipped.sam | samtools fixmate -m -O bam - ${prefix}.clipped.fixmate.bam
	samtools sort -@ 4 -O bam -o ${prefix}.clipped.sorted.bam ${prefix}.clipped.fixmate.bam
	samtools markdup -@ 4 -r -S ${prefix}.clipped.sorted.bam ${prefix}.clipped.sorted.dedup.bam
	samtools view -@ 4 -h -b -f 3 -e '[NM]<=2' ${prefix}.clipped.sorted.dedup.bam > ${prefix}.clipped.sorted.dedup.concord.bam
	samtools sort -@ 4 -O bam -o out_map2ref/${prefix}.bam ${prefix}.clipped.sorted.dedup.concord.bam

    rm ${prefix}.*.sam
    rm ${prefix}.*.bam
    """
}

process BCFTOOLS_MPILEUP {

    tag "$meta.id"
    publishDir "${params.outdir}/out_variant"

    conda (params.enable_conda ? 'bioconda::samtools bioconda::bamtools bioconda::freebayes bioconda::bedtools bioconda::vcflib bioconda::rtg-tools bioconda::bcftools conda-forge::matplotlib' : null)
    container 'kaehwan88/variantcalling_env'

    input:
    tuple val(meta), path(bam)
    path bwa_index

    output:
    tuple val(meta), path("out_variant/${prefix}.vcf"), emit: vcf

    script:
    prefix = "$meta.id"
    """
    mkdir out_variant
    FAIDX=`find -L ./ -name "*.fai" | sed 's/\\.fai\$//'`

    bamtools index -in ${bam}
    bcftools mpileup --max-depth 400 -Ob -a "AD,ADF,ADR,DP,SP,INFO/ADF,INFO/ADR" -f \${FAIDX} ${bam} | bcftools call -A -mv -Ov > ${prefix}.tmp.vcf
    bgzip ${prefix}.tmp.vcf
    tabix -p vcf ${prefix}.tmp.vcf.gz
    rtg vcffilter -q 30 -d 10 --snps-only -i ${prefix}.tmp.vcf.gz -o ${prefix}.vcf.gz
    bgzip -d ${prefix}.vcf.gz

    mv ${prefix}.vcf out_variant
    """
}

process MINIMAP2_NANOPORE {

    tag "$meta.id"

    conda (params.enable_conda ? 'bioconda::coverm=0.6.1' : null)
    container 'quay.io/biocontainers/coverm:0.6.1--h07ea13f_6'

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("out_minimap2/${prefix}.fastq"), emit: reads

    script:
    prefix = "$meta.id"
    """
    mkdir out_minimap2

    # basic mapping using minimap2
    minimap2 -a -x map-ont ${ref} ${reads} > ${prefix}.sam

    # compressing, sorting, indexing with samtools
    samtools sort -@ 16 -m 5G -o ${prefix}.tmp1.bam ${prefix}.sam
    samtools index ${prefix}.tmp1.bam

    # access only the mapped reads
    samtools view -b -F 4 ${prefix}.tmp1.bam > ${prefix}.tmp2.bam

    # extract mapped reads with at least 100% alignment & 80% identity
    coverm filter -b ${prefix}.tmp2.bam -o ${prefix}.bam --min-read-aligned-percent 100 --min-read-percent-identity 80 --threads 16

    # bam to fastq
    samtools fastq -@ 16 ${prefix}.bam > out_minimap2/${prefix}.fastq
    """
}
