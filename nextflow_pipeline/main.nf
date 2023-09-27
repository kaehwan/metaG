#!/usr/bin/env nextflow

// activate dsl2
nextflow.enable.dsl = 2

//
params.help = false
def helpMessage() {
    log.info """\
             ================================================
             ================================================
             ~ version ${workflow.manifest.version}

             Usage:
             nextflow run main.nf -profile docker --input "/path/to/sampleList.csv"

             Input arguments:
                --input             Specify the directory of fastq.gz file(s) in CSV file format.
                                    CSV headers for Illumina reads: "sample_id,short_read_1,short_read_2",
                                    CSV headers for Nanopore reads: "sample_id,long_read".

             Output arguments:
                --outdir            Specify the output directory. [Default: "./results/"]
                --mode              Specify the analysis mode: mapping, assembly, nanopore. [Default: "mapping"]

             Database arguments:
                --metaphlan4_db     Specify the location of unarchived MetaPhlAn4 database.
                --kraken2_db        Specify the location of unarchived Kraken2 database.
             """.stripIndent()
}
if (params.help) {
    helpMessage()
    exit 0
}

// print version and exit
params.version = false
if (params.version) {
    println """\
            =================================================
            G B S   M E T A G E N O M I C S   P I P E L I N E
            =================================================
            ~ version ${workflow.manifest.version}
            """.stripIndent()
    ["bash", "${projectDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// pipeline must-have parameters
if (params.input)             { ch_input = file(params.input) }                 else { exit 1, 'Input spreadsheat is not specified.' }
if (params.reference)         { ch_ref = file(params.reference) }               else { exit 1, 'Reference genome is not specified.' }
if (params.seqtype_reference) { ch_st_ref = file(params.seqtype_reference) }    else { exit 1, 'Sequence type genome is not specified.'}
if (params.metaphlan4_db)     { ch_metaphlan4_db = file(params.metaphlan4_db) } else { exit 1, 'Metaphlan4 database is not specified.' }
if (params.kraken2_db)        { ch_kraken2_db = file(params.kraken2_db) }       else { exit 1, 'Kraken2 database is not specified.' }
if (params.gtdb)              { ch_gtdb = Channel.value(file(params.gtdb)) }    else { exit 1, 'Unarchived GTDB-Tk reference data is not specified' }

/*
 * Get number of columns from sampleList.csv: sample_id, short_read_1, short_read_2
 * Get number of samples from sampleList.csv
 */
def getColNo(fileName) {
    lines = file(fileName).readLines()
    return lines[0].split('\t').size()
}

sampleSize = file(params.input).countLines() - 1

/*
 * PROCESSES
 */
include { INPUT_CHECK; INPUT_CHECK_NANOPORE                                         } from './modules/qc_input'
include { FASTP                                                                     } from './modules/trim_input'
include { BWA_INDEX as BWA_INDEX_GBS; BWA_INDEX as BWA_INDEX_ST                     } from './modules/map_reads'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GBS; SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_ST } from './modules/map_reads'
include { BWA_MEM; BWA_SNP_ALIGN; BCFTOOLS_MPILEUP; MINIMAP2_NANOPORE               } from './modules/map_reads'
include { SPADES; METAWRAP; FLYE; MEDAKA                                            } from './modules/assembly'
include { QUAST; QUAST_NANOPORE                                                     } from './modules/qc_assembly'
include { MLST                                                                      } from './modules/typing'
include { METAPHLAN4; KRAKEN2; GTDB_PREP; GTDBTK                                    } from './modules/assign_taxo'

/*
 * MAIN WORKFLOW
 */
workflow {

    /*
    ======================================
    MAPPING-BASED METAGENOMICS
    ======================================
    */

    GTDB_PREP( ch_gtdb )

    if ( params.mode == 'mapping' ) {

        // DEFINE INPUT: params.input (a spreadsheet or directory of paired-end fastq.gz files)
        INPUT_CHECK()
        ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
        // QC READS
        FASTP( ch_raw_short_reads )
        // MAPPING TO GBS REF
        BWA_INDEX_GBS( ch_ref )
        SAMTOOLS_FAIDX_GBS( ch_ref )
        BWA_MEM( FASTP.out.reads, BWA_INDEX_GBS.out.bwa_index, SAMTOOLS_FAIDX_GBS.out.samtools_faidx )
        // DE NOVO ASSEMBLY
        SPADES( BWA_MEM.out.reads )
        // QC ASSEMBLY
        QUAST( BWA_MEM.out.reads, SPADES.out.assembly )
        // GBS SEQUENCE TYPING
        MLST( SPADES.out.assembly )
        // TAXONOMICAL ASSIGNMENT
        METAPHLAN4( FASTP.out.reads, ch_metaphlan4_db )
        GTDBTK( SPADES.out.assembly, GTDB_PREP.out )
        // SNP PROFILING
        BWA_INDEX_ST( ch_st_ref )
        SAMTOOLS_FAIDX_ST( ch_st_ref )
        BWA_SNP_ALIGN( BWA_MEM.out.reads, BWA_INDEX_ST.out.bwa_index, SAMTOOLS_FAIDX_ST.out.samtools_faidx )
        BCFTOOLS_MPILEUP( BWA_SNP_ALIGN.out.bam, SAMTOOLS_FAIDX_ST.out.samtools_faidx )

    } else if ( params.mode == "assembly" ) {

        // DEFINE INPUT: params.input (a spreadsheet or directory of paired-end fastq.gz files)
        INPUT_CHECK()
        ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
        // QC READS
        FASTP( ch_raw_short_reads )
        // METAGENOME-BASED ASSEMBLY
        METAWRAP( FASTP.out.reads )
        // TAXONOMICAL ASSIGNMENT
        GTDBTK( METAWRAP.out.mag, GTDB_PREP.out )
        // GBS SEQUENCE TYPING
        MLST( METAWRAP.out.mag )

    } else if ( params.mode == "nanopore" ) {

        // DEFINE INPUT: params.input (a spreadsheet or directory of nanopore fastq files)
        INPUT_CHECK_NANOPORE()
        ch_raw_long_reads = INPUT_CHECK_NANOPORE.out.raw_long_reads
        // MAPPING TO GBS REF
        MINIMAP2_NANOPORE( ch_raw_long_reads, ch_ref )
        // GENOME ASSEMBLY & POLISHING
        FLYE( MINIMAP2_NANOPORE.out.reads )
        MEDAKA( MINIMAP2_NANOPORE.out.reads, FLYE.out.assembly )
        // QC ASSEMBLY
        QUAST_NANOPORE( MINIMAP2_NANOPORE.out.reads, MEDAKA.out.assembly )
        // TAXONOMICAL ASSIGNMENT
        GTDBTK( MEDAKA.out.assembly, GTDB_PREP.out)
        // GBS SEQUENCE TYPING
        MLST( MEDAKA.out.assembly )
    }
}



/*
 * WORKFLOW TRACING
 */
// With errors
workflow.onError {
    log.info """\
             Opps... The pipeline execution stopped with the following message: ${workflow.errorMessage}
             """
             .stripIndent()
}

// In general
workflow.onComplete {
    log.info """\
             --------------------------
             Pipeline Execution Summary
             --------------------------
             Name           : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}
             Profile        : ${workflow.profile}
             Launch Dir     : ${workflow.launchDir}
             Work Dir       : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : ""}
             Status         : ${workflow.success ? "success" : "failed"}
             Error report   : ${workflow.errorReport ?: "-"}
             """.stripIndent()

             // run a small clean-up script to remove "work/" directory after successful completion
             if (!params.debug && workflow.success){
                 ["bash", "${projectDir}/bin/clean.sh", "${workflow.runName}"].execute()
             }
}
