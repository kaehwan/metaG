/*
 * Check file extension
 */
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

/* 
 * Check input samplesheet and get read channels
 */
workflow INPUT_CHECK {
    main:
    if(hasExtension(params.input, "csv")) {
        Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                def id  = row.sample_id
                def sr1 = row.short_read_1 ? file(row.short_read_1, checkIfExists: true) : false
                def sr2 = row.short_read_2 ? file(row.short_read_2, checkIfExists: true) : false
                return [ id, sr1, sr2 ]
            }
            .set { ch_input_rows }
        ch_input_rows
            .map { id, sr1, sr2 ->
                def meta = [:]
                meta.id = id
                return [ meta, [sr1, sr2] ]
            }
            .set { ch_raw_short_reads }
    } else {
        Channel
            .fromFilePairs(params.input, size: 2)
            .ifEmpty { exit 1, "Cannot find reads in: ${params.input}\n" }
            .map { row ->
                def meta = [:]
                meta.id = row[0]
                return [ meta, row[1] ]
            }
            .set { ch_raw_short_reads }
    }

    emit:
    raw_short_reads = ch_raw_short_reads
}

workflow INPUT_CHECK_NANOPORE {
    main:
    if(hasExtension(params.input, "csv")) {
        Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                def id  = row.sample_id
                def lr  = row.long_read ? file(row.long_read, checkIfExists: true) : false
                return [ id, lr ]
            }
            .set { ch_input_rows }
        ch_input_rows
            .map { id, lr ->
                def meta = [:]
                meta.id = id
                return [ meta, [lr] ]
            }
            .set { ch_raw_long_reads }
    } else {
        Channel
            .fromFilePairs(params.input, size: 1)
            .ifEmpty { exit 1, "Cannot find reads in: ${params.input}\n" }
            .map { row ->
                def meta = [:]
                meta.id = row[0]
                return [ meta, row[1] ]
            }
            .set { ch_raw_long_reads }
    }

    emit:
    raw_long_reads = ch_raw_long_reads
}