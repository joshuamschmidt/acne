//
// Check input samplesheet and get gsfile channels
// TODO
// Check format of samplesheet. At the moment it just queues gsfile channel

// Check default hmm, gc_model file
// Check user supplied input(s).
// if gz, then decompress.....

//
// Check input samplesheet and get gsfile channels
//


include { SAMPLESHEET_CHECK } from '../modules/local/samplesheet_check'
include { TABIX_BGZIP  } from '../modules/nf-core/modules/tabix/bgzip/main'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_channel_from_input(it) }
        .set { rawfiles }

    TABIX_BGZIP(rawfiles)

    emit:
    gsfiles =  TABIX_BGZIP.out.output  // channel: [ val(meta), [ gsfile ] ]
    // TODO versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ gsfile, etc ] ]
// following example from https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/input_check.nf


def create_channel_from_input(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id  = row.id

    // add path of the gsfile to the meta map
    def gsfile_meta = []
    gsfile_meta << meta

    !file(row.gsfile).exists() ? (exit 1, "Genome Studio file does not exist!\n${row.gsfile}") : !row.gsfile.toString().endsWith(".gz") ? (exit 1, "Genome Studio file should have gz extension!\n${row.gsfile}") : gsfile_meta <<  file(row.gsfile)
    !file(row.gsfile_samples).exists() ? (exit 1, "GS samples file does not exist!\n${row.gsfile_samples}") : !row.gsfile_samples.toString().endsWith(".csv") ? (exit 1, "GS samples file should have csv extension!\n${row.gsfile}") : gsfile_meta <<  file(row.gsfile_samples)
    return gsfile_meta
}

def     ( ) {
    // create meta map
    def meta = [:]
    meta.id  = row.id


}
