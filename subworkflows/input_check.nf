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
        .map { create_gsfile_channel(it) }
        .set { rawfiles }

    TABIX_BGZIP(rawfiles)

    emit:
    gsfiles =  TABIX_BGZIP.out.output  // channel: [ val(meta), [ gsfile ] ]
    // TODO versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ gsfile, etc ] ]
// following example from https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/input_check.nf


def create_gsfile_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id  = row.id

    // add path of the gsfile to the meta map
    def gsfile_meta = []
    if (!file(row.gsfile).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Genome Studio file does not exist!\n${row.gsfile}"
    } else {
        gsfile_meta = tuple(meta, file(row.gsfile))
    }
    return gsfile_meta
}

def     ( ) {
    // create meta map
    def meta = [:]
    meta.id  = row.id


}
