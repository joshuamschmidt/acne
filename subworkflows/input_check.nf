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
            .map{ create_meta_list_from_input(it) }
            .multiMap{  it -> 
                gsfile:          it[0]
                gs_sample_sheet: it[1]
                sample_include:  it[2]
                snp_include:     it[3]
            }
            .set { input_ch }

    TABIX_BGZIP(input_ch.gsfile)

    emit:
    gsfiles =  TABIX_BGZIP.out.output          // channel: [ val(meta), [ gsfile ] ]
    gs_sample_sheet = input_ch.gs_sample_sheet // channel: [ val(meta), [ gs_sample_sheet ] ] 
    sample_include = input_ch.sample_include   // channel: [ val(meta), [ sample_include ] ]
    snp_include = input_ch.snp_include         // channel: [ val(meta), [ snp_include ] ]
    // TODO versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function turn each row in the samplesheet into a list of tuples (meta, path) 


def create_meta_list_from_input(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id  = row.id

    def gsfile          = !file(row.gsfile).exists() ? (exit 1, "Genome Studio file does not exist!\n${row.gsfile}") : !row.gsfile.toString().endsWith(".gz") ? (exit 1, "Genome Studio file should have gz extension!\n${row.gsfile}") : tuple(meta, file(row.gsfile))
    def gs_sample_sheet = !file(row.gs_sample_sheet).exists() ? (exit 1, "GS samples file does not exist!\n${row.gs_sample_sheet}") : !row.gs_sample_sheet.toString().endsWith(".csv") ? (exit 1, "GS samples file should have csv extension!\n${row.gs_sample_sheet}") : tuple(meta, file(row.gs_sample_sheet))
    def sample_include  = !row.sample_include ? tuple(meta, '') : !file(row.sample_include).exists() ? (exit 1, "sample include file does not exist!\n${row.sample_include}") : tuple(meta, file(row.sample_include)) 
    def snp_include     = !row.snp_include ? tuple(meta, '') : !file(row.snp_include).exists() ? (exit 1, "SNP include file does not exist!\n${row.snp_include}") : tuple(meta, file(row.snp_include))

    def out = [gsfile, gs_sample_sheet, sample_include, snp_include]

    return out
}


