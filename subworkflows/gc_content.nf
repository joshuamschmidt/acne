/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PFBTOBED     } from '../modules/local/pfb_to_bed'
include { BEDTOOLS_NUC } from '../modules/local/bedtools/nuc'
include { EXCTRACTGC   } from '../modules/local/extract_gc'

workflow GC_CONTENT {
  take:
    ch_pfb // channel: [ val(meta), path(pfb) ]
    ch_fasta // channel: [ path(fasta) ]
    ch_fai   // channel: [ path(fai) ]

  main:
    ch_versions = Channel.empty()

    PFBTOBED(ch_pfb)
    BEDTOOLS_NUC(PFBTOBED.out.bed, ch_fasta, ch_fai)
    EXCTRACTGC(BEDTOOLS_NUC.out.output)
  
  emit:
  gc      = EXCTRACTGC.out.gc     // channel: [ val(meta), path(gc) ]
}
