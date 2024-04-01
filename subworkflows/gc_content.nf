/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PFBTOBED              } from '../modules/local/pfb_to_bed'
include { BEDTOOLS_MAKEWINDOWS  } from '../modules/nf-core/bedtools/split_gs'
include { BEDTOOLS_NUC          } from '../modules/local/bedtools/nuc'


workflow GC_CONTENT {
  take:
    gs_file

  main:
    MAKEPFB(gs_file)
    SPLITGS(gs_file)
    Channel
        .from( MAKEPFB.out )
        .first()
        .set { gc_in_ch }
    gc_in_ch.view()
    //PENNCNV_GC(gc_in_ch, Channel.fromPath(params.gc_model))
}
