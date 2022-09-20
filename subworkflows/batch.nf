/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKEPFB        } from '../modules/local/pfb/main'
include { SPLITGS        } from '../modules/local/split/main'
include { PENNCNV_GC     } from '../modules/local/detect/main'


workflow BATCH_CALL {
  take:
    gs_file

    main:
    MAKEPFB(gs_file) | flatten | set {  pfb_ch }
    SPLITGS(gs_file)
    Channel
        .from( pfb_ch )
        .first()
        .set { gc_in_ch }
    PENNCNV_GC(gc_in_ch, Channel.fromPath(params.gc_model))
}
