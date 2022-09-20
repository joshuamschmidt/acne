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
    MAKEPFB(gs_file)
    SPLITGS(gs_file)
    Channel
        .from( MAKEPFB.out )
        .first()
        .set { gc_in_ch }
    gc_in_ch.view()
    //PENNCNV_GC(gc_in_ch, Channel.fromPath(params.gc_model))
}
