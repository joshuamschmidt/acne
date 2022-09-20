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
    MAKEPFB(gs_file) | set {  pfb_ch }
    SPLITGS(gs_file)
    PENNCNV_GC(pfb_ch, params.gc_model)
}
