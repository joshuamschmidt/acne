/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKEPFB        } from '../modules/local/pfb/main'
//include { SPLITGS        } from '../modules/local/split/main'
//include { PENNCNV_DETECT } from '../modules/local/detect/main'


workflow BATCH_CALL {
  take:
    gs_file

    main:
    MAKEPFB(gs_file)
}
