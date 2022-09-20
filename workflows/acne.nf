// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.gc_model,
    params.hmm
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Set input
//ch_input_sample = path(params.input, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARTITIONGS } from '../modules/local/partition/main'
include { BATCH_CALL  } from '../subworkflows/batch'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACNE {
  take:
    gs_file
    split
    split_n

  main:
    if (split) {
        if (!split_n) {
        log.error "must inlude n samples to split large GS project"
        exit 1
        }
        PARTITIONGS(gs_file, split_n)
        //split_channel = PARTITIONGS.gs
        BATCH_CALL(PARTITIONGS.out.gs)
        }
    //BATCH_CALL(params.runID, ch_input_sample)
    }
