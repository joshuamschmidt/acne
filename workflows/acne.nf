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
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK  } from '../subworkflows/input_check'
include { BATCH_CALL   } from '../subworkflows/batch'
include { PARTITIONGS  } from '../modules/local/partition'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACNE {

    INPUT_CHECK( ch_input )
        .view()
        .set{ ch_g_input }


    // big GS files can be partitioned for efficiency
    // 1: partition into batches of size partition_n. (affects PBF and GC model creation steps)
    // 2: also partition how many ind files are passed to the CNV caller.

    // note also that splitting uses python polars - so more efficient using uncompressed input
    if (params.partition) {
        if (!params.partition_n) {
        log.error "must inlude n samples to split large GS project"
        exit 1
        }
        // PARTITIONGS returns uncompressed
        PARTITIONGS(INPUT_CHECK.out.output, params.partition_n)
        //| flatten | BATCH_CALL
        //Channel
        //    .from(PARTITIONGS.out)
        //    .set { split_channel }
        //BATCH_CALL(split_channel)
        }
    //BATCH_CALL(params.runID, ch_input_sample)
    }
