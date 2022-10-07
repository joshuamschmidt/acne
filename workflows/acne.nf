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
include { SPLITGS      } from '../modules/local/split_gs'
include { MAKEPFB      } from '../modules/local/make_pfb'
include { PENNCNV_GC   } from '../modules/local/penn_gc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACNE {

    INPUT_CHECK( ch_input )
        //.view()
        //.set{ ch_g_input }


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
        PARTITIONGS(INPUT_CHECK.out.gsfiles, params.partition_n)
            .transpose()
            .map {
                meta, partition ->
                def meta_clone = meta.clone()
                parition_suffix = partition.baseName.split('-').last()
                meta_clone.id=meta_clone.id+'_'+parition_suffix
                [ meta_clone, partition ]
            }
            .set { ch_pre_split }

        // pfb and gcmodel files, unique to each subbatch instantiated by PARTITIONGS
        MAKEPFB( ch_pre_split )

        PENNCNV_GC( MAKEPFB.out.output, params.gc_model )
            .out.output.cross(MAKEPFB.out.output)
            .view()


        //
        SPLITGS( ch_pre_split )
            .transpose()
            .set { ch_post_split }


        //ch_post_split.cross(MAKEPFB.out.output).view()

    } else {
        SPLITGS( INPUT_CHECK.out.gsfiles )

    }
    //BATCH_CALL(params.runID, ch_input_sample)
}
