//
// VALIDATE INPUTS
//

// Check default hmm, gc_model file
// Check user supplied input(s).
// if gz, then decompress.....

/*
example include:
include { GS_INPUT as BWAMEM1_INDEX             } from '../../modules/nf-core/modules/bwa/index/main'
*/



process GS_INPUT {

  input:
}



workflow VALIDATE_INPUTS {
    take:
        gs_file                 // channel: [mandatory] input gs_file
        hmm                     // channel: [optional]  penn cnv hmm for cnv calling
        gc_file                 // channel: [optional]  gc_file for GC correction of input data

    main:

    ch_versions = Channel.empty()


  }

