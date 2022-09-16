#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// default params
params.outdir = 'run/'
params.inputFile = 'inputFile.txt'
params.runID = ''
params.gc_model = ''
params.hmm = ''
params.split = false

log.info """\
 A C N E - N F   P I P E L I N E
 ===================================
 input        : ${params.inputFile}
 runID        : ${params.runID}
 gcModel      : ${params.gc_model}
 hmm          : ${params.hmm}
 split        : ${params.split}
 outdir       : ${params.outdir}
 """



Channel
    .fromPath(params.inputFile)
    .set { input_ch }





