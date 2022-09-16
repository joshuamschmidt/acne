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

process partition {

    label 'pythonTasks'
}


process split {

    label 'pythonTasks'

    input:
    file input_file from input_ch

    output:
    path '*.input.txt' into split_input_ch

    script:
    """
    splitGS.py --input $input_file
    """
}

cnv_call_ch = split.input_ch.flatten()
    .map {
        tuple( it.name.split('.')[0], it )
    }

process makePFB {
    publishDir "$params.outdir/PFB/", pattern: "*.pfb", mode: 'copy'

    label 'pythonTasks'

    input:
    file input_file from input_ch

    output:
    file "runID.pfb" into pfb_ch

    script:
    """
    makePFB.py --input $input_file \
    --output ".regions.bed.gz"
    """
}

process callCNV {
    publishDir "$params.outdir/logs", pattern: "*.rawcnv.log", mode: 'copy'
    publishDir "$params.outdir/raw_calls", pattern: "*.rawcnv", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file) from cnv_call_ch
    file pfb_file from pfb_ch

    output:
    tuple val(sample_id), path("${sample_id}.rawcnv") into raw_call_ch
    path "${sample_id}.rawcnv.log" into call_log_ch
    
    script:

    """
    "$penncnv"/detect_cnv.pl \
    -test \
    --loh \
    --confidence \
    --hmmfile $hmm \
    --pfbfile $pfb_file \
    --logfile "$sample_id".rawcnv.log 
    --gcmodelfile $gc_model \
    --output  "$sample_id".rawcnv \
    $input_file;
    """

}
