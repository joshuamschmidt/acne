process PENNCNV_DETECT {

    //container 'genomicslab/penncnv:1.0.5'
    container 'wallen/penncnv:1.0.5'

    publishDir "$params.outdir/penn_call_logs", pattern: "*.rawcnv.log", mode: 'copy'
    publishDir "$params.outdir/raw_calls", pattern: "*.rawcnv", mode: 'copy'

    input:
    tuple val(meta), path(input), path(gc_model), path(pfb)
    path(hmm)

    output:
    tuple val(meta), path ("${input.baseName}.rawcnv"),         emit: raw_call
    tuple val(meta), path ("${input.baseName}.rawcnv.log"),     emit: raw_log
    
    script:

    """
    detect_cnv.pl \
    -test \
    --loh \
    --confidence \
    --hmmfile $hmm \
    --pfbfile $pfb \
    $input \
    --logfile ${input.baseName}.rawcnv.log \
    --gcmodelfile $gc_model \
    --output ${input.baseName}.rawcnv;
    """

}
