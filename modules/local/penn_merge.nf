process PENNCNV_MERGE {

    //container 'genomicslab/penncnv:latest'
    container 'wallen/penncnv:1.0.5'

    publishDir "$params.outdir/merge_calls", pattern: "*.mergedcnv", mode: 'copy'

    input:
    tuple val(meta), path(input), path(pfb)

    output:
    tuple val(meta), path ("${prefix}.mergedcnv"), emit: output
    
    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    clean_cnv.pl combineseg \
    --bp --fraction 0.50 \
    --signalfile $pfb \
    $input > ${prefix}.mergedcnv;
    """

}