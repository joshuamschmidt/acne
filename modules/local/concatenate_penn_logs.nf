process CONCATENATE_PENN_LOGS {

    //container 'genomicslab/penncnv:latest'
    //container 'wallen/penncnv:1.0.5'

    publishDir "$params.outdir/concatenated_calls", pattern: "*.catcnv.log", mode: 'copy'

    input:
    tuple val(meta), path('*.rawcnv.log')

    output:
    tuple val(meta), path ("${prefix}.catcnv.log"), emit: output
    
    script:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    cat *.rawcnv.log | sed -e '/detect_cnv/d' > "$prefix".catcnv.log
    """
}