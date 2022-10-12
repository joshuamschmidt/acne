process CONCATENATE_PENN_CALLS {
    tag "$meta.id"
    //container 'genomicslab/penncnv:latest'
    //container 'wallen/penncnv:1.0.5'

    publishDir "$params.outdir/concatenated_calls", pattern: "*.catcnv", mode: 'copy'

    input:
    tuple val(meta), path ('*.rawcnv')

    output:
    tuple val(meta), path ("${prefix}.catcnv"), emit: output
    
    script:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    cat *.rawcnv > "$prefix".catcnv
    """
}